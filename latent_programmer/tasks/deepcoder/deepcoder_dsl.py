# coding=utf-8
# Copyright 2022 The Google Research Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

r"""DeepCoder DSL.

Features:
  * The operators and constant lambda functions are as described in Appendix F
    of the DeepCoder paper (https://arxiv.org/pdf/1611.01989.pdf)
  * The Program class stores a program and can do the following:
      * Run the program on inputs to produce an output or program state
      * Serialize the program into a sequence of tokens
      * Deserialize a sequence of tokens into a new Program
  * The Statement class represents partial programs (a single assignment)
  * The ProgramState class captures program state at partial executions
  * Utilities for random dataset generation:
      * Generate a new random program
          * Check to make sure the random program doesn't have dead code
      * Generate random inputs for a program
      * Ways of generating compositional generalization splits
  * A flag controlling whether we use mod-10 arithmetic or not

Helpful properties to know:

  * When converting from a list of string tokens to strings (representing
    programs or program states) or vice versa, tokens are all space separated.
    There is no need to parse parentheses. Use `' '.join(tokens)` or
    `string.split()` to convert between lists of string tokens and strings.

  * A program state is given as a string in the following form, which does not
    distinguish between inputs and previously-computed local variables:

    x0 = [ 3 1 2 ] | x1 = [ 4 2 3 ] | x2 = 9

  * A program is given as a string in the following form:

    x0 = INPUT | x1 = Map (+1) x0 | x2 = Sum x1

  * In the context of running an entire program (e.g., to see if it satisfies
    the I/O examples), the program "output" is the result of its last line. In
    the context of running just one program line, we can think of the input and
    output as being entire program states.

  * Program variables (x0, x1, ...) only have types `int` or `list[int]`. Lambda
    functions may have types such as `int -> bool` or `(int, int) -> int`, but
    lambdas are only used as arguments to higher-order functions and can't be
    program variables themselves.

  * If desired, variables can be used in an order different from x0, x1, etc.
"""

import ast
import functools
import re
from typing import Any, Callable, Dict, List, Optional, Tuple, Type, Union

from absl import flags

_DEEPCODER_MOD = flags.DEFINE_integer(
    'deepcoder_mod', 0,
    'The modulo we use for DeepCoder arithmetic, or 0 to not apply any mod.')

_DEEPCODER_MAX_LIST_LENGTH = flags.DEFINE_integer(
    'deepcoder_max_list_length', 20,
    'The maximum length of a DeepCoder list input.')


def deepcoder_mod():
  """Hides the internal flag from code outside this module."""
  return _DEEPCODER_MOD.value


def deepcoder_max_list_length():
  """Hides the internal flag from code outside this module."""
  return _DEEPCODER_MAX_LIST_LENGTH.value


# Types for type specifications, e.g., `([int], bool)` has type LambdaType and
# specifies a lambda that takes 1 int and returns a bool.
InputsType = List[Union[Type[Any], 'LambdaType']]
OutputType = Type[Any]
LambdaType = Tuple[InputsType, OutputType]

# A type for the result of applying any operation.
ResultType = Union[int, List[int]]

MIN_INT = -256
MAX_INT = 255


def variable_token(index):
  return f'x{index}'

MAX_NUM_VARIABLES = 10  # How many variables are available to use.
ALL_VARIABLES = frozenset(variable_token(index)
                          for index in range(MAX_NUM_VARIABLES))


def variable_index_from_token(token):
  if not re.fullmatch(r'x\d+', token):
    raise ParseError(f'Invalid variable token: {token}')
  index = int(token[1:])
  if index < 0 or index >= MAX_NUM_VARIABLES:
    raise ParseError(f'Variable token has out-of-bounds index: {token}')
  return index


class ParseError(Exception):
  """Could not parse from a string or tokens, similar to a syntax error."""
  # This could happen if a model predicts a bad token sequence.


class RunError(Exception):
  """Could not execute the operation or program, similar to a runtime error."""
  # This is raised for type mismatch, wrong arity given to an operation,
  # reference to an invalid variable, or running a program on the wrong number
  # of inputs.
  # This is NOT raised when an operation returns None, which results in run()
  # returning None.


class DeepCoderError(Exception):
  """Something happened that suggests there's a bug."""
  # This is more severe of an error than a program failing to execute. We should
  # not see this error during normal usage even considering model inaccuracy,
  # unless there's a bug.


class Example(object):
  """A DeepCoder specification in the form of an I/O example."""

  def __init__(self, inputs, output):
    self.inputs = inputs
    self.output = output


class ProgramTask(object):
  """A DeepCoder program with I/O examples."""

  def __init__(self, program, examples):
    self.program = program
    self.examples = examples


def join_token_lists(token_lists,
                     separator_token):
  return functools.reduce(lambda a, b: a + [separator_token] + b, token_lists)


def mod_result(result):
  """If desired, apply mod to ints."""
  if result is None:
    return result
  if deepcoder_mod() > 0:
    result_type = type(result)
    if result_type == int:
      return result % deepcoder_mod()
    elif result_type == list:
      return [x % deepcoder_mod() for x in result]
    else:
      raise DeepCoderError(f'Unhandled result in mod_result: {result}')
  return result


def validate_int(i):
  """Checks that the integer is in range."""
  if deepcoder_mod() > 0:
    return 0 <= i < deepcoder_mod()
  return MIN_INT <= i <= MAX_INT


def validate_result(result):
  """Returns whether an object is valid or not."""
  # Distinguish between bool and int.
  # pylint: disable=unidiomatic-typecheck
  if type(result) == int:
    return validate_int(result)
  elif type(result) == list:
    return all(type(x) == int and validate_int(x) for x in result)
  else:
    return False
  # pylint: enable=unidiomatic-typecheck


def tokenize_result(result):
  """Returns a list of tokens for the result of an operation."""
  if isinstance(result, int):
    return [str(result)]
  elif isinstance(result, list):
    return ['['] + sum([tokenize_result(x) for x in result], []) + [']']
  else:
    raise DeepCoderError(f'Unhandled type in tokenize_result({result})')


def result_to_str(result):
  return ' '.join(tokenize_result(result))


class ProgramState(object):
  """Holds a program state (for one example)."""

  def __init__(self, state, variables):
    self.state = list(state)
    self.variables = list(variables)
    if len(state) != len(variables):
      raise DeepCoderError(
          f'`state` has length {len(state)} but `variables` has length '
          f'{len(variables)}')
    self._variable_to_result = {token: result
                                for token, result in zip(variables, state)}

  def __len__(self):
    return len(self.state)

  def get_index(self, index):
    if index < 0 or index >= len(self.state):
      raise RunError(f'Invalid index: {index}')
    return self.state[index]

  def get_variable(self, token):
    if token not in self._variable_to_result:
      raise RunError(f'Invalid variable token: {token}')
    return self._variable_to_result[token]

  def __eq__(self, other):
    return (isinstance(other, ProgramState) and
            self.state == other.state and
            self.variables == other.variables)

  def copy(self):
    return ProgramState(list(self.state), list(self.variables))

  def add_result(self, result, variable):
    if variable in self._variable_to_result:
      raise RunError(f'Cannot add new result for variable {variable} that '
                     f'already exists in state: {self}')
    self.state.append(result)
    self.variables.append(variable)
    self._variable_to_result[variable] = result

  def get_output(self):
    return self.state[-1]

  def get_output_variable(self):
    return self.variables[-1]

  def tokenize(self):
    lines = [[token, '='] + tokenize_result(result)
             for token, result in zip(self.variables, self.state)]
    return join_token_lists(lines, separator_token='|')

  def __str__(self):
    return ' '.join(self.tokenize())

  @classmethod
  def from_tokens(cls, tokens):
    return cls.from_str(' '.join(tokens))

  @classmethod
  def from_str(cls, string):
    """Creates a ProgramState from its string representation."""
    lines = [line.strip() for line in string.split('|')]
    state = []
    variables = []
    for line in lines:
      splitted = line.split('=')
      if len(splitted) != 2:
        raise ParseError(f"Expected exactly one '=': {line}")
      lhs, rhs = [part.strip() for part in splitted]
      _ = variable_index_from_token(lhs)  # Make sure lhs has the right format.
      variable = lhs
      if variable in variables:
        raise ParseError(f'Found duplicate variable: {variable}')
      # If we find spaces surrounded by digits, place a comma there. This
      # converts comma-less "[ 6 7 8 ]" into something parseable, "[ 6, 7, 8 ]".
      rhs = re.sub(r'(?<=\d) +(?=-|\d)', ', ', rhs)
      try:
        result = ast.literal_eval(rhs)
        result_is_valid = validate_result(result)
      except (ValueError, SyntaxError):
        result = None
        result_is_valid = False
      if not result_is_valid:
        raise ParseError(f'Found invalid result `{result}` from RHS {rhs!r}')
      state.append(result)
      variables.append(variable)
    return cls(state, variables)


class Function(object):
  """Base class for functionality in the DeepCoder DSL."""

  def __init__(self, token, func,
               inputs_type, output_type):
    self.token = token
    self.func = func
    self.inputs_type = inputs_type
    self.output_type = output_type
    self.arity = len(inputs_type)


class Lambda(Function):
  """A lambda function like `*2` or `+`."""


class Operation(Function):
  """Base class for first-order and higher-order operations."""

  def run(self, inputs):
    """Runs an operation on input arguments."""
    if len(inputs) != self.arity:
      raise RunError(f'Arity was {len(inputs)} but expected {self.arity}')
    try:
      result = self.func(*inputs)
    except TypeError as e:
      raise RunError(
          f'Encountered TypeError when running {self.func} on {inputs}') from e
    result = mod_result(result)
    if not validate_result(result):
      return None
    return result


class FirstOrderOperation(Operation):
  """A first-order function like `Head` or `Reverse`."""


class HigherOrderOperation(Operation):
  """A higher-order function like `Map` or `ZipWith`."""


class Statement(object):
  """One statement in a program."""

  def __init__(self, variable, operation,
               args):
    self.variable = variable
    self.operation = operation
    self.args = list(args)

  def run(self, initial_state):
    """Runs the operation and assigns it to a variable."""
    if self.variable in initial_state.variables:
      raise RunError(
          f'Statement has variable {self.variable} which already exists in the '
          f'initial state: {initial_state}')
    arg_values = []
    for arg in self.args:
      if isinstance(arg, str):
        arg_values.append(initial_state.get_variable(arg))
      elif isinstance(arg, Lambda):
        arg_values.append(arg.func)
      else:
        raise DeepCoderError(
            f'Unhandled argument {arg} for statement {self} and initial_state '
            f'{initial_state}')
    result = self.operation.run(arg_values)
    if result is None:
      return None
    result_state = initial_state.copy()
    result_state.add_result(result, self.variable)
    return result_state

  def tokenize(self):
    tokens = [self.variable, '=', self.operation.token]
    for arg in self.args:
      if isinstance(arg, str):
        tokens.append(arg)
      elif isinstance(arg, Lambda):
        tokens.append(arg.token)
      else:
        raise DeepCoderError(f'Unhandled argument: {arg}')
    return tokens

  def __str__(self):
    return ' '.join(self.tokenize())

  @classmethod
  def from_tokens(cls, tokens,
                  check_variable_name = True):
    """Parses a Statement from a list of tokens."""
    # Parse LHS variable, =, and the operation.
    if len(tokens) < 4:
      raise ParseError(f'Too few tokens: {tokens}')
    variable = tokens[0]
    if check_variable_name:
      _ = variable_index_from_token(variable)  # Check its format.
    if tokens[1] != '=':
      raise ParseError(f"Second token must be '=': {tokens}")
    operation_token = tokens[2]
    if operation_token not in TOKEN_TO_OPERATION:
      raise ParseError(f'Unknown operation {operation_token} in: {tokens}')
    operation = TOKEN_TO_OPERATION[operation_token]

    # Parse operation arguments.
    args = []
    for i, token in enumerate(tokens[3:]):
      if token in TOKEN_TO_LAMBDA:
        if isinstance(operation, FirstOrderOperation) or i > 0:
          raise ParseError(f'Did not expect lambda at token {i}: {tokens}')
        arg = TOKEN_TO_LAMBDA[token]
      else:
        if isinstance(operation, HigherOrderOperation) and i == 0:
          raise ParseError(f'Expected lambda for token {i}: {tokens}')
        arg = token
        if arg == variable:
          raise ParseError(f'Cannot use LHS variable as arg: {tokens}')
        _ = variable_index_from_token(arg)  # Check its format.
      args.append(arg)
    if len(args) != operation.arity:
      raise ParseError(
          f'Statement tokens have wrong arity for operation: {tokens}')
    return cls(variable, operation, args)

  @classmethod
  def from_str(cls, string,
               check_variable_name = True):
    return cls.from_tokens(string.split(' '), check_variable_name)


class Program(object):
  """A full DeepCoder program including input handling."""

  def __init__(self, input_variables, statements):
    self.input_variables = list(input_variables)
    self.statements = list(statements)
    self.num_inputs = len(input_variables)

    if (len(set(self.get_variables()))
        != len(statements) + len(input_variables)):
      raise RunError(
          f'A variable token is duplicated in the program with '
          f'input_variables={input_variables} and statements={statements}.')

  def get_variables(self):
    return self.input_variables + [s.variable for s in self.statements]

  def run(self, inputs):
    if len(inputs) != self.num_inputs:
      raise RunError(f'Got {len(inputs)} inputs but expected {self.num_inputs}')
    state = ProgramState(inputs, self.input_variables)
    for statement in self.statements:
      state = statement.run(state)
      if state is None:
        return None
    return state

  def __len__(self):
    return len(self.statements)

  def tokenize(self):
    lines = []
    for input_variable in self.input_variables:
      lines.append([input_variable, '=', 'INPUT'])
    for statement in self.statements:
      lines.append(statement.tokenize())
    return join_token_lists(lines, separator_token='|')

  def __str__(self):
    return ' '.join(self.tokenize())

  @classmethod
  def from_tokens(cls, tokens):
    """Parses a Program from a list of tokens."""
    # Split tokens into lines at '|'.
    lines = []
    current_line = []
    for token in tokens:
      if token == '|':
        lines.append(current_line)
        current_line = []
      else:
        current_line.append(token)
    if current_line:
      lines.append(current_line)

    # Parse statements.
    found_non_input = False
    input_variables = []
    statements = []
    for line in lines:
      if line[-1] == 'INPUT':
        if found_non_input:
          raise ParseError(f'Found INPUT after a statement: {lines}')
        variable = line[0]
        _ = variable_index_from_token(variable)  # Check its format.
        if line != [variable, '=', 'INPUT']:
          raise ParseError(f'Error while parsing INPUT statement: {line}')
        input_variables.append(variable)
      else:
        found_non_input = True
        statements.append(Statement.from_tokens(line))
    return cls(input_variables, statements)

  @classmethod
  def from_str(cls, string):
    return cls.from_tokens(string.split(' '))


def _scanl1(f, xs):
  ys = []
  for i, x in enumerate(xs):
    if i == 0:
      ys.append(x)
    else:
      ys.append(f(ys[-1], x))
  return ys

# Use the Python code from Appendix F in the DeepCoder paper.
# pylint: disable=g-explicit-length-test, unnecessary-lambda
LAMBDAS = [
    Lambda('(+1)', lambda x: x + 1, [int], int),
    Lambda('(-1)', lambda x: x - 1, [int], int),
    Lambda('(*2)', lambda x: x * 2, [int], int),
    Lambda('(/2)', lambda x: x // 2, [int], int),
    Lambda('(*(-1))', lambda x: -x, [int], int),
    Lambda('(**2)', lambda x: x ** 2, [int], int),
    Lambda('(*3)', lambda x: x * 3, [int], int),
    Lambda('(/3)', lambda x: x // 3, [int], int),
    Lambda('(*4)', lambda x: x * 4, [int], int),
    Lambda('(/4)', lambda x: x // 4, [int], int),
    Lambda('(>0)', lambda x: x > 0, [int], bool),
    Lambda('(<0)', lambda x: x < 0, [int], bool),
    Lambda('(%2==0)', lambda x: x % 2 == 0, [int], bool),
    Lambda('(%2==1)', lambda x: x % 2 == 1, [int], bool),
    Lambda('(+)', lambda x, y: x + y, [int, int], int),
    Lambda('(-)', lambda x, y: x - y, [int, int], int),
    Lambda('(*)', lambda x, y: x * y, [int, int], int),
    Lambda('(min)', lambda x, y: min(x, y), [int, int], int),
    Lambda('(max)', lambda x, y: max(x, y), [int, int], int),
]

FIRST_ORDER_OPERATIONS = [
    FirstOrderOperation(
        'Head', lambda xs: xs[0] if len(xs) > 0 else None, [list], int),
    FirstOrderOperation(
        'Last', lambda xs: xs[-1] if len(xs) > 0 else None, [list], int),
    FirstOrderOperation(
        'Take', lambda n, xs: xs[:n], [int, list], list),
    FirstOrderOperation(
        'Drop', lambda n, xs: xs[n:], [int, list], list),
    FirstOrderOperation(
        'Access', lambda n, xs: xs[n] if 0 <= n < len(xs) else None,
        [int, list], int),
    FirstOrderOperation(
        'Minimum', lambda xs: min(xs) if len(xs) > 0 else None, [list], int),
    FirstOrderOperation(
        'Maximum', lambda xs: max(xs) if len(xs) > 0 else None, [list], int),
    FirstOrderOperation(
        'Reverse', lambda xs: list(reversed(xs)), [list], list),
    FirstOrderOperation(
        'Sort', lambda xs: sorted(xs), [list], list),
    FirstOrderOperation(
        'Sum', lambda xs: sum(xs), [list], int),
]
HIGHER_ORDER_OPERATIONS = [
    HigherOrderOperation(
        'Map', lambda f, xs: [f(x) for x in xs], [([int], int), list], list),
    HigherOrderOperation(
        'Filter', lambda f, xs: [x for x in xs if f(x)],
        [([int], bool), list], list),
    HigherOrderOperation(
        'Count', lambda f, xs: len([x for x in xs if f(x)]),
        [([int], bool), list], int),
    HigherOrderOperation(
        'ZipWith', lambda f, xs, ys: [f(x, y) for (x, y) in zip(xs, ys)],
        [([int, int], int), list, list], list),
    HigherOrderOperation(
        'Scanl1', _scanl1, [([int, int], int), list], list),
]
OPERATIONS = FIRST_ORDER_OPERATIONS + HIGHER_ORDER_OPERATIONS
# pylint: enable=g-explicit-length-test, unnecessary-lambda

# Maps from tokens to lambdas/operations.
TOKEN_TO_LAMBDA = {l.token: l for l in LAMBDAS}
TOKEN_TO_OPERATION = {op.token: op for op in OPERATIONS}

# Subsets of DSL functionality for generating compositional generalization
# datasets.
LAMBDAS_ONLY_MINUS_MIN = [TOKEN_TO_LAMBDA['(-)'], TOKEN_TO_LAMBDA['(min)']]
OPERATIONS_ONLY_SCAN = [TOKEN_TO_OPERATION['Scanl1']]
OPERATIONS_NO_SCAN = [op for op in OPERATIONS if op.token != 'Scanl1']


PAD, BOS, EOS, SEP = '', '<BOS>', '<EOS>', '|'
PAD_ID, BOS_ID, EOS_ID, SEP_ID = 0, 1, 2, 3


def vocab_tables():
  """Returns id-to-token and token-to-id vocabulary mappings."""
  # These tokens should be constant unless we change the DSL.
  tokens = [PAD, BOS, EOS, SEP]
  tokens.extend(['=', 'INPUT', '[', ']'])
  tokens.extend(op.token for op in OPERATIONS)
  tokens.extend(l.token for l in LAMBDAS)

  # These tokens may change based on flags or other settings.
  tokens.extend(variable_token(index) for index in range(MAX_NUM_VARIABLES))
  if deepcoder_mod() > 0:
    min_int = 0
    max_int = deepcoder_mod() - 1
  else:
    min_int, max_int = MIN_INT, MAX_INT
  tokens.extend(str(i) for i in range(min_int, max_int + 1))

  # Construct mappings between tokens and ids.
  id_to_token = {i: token for i, token in enumerate(tokens)}
  token_to_id = {token: i for i, token in enumerate(tokens)}
  assert len(id_to_token) == len(token_to_id)

  return id_to_token, token_to_id

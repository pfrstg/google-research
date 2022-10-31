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

"""Batch Ensemble MSG Builder."""

from typing import Iterator, List, Optional
import acme
from acme import adders
from acme import core
from acme import specs
from acme.agents.jax import actors
from acme.agents.jax import builders
from acme.jax import networks as networks_lib
from acme.jax import variable_utils
from acme.utils import counting
from acme.utils import loggers
import reverb
from jrl.agents.batch_ensemble_msg import config
from jrl.agents.batch_ensemble_msg import learning
from jrl.agents.batch_ensemble_msg import networks as batch_ensemble_msg_networks


class BatchEnsembleMSGBuilder(builders.ActorLearnerBuilder):
  """Efficient MSG Builder."""

  def __init__(
      self,
      config,
      # make_demonstrations: Callable[[int], Iterator[types.Transition]],
      make_demonstrations,
      ):
    self._config = config
    self._make_demonstrations = make_demonstrations

  def make_learner(
      self,
      random_key,
      networks,
      dataset,
      logger_fn,
      replay_client = None,
      counter = None,
      checkpoint = False,
  ):
    del dataset  # Offline RL

    data_iter = self._make_demonstrations()

    return learning.BatchEnsembleMSGLearner(
        ensemble_size=self._config.ensemble_size,
        beta=self._config.beta,
        networks=networks,
        rng=random_key,
        iterator=data_iter,
        num_bc_iters=self._config.num_bc_iters,
        target_entropy=self._config.target_entropy,
        behavior_regularization_type=self._config.behavior_regularization_type,
        behavior_regularization_alpha=self._config
        .behavior_regularization_alpha,
        policy_lr=self._config.policy_lr,
        q_lr=self._config.q_lr,
        counter=counter,
        logger=logger_fn('learner'),
        num_sgd_steps_per_step=self._config.num_sgd_steps_per_step,)

  def make_actor(
      self,
      random_key,
      policy_network,
      adder = None,
      variable_source = None):
    assert variable_source is not None
    return actors.GenericActor(
        actor=policy_network,
        random_key=random_key,
        # Inference happens on CPU, so it's better to move variables there too.
        variable_client=variable_utils.VariableClient(
            variable_source, 'policy', device='cpu'),
        adder=adder,
    )

  def make_replay_tables(
      self,
      environment_spec,
      policy
  ):
    """Create tables to insert data into."""
    del policy
    return []

  def make_dataset_iterator(
      self,
      replay_client):
    """Create a dataset iterator to use for learning/updating the agent."""
    return None

  def make_adder(self,
                 replay_client):
    """Create an adder which records data generated by the actor/environment."""
    return None
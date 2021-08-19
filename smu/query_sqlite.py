# coding=utf-8
# Copyright 2021 The Google Research Authors.
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

"""Queries the SMU sqlite database.

Command line interface to extract molecules from thq SMU database.
"""

import contextlib
import enum
import sys
from typing import Sequence

from absl import app
from absl import flags
from absl import logging
from rdkit import Chem

from tensorflow.io import gfile
from smu import smu_sqlite
from smu.parser import smu_utils_lib


class OutputFormat(enum.Enum):
  pbtxt = 1
  sdf_opt = 2
  sdf_init = 3
  sdf_init_opt = 4
  # Protocol buffer serialized string representation.
  proto = 5


flags.DEFINE_string(
    'input_sqlite', None,
    'Path of sqlite file to read. Must be on the local filesystem.')
flags.DEFINE_string(
    'output_path', None,
    'Path to output file to write. If not specified, will write to stdout.')
flags.DEFINE_list('btids', [], 'List of bond topology ids to query')
flags.DEFINE_list('cids', [], 'List of conformer ids to query')
flags.DEFINE_multi_string('smiles', [], 'List of smiles to query')
flags.DEFINE_enum_class('output_format', OutputFormat.pbtxt, OutputFormat,
                        'Format for the found SMU entries')

FLAGS = flags.FLAGS


class PBTextOutputter:
  """Simple internal class to write entries to text protocol buffer."""

  def __init__(self, output_path):
    """Creates PBTextOutputter.

    Args:
      output_path: file path to write to
    """
    if output_path:
      self.outfile = gfile.GFile(output_path, 'w')
    else:
      self.outfile = sys.stdout

  def write(self, conformer):
    """Writes a conformer.

    Args:
      conformer: dataset_pb2.Conformer
    """
    self.outfile.write(str(conformer))

  def close(self):
    self.outfile.close()

class SerializedProtoOutputter:
  """Simple internal class to write serialized protocol buffers to a file."""

  def __init__(self, output_path):
    """Creates SerializedProtoOutputter.

    Args:
      output_path: file path to write to
    """
    if output_path:
      self.outfile = gfile.GFile(output_path, 'w')
    else:
      self.outfile = sys.stdout   # This will not work.

  def write(self, conformer):
    """Writes a conformer.

    Args:
      conformer: dataset_pb2.Conformer
    """
    self.outfile.write(conformer.SerializeToString())

  def close(self):
    self.outfile.close()

class SDFOutputter:
  """Simple internal class to write entries as multi molecule SDF files."""

  def __init__(self, output_path, init_geometry, opt_geometry):
    """Creates SDFOutputter.

    At least one of init_geometry and opt_geometry should be True

    Args:
      output_path: file path to write to
      init_geometry: bool, whether to write with initial_geometries
      opt_geometry: bool, whether to write with optimized_geometry
    """
    self.init_geometry = init_geometry
    self.opt_geometry = opt_geometry
    if output_path:
      self.writer = Chem.SDWriter(gfile.GFile(output_path, 'w'))
    else:
      self.writer = Chem.SDWriter(sys.stdout)

  def write(self, conformer):
    """Writes a Conformer.

    Args:
      conformer: dataset_pb2.Conformer
    """
    for mol in smu_utils_lib.conformer_to_molecules(
        conformer,
        include_initial_geometries=self.init_geometry,
        include_optimized_geometry=self.opt_geometry,
        include_all_bond_topologies=True):
      self.writer.write(mol)

  def close(self):
    self.writer.close()


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  logging.info('Opening %s', FLAGS.input_sqlite)
  db = smu_sqlite.SMUSQLite(FLAGS.input_sqlite, 'r')
  if FLAGS.output_format == OutputFormat.pbtxt:
    outputter = PBTextOutputter(FLAGS.output_path)
  elif FLAGS.output_format == OutputFormat.sdf_init:
    outputter = SDFOutputter(
        FLAGS.output_path, init_geometry=True, opt_geometry=False)
  elif FLAGS.output_format == OutputFormat.sdf_opt:
    outputter = SDFOutputter(
        FLAGS.output_path, init_geometry=False, opt_geometry=True)
  elif FLAGS.output_format == OutputFormat.sdf_init_opt:
    outputter = SDFOutputter(
        FLAGS.output_path, init_geometry=True, opt_geometry=True)
  elif FLAGS.output_format == OutputFormat.proto:
    outputter = SerializedProtoOutputter(FLAGS.output_path)
  else:
    raise ValueError(f'Bad output format {FLAGS.output_format}')

  with contextlib.closing(outputter):
    for cid in (int(x) for x in FLAGS.cids):
      conformer = db.find_by_conformer_id(cid)
      outputter.write(conformer)
    for btid in (int(x) for x in FLAGS.btids):
      conformers = db.find_by_bond_topology_id(btid)
      if not conformers:
        raise KeyError(f'Bond topology {btid} not found')
      for c in conformers:
        outputter.write(c)
    for smiles in FLAGS.smiles:
      conformers = db.find_by_smiles(smiles)
      if not conformers:
        raise KeyError(f'SMILES {smiles} not found')
      for c in conformers:
        outputter.write(c)


if __name__ == '__main__':
  app.run(main)

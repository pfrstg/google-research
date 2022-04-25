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

"""Examines extending bond lengths from our initial topology detection."""

import collections
import csv
import itertools

from absl import app
from absl import logging

from smu import dataset_pb2

from smu import smu_sqlite
from smu.geometry import bond_length_distribution
from smu.geometry import smu_molecule
from smu.geometry import topology_from_geom

def get_modified_bond_lengths(epsilon):
  orig_bond_lengths = bond_length_distribution.AllAtomPairLengthDistributions()
  orig_bond_lengths.add_from_sparse_dataframe_file(
    '20220128_bond_lengths.csv',
    bond_length_distribution.STANDARD_UNBONDED_RIGHT_TAIL_MASS,
    bond_length_distribution.STANDARD_SIG_DIGITS)

  bond_lengths = bond_length_distribution.AllAtomPairLengthDistributions()
  for atom_a, atom_b in itertools.combinations_with_replacement(
      [dataset_pb2.BondTopology.ATOM_C,
       dataset_pb2.BondTopology.ATOM_N,
       dataset_pb2.BondTopology.ATOM_O,
       dataset_pb2.BondTopology.ATOM_F], 2):
    for bond in [dataset_pb2.BondTopology.BOND_UNDEFINED,
                 dataset_pb2.BondTopology.BOND_SINGLE,
                 dataset_pb2.BondTopology.BOND_DOUBLE,
                 dataset_pb2.BondTopology.BOND_TRIPLE]:
      if not bond_length_distribution.is_valid_bond(atom_a, atom_b, bond):
        continue
      try:
        dist = orig_bond_lengths[(atom_a, atom_b)][bond]
        if bond == dataset_pb2.BondTopology.BOND_UNDEFINED:
          bond_lengths.add(
            atom_a, atom_b, bond,
            bond_length_distribution.FixedWindow(dist.min() - epsilon,
                                                 2.0,
                                                 dist.right_tail_mass))
        else:
          bond_lengths.add(
            atom_a, atom_b, bond,
            bond_length_distribution.FixedWindow(dist.min() - epsilon,
                                                 dist.max() + epsilon,
                                                 None))
      except KeyError:
        # We have a few missing cases from our original empirical dists
        # For this exercise, we will just copy the CC bond dists for this order
        bond_lengths.add(
          atom_a, atom_b, bond,
          orig_bond_lengths[(dataset_pb2.BondTopology.ATOM_C,
                             dataset_pb2.BondTopology.ATOM_C)][bond])

  return bond_lengths


def main(argv):
  db = smu_sqlite.SMUSQLite('20220128_complete_v2.sqlite')

  buffers = [0, 0.01, 0.025, 0.05]
  bond_lengths = {buf: get_modified_bond_lengths(buf) for buf in buffers}
  for dists in bond_lengths.values():
    bond_length_distribution.add_itc_h_lengths(dists)
  smiles_id_dict = db.get_smiles_id_dict()
  matching_parameters = smu_molecule.MatchingParameters()
  matching_parameters.check_hydrogen_dists = True

  count_processed = 0
  count_matched = 0

  with open('extend_bond_dists.csv', 'w') as outf:
    fields = ['conformer_id']
    for buf in buffers:
      fields.append(f'is_matched_{buf}')
      fields.append(f'num_matched_{buf}')
    writer = csv.DictWriter(outf, fields)
    writer.writeheader()

    for conformer in db:
    # for conformer in [db.find_by_conformer_id(375986006)]:
      if conformer.fate != dataset_pb2.Conformer.FATE_DISASSOCIATED:
        continue

      count_processed += 1
      if count_processed % 25000 == 0:
        logging.info(f'Processed {count_processed}, matched {count_matched}')

      row = {'conformer_id': conformer.conformer_id}
      any_match = False

      for buf in buffers:

        matches = topology_from_geom.bond_topologies_from_geom(
          conformer,
          bond_lengths=bond_lengths[buf],
          matching_parameters=matching_parameters)

        matching_bt = [smiles_id_dict[bt.smiles]
                       for bt in matches.bond_topology]
        is_matched = (
          conformer.bond_topologies[0].bond_topology_id in matching_bt)

        # if matches.bond_topology:
        #   logging.info('For %d, bt %d, got %s',
        #                conformer.conformer_id,
        #                conformer.bond_topologies[0].bond_topology_id,
        #                str(matching_bt))

        row[f'is_matched_{buf}'] = is_matched
        row[f'num_matched_{buf}'] = len(matching_bt)
        any_match = any_match or is_matched

      if any_match:
        writer.writerow(row)
        count_matched += 1
        # if count_matched > 1000:
        #   break

  print(f'Final stats: {count_matched} / {count_processed}')

if __name__ == '__main__':
  app.run(main)

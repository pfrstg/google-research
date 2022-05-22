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

"""Use the various indices in the database for fast lookups."""

from smu import smu_sqlite
from smu.parser import smu_utils_lib

db = smu_sqlite.SMUSQLite('20220128_standard_v2.sqlite')

print('There are several ways to efficiently get specific sets of molecules')

print()
print('First is a lookup by molecule id')
cid_molecule = db.find_by_molecule_id(57001)
print('Looking up 57001 returns molecule with id', cid_molecule.molecule_id,
      'and bond topology with SMILES', cid_molecule.bond_topologies[0].smiles)

try:
  db.find_by_molecule_id(999999)
except KeyError:
  print('Looking up a molecule id not in the DB raises a KeyError')

print()
print('Looking up by bond topology id will return zero or more molecules')
bt_molecules = list(db.find_by_bond_topology_id_list(
  [7984],
  which_topologies=smu_utils_lib.WhichTopologies.all))
print('Querying for bond topology id 8617 returned',
      len(bt_molecules),
      'molecules')

print('Note that the molecules returned may have multiple bond topologies,'
      'and may or may not have the requested bond topology first')
for conf in bt_molecules:
  print('    Result with molecule_id', conf.molecule_id)
  for bt in conf.bond_topologies:
    print('        has bond topology with id', bt.bond_topology_id,
          'and SMILES', bt.smiles)

print()
print(
    'Finding by SMILES is essentially equivalent to finding by bond topology id'
)
smiles_molecules = list(db.find_by_smiles_list(
  ['O=NONNNO'],
  which_topologies=smu_utils_lib.WhichTopologies.all))
print('With query O=NONNNO', 'we found', len(smiles_molecules), 'results')

print('Note that the SMILES are canonicalized internally, you do not need to')
print('So the equivalent SMILES query ONNNON=O returns the same',
      len(list(db.find_by_smiles_list(
        ['ONNNON=O'],
        which_topologies=smu_utils_lib.WhichTopologies.all))),
      'results')

print()
print('You can also find all the molecules with a given stoichiometry')
stoich_molecules = list(db.find_by_stoichiometry('cn2o3'))
print('For example, "cn2o3" finds', len(stoich_molecules), 'results')
print('The first couple of molecule ids are:',
      [c.molecule_id for c in stoich_molecules[0:5]])

print()
print('You may note that there is a "find_by_expanded_stoichiometry" method',
      'in smu_sqlite')
print('This is primarily intended to support the "topology queries" that are',
      'documented in query_sqlite.py')
print('Since these topology queries are more involved to set up, it is',
      'recommended that you use query_sqlite.py for that kind of query')

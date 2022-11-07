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
"""Basic access to different kinds of fields."""

from smu import smu_sqlite

db = smu_sqlite.SMUSQLite('20220621_standard.sqlite')

#-----------------------------------------------------------------------------
# This is an arbitrary choice of the molecule to use.
#-----------------------------------------------------------------------------
molecule = db.find_by_mol_id(57001)

print('We will examine molecule with id', molecule.mol_id)

print('The computed properties are generally in the .prop field')

print()
print('Scalar values are access by name (note the .value suffix),',
      'like this single point energy: ',
      molecule.prop.spe_comp_b5.value)

print()
print('Fields with repeated values',
      'like vib_intens and exc_ene_cc2_tzvp)',
      'use an index with [] on the repeated values')

print('The 0th and 6th vib_intens:',
      molecule.prop.vib_intens.value[0],
      molecule.prop.vib_intens.value[6])

print('Or you can iterate over all of them')
for value in molecule.prop.exc_ene_cc2_tzvp.value:
  print('Excitation energy:', value)

print('Or just ask how many exc_ene_cc2_tzvp there are:',
      len(molecule.prop.exc_ene_cc2_tzvp.value))

print()
print(
    'Some fields like elec_dip_pbe0_augpc1 have explicit x,y,z components'
)

print(molecule.prop.elec_dip_pbe0_augpc1.x,
      molecule.prop.elec_dip_pbe0_augpc1.y,
      molecule.prop.elec_dip_pbe0_augpc1.z)

print()
print('Some fields are higher order tensors, with similar named components')

print('polarizability is a rank 2 tensor')
print(molecule.prop.elec_pol_pbe0_augpc1.xx,
      molecule.prop.elec_pol_pbe0_augpc1.yy,
      molecule.prop.elec_pol_pbe0_augpc1.zz,
      molecule.prop.elec_pol_pbe0_augpc1.xy,
      molecule.prop.elec_pol_pbe0_augpc1.xz,
      molecule.prop.elec_pol_pbe0_augpc1.yz)

print('octopole moment is a rank 3 tensor')
print(molecule.prop.elec_oct_pbe0_augpc1.xxx,
      molecule.prop.elec_oct_pbe0_augpc1.yyy,
      molecule.prop.elec_oct_pbe0_augpc1.zzz,
      molecule.prop.elec_oct_pbe0_augpc1.xxy,
      molecule.prop.elec_oct_pbe0_augpc1.xxz,
      molecule.prop.elec_oct_pbe0_augpc1.xyy,
      molecule.prop.elec_oct_pbe0_augpc1.yyz,
      molecule.prop.elec_oct_pbe0_augpc1.xzz,
      molecule.prop.elec_oct_pbe0_augpc1.yzz,
      molecule.prop.elec_oct_pbe0_augpc1.xyz)

print()
print('A couple of important fields are not inside "properties"')

geometry = molecule.optimized_geometry
print('For example, the optimized geometry has an energy of',
      geometry.energy.value,
      'and positions for',
      len(geometry.atom_positions),
      'atoms and the first atom x-coordinate is',
      geometry.atom_positions[0].x)

print()
print('In addition to looking at dataset.proto for field documentation,',
      'you can just print a given molecule to see what fields are available')

print(molecule)

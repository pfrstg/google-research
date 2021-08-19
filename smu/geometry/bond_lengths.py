"""Extract bond lengths from SMU molecules."""

from typing import Dict, Tuple

import apache_beam as beam
import numpy as np

import utilities

from smu import dataset_pb2
from smu.parser import smu_utils_lib


MAX_DIST = 2.0


def bond_lengths(conformer: dataset_pb2.Conformer,
                 data: Dict[Tuple[int, int, int], np.array]):
  """Accumulate bond lengths from `conformer` into `data`.
  Args:
    conformer:
    data: A dict from tuples of (atomic number, bond type, atomic_number) to
          lists of counts of discrtized distances.
          Note there is an ordering requirement on the key, see below.
  """
  bt = conformer.bond_topologies[0]
  geom = conformer.optimized_geometry

  bonded = utilities.bonded(bt)

  natoms = len(bt.atoms)

  for a1 in range(0, natoms):
    atomic_number1 = smu_utils_lib.ATOM_TYPE_TO_ATOMIC_NUMBER[bt.atoms[a1]]
    for a2 in range(a1 + 1, natoms):
      atomic_number2 = smu_utils_lib.ATOM_TYPE_TO_ATOMIC_NUMBER[bt.atoms[a2]]
      # Do not process H-H pairs
      if atomic_number1 == 1 and atomic_number2 == 1:
        continue

      d = utilities.distance_between_atoms(geom, a1, a2)
      if d > MAX_DIST:
        continue

      discretized = int(d * utilities.DISTANCE_BINS)
      key = (min(atomic_number1, atomic_number2),
             int(bonded[a1, a2]),
             max(atomic_number1, atomic_number2))
      data[key][discretized] += 1

def identify_range(data: np.array) -> Tuple[int, int]:
  """Return a tuple containing the index of the first and
  last non zero items in `data`.
  """
  first_non_zero = np.argmax(data > 0)
  last_non_zero = data.shape[0] - np.argmax(np.flip(data) > 0)
  return (first_non_zero, last_non_zero)

def write_bond_length_distributions(data, output_fname):
  """
  """
  for key, value in data.items():
    n = np.sum(value)
    if n == 0:
      print(f"{key} no data")
      continue

    (first_non_zero, last_non_zero) = identify_range(value)
    print(f"{key} range {first_non_zero} to {last_non_zero} N={np.sum(value)}")
    if first_non_zero == last_non_zero:
      continue
    fname = f"{output_fname}.{key[0]}.{key[1]}.{key[2]}"
    with open(fname, "w") as output:
      for i in range(first_non_zero, last_non_zero + 1):
        d = i / utilities.DISTANCE_BINS
        print(f"{d:.3f} {value[i]}", file=output)

class GetBondLengthDistribution(beam.DoFn):
  """Generates a bond length distribution."""
  def process(self, conformer:dataset_pb2.Conformer):
    bt = conformer.bond_topologies[0]
    geom = conformer.optimized_geometry

    bonded = utilities.bonded(bt)

    natoms = len(bt.atoms)

    for a1 in range(0, natoms):
      atomic_number1 = smu_utils_lib.ATOM_TYPE_TO_ATOMIC_NUMBER[bt.atoms[a1]]
      for a2 in range(a1 + 1, natoms):
        atomic_number2 = smu_utils_lib.ATOM_TYPE_TO_ATOMIC_NUMBER[bt.atoms[a2]]
        # Do not process H-H pairs
        if atomic_number1 == 1 and atomic_number2 == 1:
          continue

        d = utilities.distance_between_atoms(geom, a1, a2)
        if d > MAX_DIST:
          continue

        discretized = int(d * utilities.DISTANCE_BINS)
        yield (min(atomic_number1, atomic_number2),
               int(bonded[a1, a2]),
               max(atomic_number1, atomic_number2), discretized), 1

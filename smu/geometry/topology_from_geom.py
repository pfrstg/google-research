#"""Functions related to discerning the BondTopology from the geometry."""
import itertools

from google.protobuf import text_format
from typing import Dict, List, Optional, Tuple

import numpy as np

import apache_beam as beam

import utilities

from smu import dataset_pb2
from smu.geometry import bond_length_distribution
from smu.geometry import smu_molecule
from smu.geometry import utilities
from smu.parser import smu_utils_lib
import smu_molecule

# The longest distance considered (Angstroms).
THRESHOLD = 2.0


def hydrogen_to_nearest_atom(bond_topology: dataset_pb2.BondTopology,
                             distances: np.array) -> Optional[dataset_pb2.BondTopology]:
  """Generate a BondTopology that joins each Hydrogen atom to its nearest
      heavy atom.
    Note that if we detect an H atom that cannot be attached to a heavy
    atom, we return None.
  Args:
    bond_topology:
    distances: natoms*natoms array of interatomic distances.
  Returns:
  """
  result = dataset_pb2.BondTopology()
  result.atoms[:] = bond_topology.atoms
  natoms = len(bond_topology.atoms)
  for a1 in range(0, natoms):
    if bond_topology.atoms[a1] != dataset_pb2.BondTopology.AtomType.ATOM_H:
      continue

    shortest_distance = 1.0e+30
    closest_heavy_atom = -1
    for a2 in range(0, natoms):
      if bond_topology.atoms[a2] == dataset_pb2.BondTopology.AtomType.ATOM_H:
        continue

      if distances[a1, a2] >= THRESHOLD:
        continue

      if distances[a1, a2] < shortest_distance:
        shortest_distance = distances[a1, a2]
        closest_heavy_atom = a2

    if closest_heavy_atom < 0:
      return None

    bond = dataset_pb2.BondTopology.Bond(atom_a=a1,
                                         atom_b=closest_heavy_atom,
                                         bond_type=dataset_pb2.BondTopology.BondType.BOND_SINGLE)
    result.bonds.append(bond)

  return result


def indices_of_heavy_atoms(bond_topology: dataset_pb2.BondTopology) -> List[int]:
  """Return the indices of the heavy atoms in `bond_topology`.
    Args:
    Returns:
  """
  return [
      i for i, t in enumerate(bond_topology.atoms) if t != dataset_pb2.BondTopology.AtomType.ATOM_H
  ]

def report_atoms_too_close(distances:np.array, cutoff, bt:dataset_pb2.BondTopology):
  """Report atom pairs in `distances` closer than `cutoff`.
  Args:
    distances: pairwise distances between the atoms in `bt`.
    cutoff: report distances closer than `cutoff`.
    bt: BondTopology
  """
  print("Atoms too close")
  natoms = np.shape(distances)[0]
  for i in range(0,natoms):
    for j in range(i + 1,natoms):
      if distances[i, j] < cutoff:
        s1 = smu_utils_lib.ATOM_TYPE_TO_RDKIT[bt.atoms[i]]
        s2 = smu_utils_lib.ATOM_TYPE_TO_RDKIT[bt.atoms[j]]
        print(f"  too short {s1} {s2} {distances[i,j]}")

def atoms_too_close(distances:np.array, cutoff) -> bool:
  """Return true if any values in `distances` are less than `cutoff`.
  The diagonals are all zero, so will be < cutoff, so if there are
  more than the size of `distances` less than cutoff, those must be
  off diagonals.
  """
  return np.sum(np.less(distances, cutoff)) > np.shape(distances)[0]


def bond_topologies_from_geom(
    bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions,
    conformer_id: int,
    bond_topology: dataset_pb2.BondTopology, geometry: dataset_pb2.Geometry,
    single_point_energy_pbe0d3_6_311gd: dataset_pb2.Properties.ScalarMolecularProperty,
    matching_parameters: smu_molecule.MatchingParameters) -> dataset_pb2.TopologyMatches:
  """Return all BondTopology's that are plausible.

    Given a molecule described by `bond_topology` and `geometry`, return all possible
    BondTopology that are consistent with that.
    Note that `bond_topology` will be put in a canonical form.
    Args:
      bond_length_distribution:
      bond_topology:
      geometry:
      single_point_energy_pbe0d3_6_311gd:
    Returns:
      TopologyMatches
  """
  result = dataset_pb2.TopologyMatches()    # To be returned.
  result.conformer_id = conformer_id
  result.single_point_energy_pbe0d3_6_311gd = single_point_energy_pbe0d3_6_311gd.value
  if len(bond_topology.atoms) <= 1:
    return result    # empty.
  # Return empty result if there is no geometry data. Should not happen.
  if len(geometry.atom_positions) == 0:
    return result

  utilities.canonical_bond_topology(bond_topology)
  distances = utilities.distances(geometry)
  if atoms_too_close(distances, 0.90):
    report_atoms_too_close(distances, 0.90, bond_topology)
    return result

  # First join each Hydrogen to its nearest heavy atom, thereby
  # creating a starting BondTopology from which all others can grow
  starting_bond_topology = hydrogen_to_nearest_atom(bond_topology, distances)
  if starting_bond_topology is None:
    return result

# heavy_atoms = [a for a in bond_topology.atoms if a != dataset_pb2.BondTopology.AtomType.ATOM_H]
# heavy_atom_indices = indices_of_heavy_atoms(bond_topology)
  heavy_atom_indices = [
      i for i, t in enumerate(bond_topology.atoms) if t != dataset_pb2.BondTopology.AtomType.ATOM_H
  ]

  # For each atom pair, a list of possible bond types.
  # Key is a tuple of the two atom numbers, value is an np.array
  # with the score for each bond type.

  bonds_to_scores: Dict[Tuple[int, int], np.array] = {}
  for i, j in itertools.combinations(heavy_atom_indices, 2):    # All pairs.
    dist = distances[i, j]
    if dist > THRESHOLD:
      continue
    btypes = np.zeros(4, np.float32)
    for btype in range(0, 4):
      btypes[btype] = bond_lengths.pdf_length_given_type(bond_topology.atoms[i],
                                                         bond_topology.atoms[j], btype, dist)

    if np.count_nonzero(btypes) > 0:
      bonds_to_scores[(i, j)] = btypes

  if not bonds_to_scores:    # Seems unlikely.
    return result

  mol = smu_molecule.SmuMolecule(starting_bond_topology, bonds_to_scores, matching_parameters)

  search_space = mol.generate_search_state()
  for s in itertools.product(*search_space):
    bt = mol.place_bonds(list(s))
    if not bt:
      continue

    utilities.canonical_bond_topology(bt)
    if utilities.same_bond_topology(bond_topology, bt):
      bt.is_starting_topology = True
    bt.smiles = smu_utils_lib.compute_smiles_for_bond_topology(bt,
                                                               include_hs=True,
                                                               labeled_atoms=True)
    compute_probability(bt, bond_lengths, distances)
    result.bond_topology.append(bt)

  if len(result.bond_topology) > 1:
    result.bond_topology.sort(key=lambda bt: bt.goodness_of_fit, reverse=True)

  return result

def compute_probability(bt: dataset_pb2.BondTopology,
                        bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions,
                        distances: np.array):
  """Set the topology_probability attribute in `bt` based on the bonding
  pattern in `bt`.
  Args:
    bt:
    bond_lengths: bond length distributions
    distances: matrix of interatomic distances
  """
  result = 1.0
  for bond in bt.bonds:
    at1 = bt.atoms[bond.atom_a]
    at2 = bt.atoms[bond.atom_b]

    pdfs = bond_lengths.probability_of_bond_types(at1, at2, distances[bond.atom_a, bond.atom_b])
    if bond.bond_type in pdfs:
      result *= pdfs[bond.bond_type]
    else:
      result = 0.0

  bt.topology_probability = result


class TopologyFromGeom(beam.DoFn):
  """Beam class for extracting BondTopology from Conformer protos."""

  def __init__(self, bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions):
    super().__init__()
    self._bond_lengths = bond_lengths

  def process(self, conformer: dataset_pb2.Conformer):
    """Called by Beam.
      Returns a TopologyMatches for the plausible BondTopology's in `conformer`.
    Args:
      conformer:
    Yields:
      dataset_pb2.TopologyMatches
    """
    matching_parameters = smu_molecule.MatchingParameters()
    yield bond_topologies_from_geom(self._bond_lengths, conformer.conformer_id,
                                    conformer.bond_topologies[0],
                                    conformer.optimized_geometry,
                                    conformer.properties.single_point_energy_pbe0d3_6_311gd,
                                    matching_parameters)

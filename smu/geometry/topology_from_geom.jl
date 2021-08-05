using Combinatorics

using dataset_pb2

using BondLengthDistributions
using SmuMolecules
using SmuUtilities


# Atoms separated by > THRESHOLD are not further examined.
const THRESHOLD = 2.0
# In order to attach a Hydrogen to a nearby heavy atom, the
# two must be within HYDROGEN_THRESHOLD
const HYDROGEN_THRESHOLD = 1.4

"""Generate a BondTopology that joins each Hydrogen atom to its nearest
    heavy atom.
Args:
  bond_topology:
  distances:
Returns:
  A newly created BondTopology, based on `bond_topology`, but with 
  newly added bonds to Hydrogens - if found.
"""
function hydrogen_to_nearest_atom(bond_topology::dataset_pb2.BondTopology,
                             distances::Array{Float32,2})::Union{dataset_pb2.BondTopology, Nothing}
  result = dataset_pb2.BondTopology()
  # Not sure if this is safe or not.
  setproperty!(result, :atoms, bond_topology.atoms)

  single_bond = 1

  natoms = length(bond_topology.atoms)
  for a1 in 1:natoms
    bond_topology.atoms[a1] == dataset_pb2.BondTopology_AtomType.ATOM_H || continue

    shortest_distance = typemax(Float32)
    closest_heavy_atom = -1
    for a2 in 1:natoms
      bond_topology.atoms[a2] == dataset_pb2.BondTopology_AtomType.ATOM_H && continue
      a1 == a2 && continue

      distances[a1, a2] >= HYDROGEN_THRESHOLD && continue

      if distances[a1, a2] < shortest_distance
        shortest_distance = distances[a1, a2]
        closest_heavy_atom = a2
      end
    end

    closest_heavy_atom < 0 && return nothing

    add_bond!(a1 - 1, closest_heavy_atom - 1, single_bond, result)
  end

  @debug("With H's attached have $(length(result.bonds)) bonds")
  return result
end

# Return a Vector of all atom numbers that are within `threshold` of atom `atom`
# as defined by `geometry`.
# Args:
#  bt: BondTopology, only used for identifying H and non H atoms
#  geometry: the coordinates of each atom
#  atom: the atom for which we are determining neighbours.
#  threshold: only capture other atoms within this distance of `atom`.
#  distances: precomputed interatomic distance matrix.
function plausible_connections(bt::dataset_pb2.BondTopology,
                               geometry::dataset_pb2.Geometry,
                               atom::Integer,
                               threshold::Float32,
                               distances::Array)::Vector{Integer}
  result = Vector{Integer}()
  natoms = length(geometry.atom_position)
  for i in 1:natoms
    for j in (i+1):natoms
      distances[i, j] > threshold && continue
      push!(result, j)
    end
  end
  result
end

# Enable sorting of BondTopology's
import Base: isless
Base.isless(bt1::BondTopology, bt2::BondTopology) = isless(bt1.score, bt2.score)

# Convenenience function for filtering scores.
nonzero = x->x>0.0

# Convenenience function for selecting non hydrogen atoms
not_hydrogen = x->x != BondTopology_AtomType.ATOM_H

# Return a Dict relating possibly bonded atoms to a Vector of scores.
# Only the atom pairs implied by `atoms_to_examine` are considered.
# The keys to the Dict are tuples of atom numbers within `bond_topology`.
# Precomputed interatomic distances are in `distances`.
# The values of the Dict are Vectors of scores for each bond type if
# those two atoms were bonded with a bond corresponding to the index.
function atom_pairs_to_scores(bond_lengths::AllBondLengthDistributions,
                              bond_topology::dataset_pb2.BondTopology,
                              distances::Matrix,
                              atoms_to_examine::Vector)::Dict{Tuple{Int32,Int32}, Vector{Float32}}

  bonds_to_scores = Dict{Tuple{Int32,Int32}, Vector{Float32}}()  # to be returned.

  for (i, j) in Combinatorics.combinations(atoms_to_examine, 2)  # All pairs
    @debug("Between atoms $i and $j distance $(distances[i,j])")
    (dist = distances[i, j]) > THRESHOLD && continue

    scores = [pdf(bond_lengths, bond_topology.atoms[i], bond_topology.atoms[j], btype, dist)
              for btype in 0:4]

    @debug("Scores btw $(i) and $(j) $(scores)")
    any(nonzero, scores) && (bonds_to_scores[(i, j)] = scores)
  end
  bonds_to_scores
end


"""Return all BondTopology's that are plausible.

  Given a molecule described by `bond_topology` and `geometry`, return all possible
  BondTopology that are consistent with that.
  Note that `bond_topology` will be put in a canonical form.
  Args:
    bond_length_distribution:
    bond_topology:
    geometry:
  Returns:
    TopologyMatches
"""
function bond_topologies_from_geom(
    bond_lengths::AllBondLengthDistributions,
    bond_topology::dataset_pb2.BondTopology,
    geometry::dataset_pb2.Geometry,
    matching_parameters::MatchingParameters)::dataset_pb2.TopologyMatches
  result = dataset_pb2.TopologyMatches()    # To be returned.
  @debug("have $(length(bond_topology.atoms)) atoms and $(length(bond_topology.bonds)) bonds")
  @debug(bond_topology)
  length(bond_topology.atoms) <= 1 && return result  # return an empty result
  length(geometry.atom_positions) == 0 && return result
# hasproperty(bond_topology, :smiles) && println(bond_topology.smiles)

  canonical_bond_topology!(bond_topology)
  distances = SmuUtilities.distances(geometry)
  @debug(distances)
  longest_supposedly_bonded_distance(bond_topology, distances) > THRESHOLD && return result

  # First join each Hydrogen to its nearest heavy atom, thereby
  # creating a starting BondTopology from which all others can grow
  starting_bond_topology = hydrogen_to_nearest_atom(bond_topology, distances)
  starting_bond_topology === nothing && @debug("Hydrogens not attached")
  starting_bond_topology === nothing && return result
  @debug("Hydrogens joined OK, have $(length(starting_bond_topology.bonds)) bonds")

  heavy_atom_indices = findall(not_hydrogen, bond_topology.atoms)
  @debug("heavy_atom_indices $(heavy_atom_indices)")
  length(heavy_atom_indices) < 2 && return result

  # For each atom pair, a list of possible bond types.
  # Key is a tuple of the two atom numbers, value is a Vector
  # with the score for each bond type.

  bonds_to_scores = atom_pairs_to_scores(bond_lengths, bond_topology, distances, heavy_atom_indices)

  @debug("bonds_to_scores $(length(bonds_to_scores))")
  isempty(bonds_to_scores) && return result

  found_topologies = Vector{BondTopology}()  # Will be set into `result`.

  mol = SmuMolecule(starting_bond_topology, bonds_to_scores, matching_parameters)

  search_space = generate_search_state(mol)
  @debug("search_space $search_space")
  for s in Iterators.product(search_space...)
    @debug("Placing state $s")
    (bt = place_bonds!(s, mol)) == nothing && continue
    is_single_fragment(bt) || @debug("not single fragment")
    is_single_fragment(bt) || continue

    canonical_bond_topology!(bt)
    bt.is_starting_topology = same_bond_topology(bond_topology, bt)
    bt.smiles = to_smiles(bt, SmilesOptions(true))
    push!(found_topologies, bt)
  end

  @debug("found any topologies $(length(found_topologies))")
  isempty(found_topologies) && return result

  length(found_topologies) > 1 && sort!(found_topologies, rev=true)

  setproperty!(result, :bond_topology, found_topologies)

  return result
end

# Add a new bond to `bt`, involving atoms `a1` and `a2` of type `btype`.
# bt should be one of dataset_pb2.BondTopology_BondType
function new_bond(a1::Integer, a2::Integer, bt::Integer)::dataset_pb2.BondTopology_Bond
  dataset_pb2.BondTopology_Bond(atom_a=a1, atom_b=a2, bond_type=bt)
end

# Matching is complete if every atom needs no more bonds.
function matching_complete(bonds_needed::Vector)::Bool
  findfirst(x->x > 0, bonds_needed) === nothing
end

# If there are atoms with just one possible other atom connected,
# add those bonds to `bt`.
# The idea is that this should be called repeatedly, until no
# change is observed.
# Args:
#  bt: a BondTopology to which bonds may be added.
#  bonds_needed: for each atom, how many bonds are needed
#  maybe_connected: for each atom, a list of the other atoms
#     that might be connected.
# Returns:
#  the number of bonds added.
function make_unambiguous_bonds!(bt::dataset_pb2.BondTopology,
                                bonds_needed::Vector{Integer},
                                maybe_connected::Vector)
  natoms = length(bt.atoms)
  result = 0
  for i in 1:natoms
    length(bonds_needed[i]) == 0 && continue  # Probably a Hydrogen
    length(maybe_connected[i]) > 1 && continue
    other_atom = maybe_connected[i][1]
    # Fatal error if the other atom has no open valence.
    bonds_needed[other_atom] == 0 && return result

    subtract = 0  # From bonds_needed[other_atom]
    if bonds_needed[i] == 1 # unambiguous, single bond only
      push!(bt.bonds, new_bond(i - 1, other_atom - 1, single_bond))
      subtract = 1
    elseif bonds_needed[i] == 2 && bonds_needed[other_atom] >= 2
      push!(bt.bonds, new_bond(i - 1, other_atom - 1, double_bond))
      subtract = 2
    elseif bonds_needed[i] == 3 && bonds_needed[other_atom] >= 3
      push!(bt.bonds, new_bond(i - 1, other_atom - 1, triple_bond))
      subtract = 3
    else
      continue
    end
    bonds_needed[i] = 0
    resize!(maybe_connected[i], 0)
    bonds_needed[other_atom] -= subtract
    filter!(x -> x ≠ i, maybe_connected[other_atom])
    result += 1
  end
  result
end

"""Return all BondTopology's that are plausible given a geometry and a set of bond
   length distributions.

  Only the coordinates in `geometry` are considered, the 'known' BondTopology in
  `bond_topology` is only used for checking.
  Note that `bond_topology` will be put in a canonical form.
  Args:
    bond_length_distribution:
    bond_topology:
    geometry:
  Returns:
    TopologyMatches
"""
function bond_topologies_from_geom_unconstrained(
    bond_lengths::AllBondLengthDistributions,
    bond_topology::dataset_pb2.BondTopology,
    geometry::dataset_pb2.Geometry,
    matching_parameters::MatchingParameters)::dataset_pb2.TopologyMatches
  result = dataset_pb2.TopologyMatches()    # To be returned.
  @debug("have $(length(bond_topology.atoms)) atoms and $(length(bond_topology.bonds)) bonds")
  @debug(bond_topology)
  length(bond_topology.atoms) <= 1 && return result  # return an empty result
  length(geometry.atom_positions) == 0 && return result

  canonical_bond_topology!(bond_topology)

  distances = SmuUtilities.distances(geometry)

  # For each atom, a list of other atoms that might be connected.
  maybe_connected = [plausible_connections(bond_topology, geometry, i, threshold, distances) for i in 1:natoms]

  # Any bond involving these types can only be single.
  # single_valent = [dataset_pb2.BondTopology_AtomType.ATOM_H, dataset_pb2.BondTopology_AtomType.ATOM_F, dataset_pb2.BondTopology_AtomType.ATOM_OMINUS]

  hydrogen = dataset_pb2.BondTopology_AtomType.ATOM_H

  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  double_bond = dataset_pb2.BondTopology_BondType.BOND_DOUBLE
  triple_bond = dataset_pb2.BondTopology_BondType.BOND_TRIPLE

  # First join each Hydrogen to its nearest heavy atom, thereby
  # creating a starting BondTopology from which all others can grow
  starting_bond_topology = hydrogen_to_nearest_atom(bond_topology, distances)
  starting_bond_topology === nothing && @debug("Hydrogens not attached")
  starting_bond_topology === nothing && return result
  @debug("Hydrogens joined OK, have $(length(starting_bond_topology.bonds)) bonds")

  # For each atom, the number of connections it needs.
  bonds_needed = [smu_atom_type_to_max_con(i) for i in bond_topology.atoms]

  # Decrement the bonds_needed attribute for each atom involved with a bond to H
  for bond in bt.bonds
    a1 = bond.atom_a + 1
    a2 = bond.atom_b + 1
    if bt.atoms[a1] == hydrogen || bt.atoms[a2] == hydrogen
      bonds_needed[a1] -= 1
      bonds_needed[a2] -= 1
    end
  end

  # Nothing can have negative connectivity.
  if findfirst(x->x < 0, connections_needed) !== nothing
    @warn("Negative connection after H placement in $(bond_topology.conformer_id)")
    return result   # empty.
  end

  # Repeatedly examine atoms with only one possible other atom connected.
  while true
    make_unambiguous_bonds(starting_bond_topology, bonds_needed, maybe_connected) > 0 || break
  end

  if matching_complete(bonds_needed)
    # actually do something here.
    return result
  end

  # Identify those atoms that need more connections
  varying_atoms = [i for i in 1:natoms if bonds_needed[i] > 0]

  length(varying_atoms) < 2 && return result  # Cannot make any bonds.

  # For each varying pair of atoms, a vector of scores per bond type.
  bonds_to_scores = atom_pairs_to_scores(bond_lengths, bond_topology, distances, varying_atoms)

  # None of the varying atoms are within range of each other.
  isempty(bonds_to_scores) && return result

  # If only one possibility, make those bonds. All single bonds for now.
  for i in 1:natoms
    bt.atoms[i] == hydrogen && continue
    length(maybe_connected[i]) == 0 && continue
    length(maybe_connected[i]) > 1 && continue
    push!(starting_bond_topology.bonds, dataset_pb2.BondTopology_Bond(atom_a=i-1, atom_b=maybe_connected[i][1], bond_type=single_bond))
    resize!(plausible_connections[i], 0)
  end

  nbonds = length(bt.bonds)

  # For each bond, the kinds of bonds to be tried. If empty, the bond in the BondTopology will be the
  # only possibility.
  bonds_to_be_tried = [Vector{Int32}() for i in 1:nbonds]

  for bond in bt.bonds
    a1 = bond.atom_a + 1
    a2 = bond.atom_b + 1
    a1 in single_valent && continue
    a1 in single_valent && continue
  end

  for i in 1:natoms
    bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_H && continue
    bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_F && continue
    bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_OMINUS && continue
  end

end

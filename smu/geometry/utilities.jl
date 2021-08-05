# utility functions needed for discerning topology from geometry

using Combinatorics

using dataset_pb2

export bohr_to_angstroms, angstroms_to_bohr, distance_between_atoms, bonded, distances
export add_atom!, add_bond!
export smu_atype_to_atomic_number, number_to_smu_btype, smu_btype_to_number, smu_atom_type_to_max_con
export atomic_number_to_smu_atom_type
export add_atom_by_atomic_number!, add_atom_by_type!
export distances, bonded
export canonical_bond_topology!, same_bond_topology
export is_single_fragment
export prettyprint
export longest_supposedly_bonded_distance
export to_smiles
export SmilesOptions
export remove_atoms!
export are_bonded

function bohr_to_angstroms(distance)
  return distance * 0.529177249
end

function angstroms_to_bohr(distance)
  return distance / 0.529177249
end

function distance_between_atoms(geom::dataset_pb2.Geometry, a1, a2)::Float32
  return bohr_to_angstroms(sqrt(
    (geom.atom_positions[a1].x - geom.atom_positions[a2].x) *
    (geom.atom_positions[a1].x - geom.atom_positions[a2].x) +
    (geom.atom_positions[a1].y - geom.atom_positions[a2].y) *
    (geom.atom_positions[a1].y - geom.atom_positions[a2].y) +
    (geom.atom_positions[a1].z - geom.atom_positions[a2].z) *
    (geom.atom_positions[a1].z - geom.atom_positions[a2].z)
  ))
end

_smu_atype_to_atomic_number = Dict{Int32, Int32}(
  BondTopology_AtomType.ATOM_H => 1,
  BondTopology_AtomType.ATOM_C => 6,
  BondTopology_AtomType.ATOM_N => 7,
  BondTopology_AtomType.ATOM_NPOS => 7,
  BondTopology_AtomType.ATOM_O => 8,
  BondTopology_AtomType.ATOM_ONEG => 8,
  BondTopology_AtomType.ATOM_F => 9
)

function smu_atype_to_atomic_number(smu_type::T)::Int32 where {T<:Integer}
  return _smu_atype_to_atomic_number[smu_type]
end

"""Given an atomic number and charge, return the BondTopology_AtomType
Args:
  atomic_number:
  charge:
Returns:
  BondTopology_AtomType
"""
function atomic_number_to_smu_atom_type(atomic_number::Integer, charge::Integer=0)::Int32
  if atomic_number == 6
    return BondTopology_AtomType.ATOM_C
  elseif atomic_number == 7
    if charge == 0
      return BondTopology_AtomType.ATOM_N
    else
      return BondTopology_AtomType.ATOM_NPOS
    end
  elseif atomic_number == 8
    if charge == 0
      return BondTopology_AtomType.ATOM_O
    else
      return BondTopology_AtomType.ATOM_ONEG
    end
  elseif atomic_number == 9
    return BondTopology_AtomType.ATOM_F
  elseif atomic_number == 1
    return BondTopology_AtomType.ATOM_H
  else
    @error("Unrecognized atomic number $(atomic_number)")
  end
end

_number_to_smu_btype = Dict{Int32, Int32}(
  0 => BondTopology_BondType.BOND_UNDEFINED,
  1 => BondTopology_BondType.BOND_SINGLE,
  2 => BondTopology_BondType.BOND_DOUBLE,
  3 => BondTopology_BondType.BOND_TRIPLE
)

function number_to_smu_btype(n::T)::Int32 where {T<:Integer}
  return _number_to_smu_btype[n]
end

_smu_btype_to_number = Dict{Int32, Int32}(
  BondTopology_BondType.BOND_UNDEFINED => 0,
  BondTopology_BondType.BOND_SINGLE => 1,
  BondTopology_BondType.BOND_DOUBLE => 2,
  BondTopology_BondType.BOND_TRIPLE => 3
)

_smu_atom_type_to_max_con = Dict{Int32, Int32}(
  BondTopology_AtomType.ATOM_H => 1,
  BondTopology_AtomType.ATOM_C => 4,
  BondTopology_AtomType.ATOM_N => 3,
  BondTopology_AtomType.ATOM_NPOS => 3,
  BondTopology_AtomType.ATOM_O => 2,
  BondTopology_AtomType.ATOM_ONEG => 1,
  BondTopology_AtomType.ATOM_F => 1
)

# Return the maximum number of bonds for `smu_atype`.
# Unfortunately the type declaration on the next line does not work.
# function smu_atom_type_to_max_con(smu_atype::BondTopology_AtomType)::Int32

function smu_atom_type_to_max_con(smu_atype::T)::Int32 where {T <:Integer}
  return _smu_atom_type_to_max_con[smu_atype]
end

function smu_btype_to_number(smu_btype::T)::Int32 where {T<:Integer}
  return _smu_btype_to_number[smu_btype]
end


function add_bond!(a1::T, a2::T, btype, bond_topology::dataset_pb2.BondTopology) where {T<:Integer}
  smu_btype = number_to_smu_btype(btype)

  new_bond = BondTopology_Bond(atom_a=a1, atom_b=a2, bond_type=smu_btype)

  if !hasproperty(bond_topology, :bonds)
    setproperty!(bond_topology, :bonds, [new_bond])
  else
    push!(bond_topology.bonds, new_bond)
  end
end

function add_atom_by_atomic_number!(bond_topology::dataset_pb2.BondTopology, atomic_number::T, charge::T=0) where {T<:Integer}
  smu_atype = atomic_number_to_smu_atom_type(atomic_number, charge)
  add_atom_by_type!(bond_topology, smu_atype)
end

function add_atom_by_type!(bond_topology::dataset_pb2.BondTopology, smu_atype::T) where {T <:Integer}
  if !hasproperty(bond_topology, :atoms)
    setproperty!(bond_topology, :atoms, [smu_atype])
  else
    push!(bond_topology.atoms, smu_atype)
  end
end

# Remove all atoms that match `f` from `bt`.
function remove_atoms!(f::Function, bt::dataset_pb2.BondTopology)
   natoms = length(bt.atoms)
   atom_xref = fill(-1, natoms)
   ndx = 0
   for i in 1:natoms
     if ! f(bt.atoms[i])  # If matches, atom is removed. Not matches means keep.
       ndx += 1
       i > ndx && (bt.atoms[ndx] = bt.atoms[i])
       atom_xref[i] = ndx
     end
   end
 
   ndx == natoms && return  # Nothing matched, nothing to remove, done.

   resize!(bt.atoms, ndx)

   nbonds = length(bt.bonds)
   ndx = 0
   for i in 1:nbonds
     b = bt.bonds[i]
     a1 = b.atom_a
     a2 = b.atom_b
     if atom_xref[a1 + 1] > 0 && atom_xref[a2 + 1] > 0
       ndx += 1
       b.atom_a = atom_xref[a1 + 1] - 1
       b.atom_b = atom_xref[a2 + 1] - 1
       i > ndx && (bt.bonds[ndx] = bt.bonds[i])
     end
   end

  resize!(bt.bonds, ndx)
end

"""Return an int array of the bonded atoms in `bond_topology`.
The returned array is square[natoms,natoms]. If there is no bond
between the atoms, the value is zero. If the atoms are bonded,
the entry will be the bond type.
Args:
  bond_topology:
Returns:
  Vector of BondType's
"""
function bonded(bond_topology::dataset_pb2.BondTopology)::Array{Int32, 2}
  natoms = length(bond_topology.atoms)
  connected = zeros(Int32, (natoms, natoms))  # to be returned
  hasproperty(bond_topology, :bonds) || return connected
  # Note need to convert to 1 indexing for Julia
  for bond in bond_topology.bonds
    a1 = bond.atom_a + 1
    a2 = bond.atom_b + 1
    connected[a1, a2] = connected[a2, a1] = bond.bond_type
  end
  connected
end

"""Given a connection matrix `bonded`, return a vector of
vectors, one for each atom, containing the indices of the
bonded atoms.
Args:
  bonded:
Returns:
"""
function connections(bonded::Array{Int32, 2})::Vector{Vector{Int32}}
  natoms = size(bonded, 1)
  result = Vector{Vector{Int32}}(undef, (natoms))
  for i in 1:natoms
    nbrs = Vector{Int32}()
    for j in 1:natoms
      if bonded[i,j] > 0
        push!(nbrs, j)
      end
    end
    result[i] = nbrs
  end
  return result
end

function connections(bt::dataset_pb2.BondTopology)::Vector{Vector{Int32}}
  result = fill(Vector{Int32}(), length(bt.atoms))
  hasproperty(bt, :bonds) || return result
  for b in bt.bonds
    push!(result[b.atom_a], b.atom_b)
    push!(result[b.atom_b], b.atom_a)
  end
  result
end

"""Return a float array of the interatomic distances in `geometry`.
Args:
  geometry: 
Returns:
  A 2-d array of distances (natoms, natoms)
"""
function distances(geometry::dataset_pb2.Geometry)::Array{Float32, 2}
  natoms = length(geometry.atom_positions)
  result = zeros(Float32, (natoms, natoms))
  for atoms in combinations(1:natoms, 2)
    i = atoms[1]
    j = atoms[2]
    result[i, j] = result[j, i] = distance_between_atoms(geometry, i, j)
  end
  return result
end

"""Return the longest distance between a bonded pair of atoms.

Args:
  bt:: BondTopology where the atom pairs are in the .bonds attribute.
  distances: (natoms,natoms) array of precomputed distances.
Returns:
  The longest distance associated with a bonded pair.
"""
function longest_supposedly_bonded_distance(bt::BondTopology, distances)::Float32
  result = zero(Float32)
  for bond in bt.bonds
    dist = distances[bond.atom_a + 1, bond.atom_b + 1]
    dist <= result && continue

    result = dist
  end
  result
end

"""An expensive way to decide if two atoms are bonded.
  Note that a1 and a2 are atom numbers within a BondTopology, and so start at 0
"""
function are_bonded(bt::BondTopology, a1::Integer, a2::Integer)::Bool
  for bond in bt.bonds
    if a1 == bond.atom_a && a2 == bond.atom_b
      return true
    end
    if a1 == bond.atom_b && a2 == bond.atom_a
      return true
    end
  end
  false
end

function shortest_not_bonded_distance(bt::BondTopology, distances)::Float32
  result = typemax(Float32)
  natoms = length(bt.atoms)
  for i in 1:natoms
    for j in (i+1):natoms
      dist = distances[i,j]
      dist >= result && continue
      are_bonded(bt, i - 1, j - 1) && continue
      result = dist
    end
  end
  result
end


"""Transform the bonds attribute of `bond_topology` to a canonical form.

Args:
  bond_topology:
Returns:
  BondTopology
"""
function canonical_bond_topology!(bond_topology)
  if length(bond_topology.bonds) < 2
    return
  end

  for bond in bond_topology.bonds
    if bond.atom_a > bond.atom_b
      bond.atom_a, bond.atom_b = bond.atom_b, bond.atom_a
    end
  end

  sort!(bond_topology.bonds, by = b -> (b.atom_a, b.atom_b))
end

Base.:(==)(b1::BondTopology_Bond, b2::BondTopology_Bond)::Bool = b1.atom_a == b2.atom_a &&
                b1.atom_b == b2.atom_b && b1.bond_type == b2.bond_type

"""Return True if bt1 == bt2.
Note that there is no attempt to canonialise the protos here.
Args:
  bt1, bt2: BondTopology's
Returns:
  true if bt1 and bt2 are the same.
"""
function same_bond_topology(bt1, bt2)::Bool
  bt1.bonds == bt2.bonds || return false

  bt1.atoms == bt2.atoms
end

"""Recusrively visit nodes in the graph defined by `nbrs`.
Args:
  nbrs:
  atom:
  visited:
Returns:
  The number of nodes visited - including `atom`.
"""
function visit(nbrs::Vector, atom::T, visited::Vector{Bool})::Int32 where {T<:Integer}
  visited[atom] = true
  result = 1    # To be returned.
  for nbr in nbrs[atom]
    if visited[nbr] > 0
      continue
    end
    result += visit(nbrs, nbr, visited)
  end

  return result
end

"""Return True if `bond_topology` is a single fragment.
Args:
  bond_topology:
Returns:
  True if `bond_topology` is a single fragment.
"""
function is_single_fragment(bond_topology::BondTopology)::Bool
  natoms = length(bond_topology.atoms)
  nbonds = length(bond_topology.bonds)
  # Some special cases are easy.
  natoms == 1 && return true
  if natoms == 2 && nbonds == 1
    return true
  end
  if natoms == 3 && nbonds == 2
    return true
  end
  if natoms == nbonds && natoms <= 4
    return true
  end

  attached = connections(bonded(bond_topology))
  # Any atom with zero neighbors means a detached atom.
  any(n->length(n) == 0, attached) && return false

  visited = zeros(Bool, natoms)
  # Mark anything with a single connection as visited.
  # Record the index of an atom that has multiple connections.
  a_multiply_connected_atom = -1
  for i in 1:natoms
    if bond_topology.atoms[i] == BondTopology_AtomType.ATOM_H
      visited[i] = 1
      continue
    end

    if length(attached[i]) > 1
      a_multiply_connected_atom = i
      continue
    end

    # A singly connected heavy atom. Mark visited if not part of a two atom fragment.
    if length(attached[attached[i][1]]) > 1
      visited[i] = 1
    end
  end

  # Not sure this can happen.
  a_multiply_connected_atom < 0 && return false

  number_visited = sum(visited) + visit(attached, a_multiply_connected_atom, visited)
  return number_visited == natoms
end


function bond_as_string(btype)::String
  if btype == dataset_pb2.BondTopology_BondType.BOND_SINGLE
    return ""
  elseif btype == dataset_pb2.BondTopology_BondType.BOND_DOUBLE
    return "="
  elseif btype == dataset_pb2.BondTopology_BondType.BOND_TRIPLE
    return "#"
  end
end

import Base.show
function show(io::IO, bt::BondTopology)
  natoms = length(bt.atoms)
  println(io, "bt with $natoms atoms and $(length(bt.bonds)) bonds")
  for i in 1:natoms
    print(io, " $(i-1) $(bt.atoms[i]) type $(smu_atype_to_atomic_number(bt.atoms[i])) ->")
    for b in bt.bonds
      if (i-1) == b.atom_a
        print(io, " $(bond_as_string(b.bond_type))$(b.atom_b)")
      elseif (i-1) == b.atom_b
        print(io, " $(bond_as_string(b.bond_type))$(b.atom_a)")
      end
    end
    println(io)
  end
end

function prettyprint(bt::BondTopology)
  natoms = length(bt.atoms)
  println("bt with $natoms atoms and $(length(bt.bonds)) bonds")
  for i in 1:natoms
    print(" $(i-1) $(bt.atoms[i]) type $(smu_atype_to_atomic_number(bt.atoms[i])) ->")
    for b in bt.bonds
      if (i-1) == b.atom_a
        print(" $(bond_as_string(b.bond_type))$(b.atom_b)")
      elseif (i-1) == b.atom_b
        print(" $(bond_as_string(b.bond_type))$(b.atom_a)")
      end
    end
    println()
  end
end

"""Class to yield as many as 9 ring numbers via calls to get_ring_number.
  Ring numbers are not reused.
"""
mutable struct RingNumbers
  available::Vector{Bool}
  RingNumbers() = new(fill(true, 9))
end

function reset!(ring_numbers::RingNumbers)
  ring_numbers.available .= true
end

""" Return an unused ring number from `ring_numbers`.

  Args:
    ring_numbers: holds a list of available ring numbers.
"""
function get_ring_number(ring_numbers::RingNumbers)
  result = findfirst(ring_numbers.available)
  @assert result !== nothing "No ring numbers available"
  ring_numbers.available[result] = false
  result
end

mutable struct Connection
  atom::Integer
  btype::Integer
  Connection() = new(-1, 0)
  Connection(a::Integer, b::Integer) = new(a, b)
end

mutable struct ConnectionTable
  connections::Vector{Vector{Connection}}
  ConnectionTable() = ConnectionTable(Vector{Vector{Connection}()}())
  function ConnectionTable(natoms::Integer)
    res = new()
#   res.connections = Vector{Vector{Connection}}(undef, natoms)
    res.connections = [Vector{Connection}() for i in 1:natoms]
    res
  end
end

#natoms(c::ConnectionTable) = length(c.connections)

#"""Return the lowest atomic connectivity in `ctable`.
#"""
#function min_ncon(ctable::ConnectionTable)::Integer
#  reduce(x,y) -> min(x, length(y)), ctable.connections, init=typemax(Int32))
#end

mutable struct RingClosure
  ring_number::Integer
  bond_type::Integer
end

mutable struct AtomState
  complete::Bool
  # The atomic sniles
  
  smiles::String
  ring_closures::Vector{RingClosure}
  AtomState() = AtomState(false, "", Vector{RingClosure}())
end

function smiles_atom_symbol_not_used(bt::dataset_pb2.BondTopology, atom::Integer,
                            hydrogen::Vector)::String
  if bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_H
    return "[H]"
  elseif bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_C
    return "C"
  elseif bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_N
    return "N"
  elseif bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_NPOS
    return "[NH$(hydrogen[atom])+]"
  elseif bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_O
    return "O"
  elseif bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_ONEG
    return "[O-]"
  elseif bt.atoms[atom] == dataset_pb2.BondTopology_AtomType.ATOM_F
    return "F"
  end
end

#function smiles_atom_symbol(bt::dataset_pb2.BondTopology, atom::Integer,
#                            hydrogen::Vector,
#                            geometry::dataset_pb2.Geometry)::String
#  s = smiles_atom_symbol(bt, atom, hydrogen)
#  atom_pos = geometry.atom_pos[atom - 1]
#  return "$(s){{$(atom_pos.x),$(atom_pos.y),$(atom_pos.z)}}"
#end

function append_coordinates!(atom_pos::dataset_pb2.Geometry_AtomPos,
                            smiles::String)
  "$smiles{{$(atom_pos.x),$(atom_pos.y),$(atom_pos.z)}}"
end

"""Compute the implicit hydrogen count for N+ atoms.
  Maybe expand to all atoms sometime.
"""
function compute_hydrogen_count(bt::dataset_pb2.BondTopology)::Vector{Integer}
  natoms = length(bt.atoms)
  result = zeros(Integer, natoms)
  sums = sum(bonded(bt), dims=2)  # For each atom, total bonds.
  for i in 1:natoms
    if bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_H
      result[i] = 0
    elseif bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_C
      result[i] = 4 - sums[i]
    elseif bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_N
      result[i] = 3 - sums[i]
    elseif bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_NPOS
      result[i] = 4 - sums[i]
    elseif bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_O
      result[i] = 2 - sums[i]
    elseif bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_ONEG
      result[i] = 1 - sums[i]
    elseif bt.atoms[i] == dataset_pb2.BondTopology_AtomType.ATOM_F
      result[i] = 1 - sums[i]
    end
  end
  result
end

mutable struct RingFinderArrays
  # Whether or not each atom has been visited by the algorithm.
  visited::Vector{Integer}
  # The result of the calculation.
  ring_membership::Vector{Integer}
  # Is the bond between two atoms a ring bond or not?
  ring_bond_count::Vector{Integer}
  # The fragment membership of each atom
  fragment_membership::Vector{Integer}
  function RingFinderArrays(natoms)
    res = new()
    res.visited = fill(Integer, natoms)
    res.ring_membership = fill(Integer, natoms)
    res.ring_bond_count = fill(Integer, natoms * natoms)
    res.fragment_membership = fill(Integer, natoms)
    res
  end
end

function update_ring_bond_count(data, a1, a2)
  natoms = length(data.visited)
  data.ring_bond_count[a1 * natoms + a2] += 1
  data.ring_bond_count[a2 * natoms + a1] += 1
end

function ring_finder(ctable::ConnectionTable, 
                      previous_atom,
                      current_atom,
                      data::RingFinderArrays,
                      counter::Integer)::Integer
  data.visited[current_atom] = counter
  found_ring = 0
  for j in ctable.connections[current_atom]
    j == previous_atom && continue
    data.visited[j] > counter && continue
    if data.visited[j]
      data.ring_membership[j] += 1
      update_ring_bond_count(ring_bond_count, matoms, current_atom, j);
      found_ring += 1
    else
      tmp = RingFinder(ctable, current_atom, j, data, counter + 1);
      if tmp 
        found_ring += tmp
        update_ring_bond_count(data.ring_bond_count, current_atom, j);
      end
    end
  end

  # No evidence of any ring artifacts.
  found_ring == 0 && ring_membership[current_atom] == 0 && return 0

  if found_ring && ring_membership[current_atom] == 0
    data.ring_membership[current_atom] = found_ring
    return found_ring
  end

  # If we found all the rings that terminate here return 0;
  found_ring == data.ring_membership[current_atom] && return 0

  rc = found_ring - data.ring_membership[current_atom];
  data.ring_membership[current_atom] = found_ring;
  rc
end

# Given a connectivity, return info about fragments and rings.
# The tuple returned contains:
#   whether or not each atom is in a ring
#   the number of fragments
#   the fragment membership of each atom
function in_ring(ctable::ConnectionTable)::Tuple{Vector{Bool}, Integer, Vector{Integer}}
  natoms = length(ctable.connections)
  data = RingFinderArrays(natoms)

  counter = 1
  fragment_number = 0
  starting_atom = 1
  while true
    ring_finder(ctable, -1, starting_atom, data, counter);
    counter_start += update_fragment_information(data, starting_atom, fragment_number)
    fragment_number += 1
    starting_atom = findfirst(x->x == 0, data.visited)
    starting_atom || break
  end

  is_ring = map(x->x > 0, data.ring_membership)
  (result, fragment_number, data.fragment_membership)
end

# Given that a number of atoms have been newly identified as
# in a fragment, those with their data.visited attribute >= `min_visited_value`
# update the data.fragment_membership data with the fact that
# these atoms are all part of `fragment_number`.
function update_fragment_information(data::RingFinderArrays,
                fragment_number::Integer,
                min_visited_value::Integer)::Integer
  atoms_in_fragment = 0;
  natoms = length(data.visited)
  for i in 1:natoms
    data.visited[i] < min_visited_value && continue

    data.fragment_membership[i] = fragment_number
    atoms_in_fragment += 1
  end

  atoms_in_fragment
end

# Options describing optional aspects of smiles generation.
mutable struct SmilesOptions
  # atom map number = atom number
  numbered_atoms::Bool
  # 3d smiles
  append_coordinates::Bool
end
SmilesOptions() = SmilesOptions(false, false)
SmilesOptions(b) = SmilesOptions(b, false)

# If there are ring closures specified in `connection`, append that
# information to `smiles`.
# Args:
#  connection:
#  smiles:
function append_ring_closures_not_used_any_more(connection::Connection, smiles::IOStream)::String
  isempty(connection.ring_closures) && return smiles
  for rc in connection.ring_closures
    if rc.bond_type == 1
    elseif rc.bond_type == 2
      write(smiles, '=')
    elseif rc.bond_type == 3
      write(smiles, '#')
    end
    write(smiles, string(connection.atom))
  end
end

# Given a set of RingClosure's in `rings`, append to `smiles`
# those ring openings or closings.
# Args:
#  rings:
#  smiles:
function append_ring_info!(rings, smiles)
  for r in rings
    maybe_add_bond_type(r.bond_type, smiles)
    write(smiles, '0' + r.ring_number)
  end
end

atype_to_symbol = Dict{Int32, String}( dataset_pb2.BondTopology_AtomType.ATOM_C => "C",
                    dataset_pb2.BondTopology_AtomType.ATOM_N => "N",
                    dataset_pb2.BondTopology_AtomType.ATOM_NPOS => "N+",
                    dataset_pb2.BondTopology_AtomType.ATOM_O => "O",
                    dataset_pb2.BondTopology_AtomType.ATOM_ONEG => "O-",
                    dataset_pb2.BondTopology_AtomType.ATOM_F => "F",
                    dataset_pb2.BondTopology_AtomType.ATOM_H => "H" )
needs_square_brackets = Dict{Int32, Bool}( dataset_pb2.BondTopology_AtomType.ATOM_C => false,
                    dataset_pb2.BondTopology_AtomType.ATOM_N => false,
                    dataset_pb2.BondTopology_AtomType.ATOM_NPOS => true,
                    dataset_pb2.BondTopology_AtomType.ATOM_O => false,
                    dataset_pb2.BondTopology_AtomType.ATOM_ONEG => true,
                    dataset_pb2.BondTopology_AtomType.ATOM_F => false,
                    dataset_pb2.BondTopology_AtomType.ATOM_H => true )  # RDKit needs [H]

# Return a smiles token for atom type `atype` with `hcount` hydrogens.
function smiles_symbol(atype, hcount)::String
  asymbol = atype_to_symbol[atype]
  needs_square_brackets[atype] || return asymbol
  if hcount == 0
    return "[$(asymbol)]"
  elseif hcount == 1
    return "[$(asymbol)H]"
  else
    return "[$(asymbol)H$hcount]"
  end
end

# Return a smiles token for atom number `current_atom` in `bt` where
# `current_atom` has `hcount[current_atom]` implicit Hydrogens.
function smiles_symbol(bt::dataset_pb2.BondTopology,
                       current_atom,
                       hcount)::String
  if current_atom == 1  # No atom map symbol written.
    return smiles_symbol(bt.atoms[current_atom], hcount[current_atom])
  end
  asymbol = atype_to_symbol[bt.atoms[current_atom]]

  if hcount[current_atom] == 0
    hs = ""
  elseif hcount[current_atom] == 1
    hs = 'H'
  else
    hs = "H$(hcount[current_atom])"
  end
  "[$asymbol$hs:$(current_atom-1)]"
end

# If btype is anything other than a single bond, add the corresponding
# bond symbol to `smiles`.
function maybe_add_bond_type(btype, smiles::IOBuffer)
  btype < 2 && return
  if btype == 2
    write(smiles, '=')
  elseif btype == 3
    write(smiles, '#')
  end
end

# Recursive inner function for generating the smiles.
# As atoms are processed, the corresponding smiles parts are written to `smiles`.
function to_smiles(bt, ctable, atom_in_smiles, asymbol, incoming_connection, options::SmilesOptions, smiles::IOBuffer)
  current_atom = incoming_connection.atom
  atom_in_smiles[current_atom].complete = true
  write(smiles, asymbol[current_atom])
  append_ring_info!(atom_in_smiles[current_atom].ring_closures, smiles)
  append_ring_info!(atom_in_smiles[current_atom].ring_openings, smiles)
  @debug("From atom $(current_atom) branches $(atom_in_smiles[current_atom].branches)")
  nbranches = length(atom_in_smiles[current_atom].branches)
  for (i, c) in enumerate(atom_in_smiles[current_atom].branches)
    i == nbranches || write(smiles, '(')
    maybe_add_bond_type(c.btype, smiles)
    to_smiles(bt, ctable, atom_in_smiles, asymbol, c, options, smiles)
    i == nbranches || write(smiles, ')')
  end
end

# Return a Vector of atomic symbols for the atoms in `bt`, governed by
# `options`.
# Args:
#  bt: BondTopology
#  options: specifies options controlling smiles formation.
function atomic_symbols(bt::dataset_pb2.BondTopology,
                        options::SmilesOptions)::Vector{String}
  natoms = length(bt.atoms)
  hcount = compute_hydrogen_count(bt)
  if options.numbered_atoms
    return [smiles_symbol(bt, i, hcount) for i in 1:natoms]
  else
    return [smiles_symbol(bt.atoms[i], hcount[i]) for i in 1:natoms]
  end
end

# Return a Vector of atomic symbols for the atoms in `bt`, governed by
# `options`. In this case it is known that 3D smiles are being generated
# and so `geometry` contains those coordinates.
# Args:
#  bt: BondTopology
#  geometry:: contains coordinates of the atoms in `bt`.
#  options: specifies options controlling smiles formation.
function atomic_symbols(bt::dataset_pb2.BondTopology,
                        geometry::dataset_pb2.Geometry,
                        options::SmilesOptions)::Vector{String}
  natoms = length(bt.atoms)
  hcount = compute_hydrogen_count(bt)
  if options.numbered_atoms
    result = [smiles_symbol(bt, i, hcount) for i in 1:natoms]
  else
    result = [smiles_symbol(bt.atoms[i], hcount[i]) for i in 1:natoms]
  end

  if options.append_coordinates  # should always be true if this version is called.
    [result[i] = append_coordinates!(geometry.atom_positions[i], result[i]) for i in 1:natoms]
  end
  result
end

# Return a smiles for `bt`, given that the smiles of each atom has been precomputed and
# is in `asymbol`.
# Args:
#  bt: BondTopology
#  asymbol: smiles symbol of each atom in `bt`
#  options: options controlling smiles generation.
#           given that the smiles symbols are precomputed, this is not really needed here.
#           Leave for future use.
# Returns:
#  valid smiles.
function to_smiles(bt::dataset_pb2.BondTopology, 
                   asymbol::Vector{String},
                   options::SmilesOptions)::String
  natoms = length(bt.atoms)
  natoms == 1 && return asymbol[1]

  @debug(bt)
  @debug("Building smiles for $natoms atoms")

  ctable = ConnectionTable(natoms)
  if hasproperty(bt, :bonds)
    for bond in bt.bonds
      a1 = bond.atom_a + 1
      a2 = bond.atom_b + 1
      btype = smu_btype_to_number(bond.bond_type)
      @debug("Adding bond between atoms $a1 and $a2, type $btype")
      push!(ctable.connections[a1], Connection(a2, btype))
      push!(ctable.connections[a2], Connection(a1, btype))
    end
  end

  atom_in_smiles = build_smiles_order(ctable)
  smiles = IOBuffer()
  start_atom = 1
  while true
    to_smiles(bt, ctable, atom_in_smiles, asymbol, Connection(start_atom, 0), options, smiles)
    start_atom = findfirst(x->x.complete == false, atom_in_smiles)
    start_atom === nothing && break
    write(smiles, '.')
  end
  join(String(take!(smiles)))
end

# Generate a smiles for `bt` governed by `options`.
function to_smiles(bt::dataset_pb2.BondTopology, options::SmilesOptions)::String
  hasproperty(bt, :atoms) || return ""
  natoms = length(bt.atoms)
  natoms == 0 && return ""
  asymbols = atomic_symbols(bt, options)
  to_smiles(bt, asymbols, options)
end

# Generate a 3D smiles for `bt` governed by `options`, where the coordinates
# are in `geometry`.
function to_smiles(bt::dataset_pb2.BondTopology,
                   geometry::dataset_pb2.Geometry,
                   options::SmilesOptions)::String
  hasproperty(bt, :atoms) || return ""
  natoms = length(bt.atoms)
  natoms == 0 && return ""

  asymbols = atomic_symbols(bt, geometry, options)
  to_smiles(bt, asymbols, options)
end

# As the smiles is formed, transient information about each atom
# must be stored.
mutable struct AtomInSmiles
  atom_number::Integer
  ring_openings::Vector{RingClosure}
  ring_closures::Vector{RingClosure}
  branches::Vector{Connection}
  complete::Bool
  function AtomInSmiles()
    res = new()
    res.ring_openings = Vector{RingClosure}()
    res.ring_closures = Vector{RingClosure}()
    res.complete = false
    res
  end
end

# Return a Vector of AtomInSmiles ordered by a suitable path through the
# molecule.
function build_smiles_order(ctable::ConnectionTable)::Vector{AtomInSmiles}
  natoms = length(ctable.connections)
  res = [AtomInSmiles() for i in 1:natoms]
  visited = fill(0, natoms)
  first_atom = 1
  ring_number = RingNumbers()
  counter = 1
  while true
    atoms_found = build_smiles_order(ctable, -1, Connection(first_atom, 0), ring_number, visited, counter, res)
    first_atom = findfirst(x->x==false, visited)
    @debug("build_smiles_order next first_atom $(first_atom)")
    first_atom === nothing && break
    counter += atoms_found
    reset!(ring_number)
  end
  return res
end

# Recursively called function that incrementally fills in `result` which is a Vector
# of AtomInSmiles objects, indicting a path through `ctable` that can be used for
# a smiles.
# Args:
#  ctable:
#  previous_atom:
#  incoming_connection:
#  ring_numbers:
#  visited:
#  counter:
#  result:
# Returns:
#  The number of atoms we traverse - both here and recursively.
function build_smiles_order(ctable::ConnectionTable,
                previous_atom,
                incoming_connection,
                ring_numbers::RingNumbers,
                visited::Vector,
                counter::Integer,
                result)::Integer
  rc = 1
  current_atom = incoming_connection.atom
  visited[current_atom] = counter
  result[current_atom].atom_number = current_atom
  result[current_atom].branches = []
  @debug("build_smiles_order: atom $(current_atom) connections $(ctable.connections[current_atom])")
  for c in ctable.connections[current_atom]
    a = c.atom
    a == previous_atom && continue
    @debug("From $(current_atom) (counter $(counter)) to $(a) visited $(visited[a])")
    visited[a] > counter && continue  # Has been processed recusively.
    if visited[a] == 0
      push!(result[current_atom].branches, c)
      rc += build_smiles_order(ctable, current_atom, c, ring_numbers, visited, counter + 1, result)
    else  # Ring closure found
      ring_number = get_ring_number(ring_numbers)
      @debug("Got ring closure $(ring_number) btw $a and $current_atom")
      push!(result[a].ring_openings, RingClosure(ring_number, c.btype))
      push!(result[current_atom].ring_closures, RingClosure(ring_number, c.btype))
    end
  end

  rc
end

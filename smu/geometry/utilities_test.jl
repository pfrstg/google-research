#module UtilitiesTest

#export test_distance_between

using OffsetArrays
using ProtoBuf
using Test

using dataset_pb2

using SmuUtilities
using TopologyFromGeometry

export all_tests

function test_bohr_angstroms()
  @test isapprox(angstroms_to_bohr(SmuUtilities.bohr_to_angstroms(1.0)), 1.0)
end

function test_add_atom()
  bond_topology = BondTopology()
  add_atom_by_atomic_number!(bond_topology, 6)
  @test length(bond_topology.atoms) == 1
  @test bond_topology.atoms[1] == BondTopology_AtomType.ATOM_C

  add_atom_by_atomic_number!(bond_topology, 7)
  @test length(bond_topology.atoms) == 2
  @test bond_topology.atoms[2] == BondTopology_AtomType.ATOM_N

  add_atom_by_atomic_number!(bond_topology, 7, 1)
  @test length(bond_topology.atoms) == 3
  @test bond_topology.atoms[3] == BondTopology_AtomType.ATOM_NPOS

  # Do it systematically
  atype = [[1, 0], [6, 0], [7, 0], [7, 1], [8, 0], [8, 1], [9, 0]]
  expected = [
    BondTopology_AtomType.ATOM_H,
    BondTopology_AtomType.ATOM_C,
    BondTopology_AtomType.ATOM_N,
    BondTopology_AtomType.ATOM_NPOS,
    BondTopology_AtomType.ATOM_O,
    BondTopology_AtomType.ATOM_ONEG,
    BondTopology_AtomType.ATOM_F,
  ]
  @test length(atype) == length(expected)

  bt = BondTopology()
  for i in zip(atype, expected)
    add_atom_by_atomic_number!(bt, i[1][1], i[1][2])
    @test bt.atoms[end] == i[2]
  end
  @test length(bt.atoms) == length(expected)
  true
end

function test_add_single_bond()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 1, bond_topology)
  @test length(bond_topology.bonds) == 1
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_SINGLE
  true
end

function test_add_double_bond()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 2, bond_topology)
  @test length(bond_topology.bonds) == 1
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_DOUBLE
  true
end

function test_add_triple_bond()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 3, bond_topology)
  @test length(bond_topology.bonds) == 1
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_TRIPLE
  true
end

function test_add_multiple_bonds()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 1, bond_topology)
  add_bond!(1, 2, 1, bond_topology)
  @test length(bond_topology.bonds) == 2
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_SINGLE
  @test bond_topology.bonds[2].bond_type == BondTopology_BondType.BOND_SINGLE
  true
end

function test_distance_between()
  geom = Geometry()
  atoms = Vector{Geometry_AtomPos}()
  a1 = Geometry_AtomPos(x=0.0, y=0.0, z=0.0)
  push!(atoms, a1)
  a2 = Geometry_AtomPos(x=0.0, y=0.0, z=angstroms_to_bohr(1.0))
  push!(atoms, a2)
  setproperty!(geom, :atom_positions, atoms)
  @test isapprox(distance_between_atoms(geom, 1, 2), 1.0)
  setproperty!(geom.atom_positions[2], :x, angstroms_to_bohr(1.0))
  setproperty!(geom.atom_positions[2], :y, angstroms_to_bohr(1.0))
  setproperty!(geom.atom_positions[2], :z, angstroms_to_bohr(1.0))
  @test isapprox(distance_between_atoms(geom, 1, 2), sqrt(3.0))

  distances = SmuUtilities.distances(geom)
  @test distances[1,1] == 0.0
  @test distances[2,2] == 0.0
  @test isapprox(distances[1,2], sqrt(3.0))
  @test isapprox(distances[2,1], sqrt(3.0))
  true
end

function test_bonded()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  connection_matrix = SmuUtilities.bonded(bond_topology)
  @test length(atoms) == size(connection_matrix, 1)
  @test connection_matrix[1,2] > 0
  @test connection_matrix[2,3] > 0
  @test count(!iszero, connection_matrix) == 4  # Bonds are placed twice.
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=0, atom_b=2, bond_type=BondTopology_BondType.BOND_DOUBLE))
  connection_matrix = SmuUtilities.bonded(bond_topology)
  @test count(!iszero, connection_matrix) == 6  # Bonds are placed twice.
  @test connection_matrix[1,2] == connection_matrix[2,3]
  @test connection_matrix[1,2] != connection_matrix[3,1]
  true
end

function test_connections()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  connections = SmuUtilities.connections(SmuUtilities.bonded(bond_topology))
  @test length(connections) == length(atoms)
  @test length(connections[1]) == 1
  @test length(connections[2]) == 2
  @test length(connections[3]) == 1
  @test connections[1][1] == 2
  @test connections[2][1] == 1
  @test connections[2][2] == 3
  @test connections[3][1] == 2
  true
end

function test_canonical_bond_topology()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=2, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=0, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  SmuUtilities.canonical_bond_topology!(bond_topology)

  expected = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(expected, :atoms, atoms)
  setproperty!(expected, :bonds, bonds)
  @test bond_topology == expected
  @test SmuUtilities.same_bond_topology(bond_topology, expected)

  true
end

function test_single_fragment1()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)  # atom 0
  push!(atoms, BondTopology_AtomType.ATOM_C)  # atom 1
  push!(atoms, BondTopology_AtomType.ATOM_C)  # atom 2
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  @test !SmuUtilities.is_single_fragment(bond_topology)

  push!(bond_topology.bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test SmuUtilities.is_single_fragment(bond_topology)
  
  push!(bond_topology.atoms, BondTopology_AtomType.ATOM_C)  # atom 3
  @test !SmuUtilities.is_single_fragment(bond_topology)

  push!(bond_topology.atoms, BondTopology_AtomType.ATOM_C)  # atom 4
  @test !SmuUtilities.is_single_fragment(bond_topology)
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=3, atom_b=4, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test !SmuUtilities.is_single_fragment(bond_topology)

  push!(bond_topology.bonds, BondTopology_Bond(atom_a=2, atom_b=3, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test SmuUtilities.is_single_fragment(bond_topology)
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=0, atom_b=4, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test SmuUtilities.is_single_fragment(bond_topology)

  true
end

function test_append_bond_type()
  buffer = IOBuffer()
  expected = OffsetArray(['=', '#'], 2:3)
  SmuUtilities.maybe_add_bond_type(0, buffer)
  @test isempty(String(take!(buffer)))

  SmuUtilities.maybe_add_bond_type(1, buffer)
  @test isempty(String(take!(buffer)))

  for b in 2:3
    SmuUtilities.maybe_add_bond_type(b, buffer)
    smiles = String(take!(buffer))
    @test length(smiles) == 1
    @test smiles[1] == expected[b]
  end

  true
end

function test_append_ring_info()
  for ring_number in 1:9
    for btype in 1:3
      buffer = IOBuffer()
      SmuUtilities.append_ring_info!([SmuUtilities.RingClosure(ring_number, btype)], buffer)
      smiles = String(take!(buffer))
      if btype == 1
        @test length(smiles) == 1
      else
        @test length(smiles) == 2
        if btype == 2
          @test smiles[1] == '='
        else
          @test smiles[1] == '#'
        end
      end
      @test smiles[end] == '0' + ring_number
    end
  end

  true
end

function test_symbol_With_atom_map()
  bt = dataset_pb2.BondTopology()
  natoms = 99
  atoms = [dataset_pb2.BondTopology_AtomType.ATOM_C for i in 1:natoms]
  setproperty!(bt, :atoms, atoms)

  hcount = zeros(Integer, natoms)

  for i in 1:99
    s = SmuUtilities.smiles_symbol(bt, i, hcount)
    if i == 1
      @test s == "C"
    else
      @test s == "[C:$(i-1)]"
    end
  end

  true
end

function test_single_atom_smiles()
  atom_smiles = Dict(dataset_pb2.BondTopology_AtomType.ATOM_C => "C",
                     dataset_pb2.BondTopology_AtomType.ATOM_N => "N",
                     dataset_pb2.BondTopology_AtomType.ATOM_NPOS => "[N+H4]",
                     dataset_pb2.BondTopology_AtomType.ATOM_O => "O",
                     dataset_pb2.BondTopology_AtomType.ATOM_ONEG => "[O-H]",
                     dataset_pb2.BondTopology_AtomType.ATOM_F => "F")
  for (atype, expected_smiles) in atom_smiles
    bt = dataset_pb2.BondTopology()
    setproperty!(bt, :atoms, [])
    push!(bt.atoms, atype)
    smi = to_smiles(bt, SmilesOptions())
    @test smi == expected_smiles
  end

  true
end

function test_two_atom_smiles()
  bt = dataset_pb2.BondTopology()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  setproperty!(bt, :atoms, [carbon, carbon])
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  bond = dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond)
  setproperty!(bt, :bonds, [bond])
  smi = to_smiles(bt, SmilesOptions())
  @test smi == "CC"

  true
end

function test_two_atom_smiles_numbered()
  bt = dataset_pb2.BondTopology()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  setproperty!(bt, :atoms, [carbon, carbon])
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  bond = dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond)
  setproperty!(bt, :bonds, [bond])
  smi = to_smiles(bt, SmilesOptions(true))
  @test smi == "C[CH3:1]"

  true
end

function test_long_chain()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  for natoms in 1:8
    bt = dataset_pb2.BondTopology()
    setproperty!(bt, :atoms, [carbon for i in 1:natoms])
    bonds = []
    for i in 2:natoms
      push!(bonds, dataset_pb2.BondTopology_Bond(atom_a=(i-1-1), atom_b=(i - 1), bond_type=single_bond))
    end
    setproperty!(bt, :bonds, bonds)
    smi = to_smiles(bt, SmilesOptions())
    @test smi == "C" ^ natoms
  end

  true
end

function test_neopentane()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  bt = dataset_pb2.BondTopology()
  natoms = 5
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [dataset_pb2.BondTopology_Bond(atom_a = 0, atom_b = (i-1), bond_type=single_bond) for i in 2:natoms])
  smi = to_smiles(bt, SmilesOptions())
  @test smi == "C(C)(C)(C)C"
  true
end

function test_rings()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  for rsize in 3:8
    bt = dataset_pb2.BondTopology()
    setproperty!(bt, :atoms, [carbon for i in 1:rsize])
    setproperty!(bt, :bonds, [dataset_pb2.BondTopology_Bond(atom_a=(i-1), atom_b=(i%rsize), bond_type=single_bond) for i in 1:rsize])
    smi = to_smiles(bt, SmilesOptions())
    @test smi == "C1" * "C"^(rsize-2) * "C1"
  end

  true
end

function test_tert_butylcoclopropane()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  natoms = 7
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=4, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=5, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=5, atom_b=6, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=6, bond_type=single_bond))

  smi = to_smiles(bt, SmilesOptions())
  @test smi == "CC(C)(C)C1CC1"

  true
end

function test_not_sure_name1()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  natoms = 5
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=4, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=4, bond_type=single_bond))
  smi = to_smiles(bt, SmilesOptions())
  @test smi == "C12CC1CC2"

  true
end

function test_quadricyclane()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  natoms = 7
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=5, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=3, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=6, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=4, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=5, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=6, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=5, atom_b=6, bond_type=single_bond))
  smi = to_smiles(bt, SmilesOptions())
  @test smi == "C2C1C3C1C4C2C34"

  true
end

function test_multiple_bonds()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  double_bond = dataset_pb2.BondTopology_BondType.BOND_DOUBLE
  triple_bond = dataset_pb2.BondTopology_BondType.BOND_TRIPLE
  natoms = 5
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=double_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=4, bond_type=triple_bond))
  smi = to_smiles(bt, SmilesOptions())
  @test smi == "CC=CC#C"
  true
end

function test_disconnected()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  double_bond = dataset_pb2.BondTopology_BondType.BOND_DOUBLE
  triple_bond = dataset_pb2.BondTopology_BondType.BOND_TRIPLE
  natoms = 8
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=5, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=5, atom_b=6, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=5, atom_b=7, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=6, atom_b=7, bond_type=single_bond))

  smi = to_smiles(bt, SmilesOptions())
  @test smi == "CC1CC1.CC1CC1"

  true
end

function test_cubane()
  carbon = dataset_pb2.BondTopology_AtomType.ATOM_C
  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  double_bond = dataset_pb2.BondTopology_BondType.BOND_DOUBLE
  triple_bond = dataset_pb2.BondTopology_BondType.BOND_TRIPLE
  natoms = 8
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [carbon for i in 1:natoms])
  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=7, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=6, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=5, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=4, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=5, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=7, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=5, atom_b=6, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=6, atom_b=7, bond_type=single_bond))

  smi = to_smiles(bt, SmilesOptions())
  @test smi == "C14C3C2C1C5C2C3C45"

  true
end

function test_all_atom_types()
  natoms = 8
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [])
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATYPE_F)
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATYPE_C)
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATYPE_O)
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATYPE_N)
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATYPE_ONEG)
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATYPE_NPOS)

  setproperty!(bt, :bonds, [])

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=1, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=4, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=5, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=6, bond_type=single_bond))

  smi = to_smiles(bt, SmilesOptions())
  @test smi == "FCONC([O-])[NH2+]"

  true
end

function test_observed_problem1()
  natoms = 13
  bt = dataset_pb2.BondTopology()
  setproperty!(bt, :atoms, [])
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 0
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 1
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 2
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 3
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_N)  # 4
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_O)  # 5
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_F)  # 6
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 7
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 8
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 9
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 10
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 11
  push!(bt.atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 12

  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE
  double_bond = dataset_pb2.BondTopology_BondType.BOND_DOUBLE

  setproperty!(bt, :bonds, [])
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=4, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=6, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=0, atom_b=7, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=2, bond_type=double_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=5, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=1, atom_b=8, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=3, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=2, atom_b=9, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=4, bond_type=single_bond))
  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=3, atom_b=10, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=4, atom_b=11, bond_type=single_bond))

  push!(bt.bonds, dataset_pb2.BondTopology_Bond(atom_a=5, atom_b=12, bond_type=single_bond))

  smi = to_smiles(bt, SmilesOptions())

  @test smi == "C1(C(C(=C(O[H])[H])[H])(N1[H])[H])(F)[H]"

  true
end

# Molecule is propane
#
#  H2                H9
#    \       H5     /
#     \      |     /
#  H1 - C0 - C4 - C7 - H8
#      /     |     \
#     /      H6     \
#    H3              H10

function test_remove_atom_from_bond_topology()
  bt = dataset_pb2.BondTopology()
  atoms = []
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 0
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 1
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 2
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 3
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 4
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 5
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 6
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_C)  # 7
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 8
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 9
  push!(atoms, dataset_pb2.BondTopology_AtomType.ATOM_H)  # 10
  setproperty!(bt, :atoms, atoms)

  single_bond = dataset_pb2.BondTopology_BondType.BOND_SINGLE

  bonds = []

  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 0, atom_b = 1, bond_type=single_bond))
  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 0, atom_b = 2, bond_type=single_bond))
  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 0, atom_b = 3, bond_type=single_bond))

  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 0, atom_b = 4, bond_type=single_bond))

  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 4, atom_b = 5, bond_type=single_bond))
  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 4, atom_b = 6, bond_type=single_bond))

  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 4, atom_b = 7, bond_type=single_bond))

  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 7, atom_b = 8, bond_type=single_bond))
  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 7, atom_b = 9, bond_type=single_bond))
  push!(bonds, dataset_pb2.BondTopology_Bond(atom_a = 7, atom_b = 10, bond_type=single_bond))

  setproperty!(bt, :bonds, bonds)
  @test length(bt.atoms) == 11
  @test length(bt.bonds) == 10

  canonical_bond_topology!(bt)
  remove_atoms!(x->x == dataset_pb2.BondTopology_AtomType.ATOM_N, bt)  # Should be a no-op
  @test length(bt.atoms) == 11
  @test length(bt.bonds) == 10

  remove_atoms!(x->x == dataset_pb2.BondTopology_AtomType.ATOM_H, bt)
  @test length(bt.atoms) == 3
  @test length(bt.bonds) == 2
  @test are_bonded(bt, 0, 1)
  @test are_bonded(bt, 1, 2)

  true
end

function all_tests()
  @test test_add_atom()
  @test test_add_single_bond()
  @test test_add_double_bond()
  @test test_add_triple_bond()
  @test test_distance_between()
  @test test_bonded()
  @test test_connections()
  @test test_canonical_bond_topology()
  @test test_single_fragment1()
  @test test_single_atom_smiles()
  @test test_two_atom_smiles()
  @test test_two_atom_smiles_numbered()
  @test test_long_chain()
  @test test_neopentane()
  @test test_rings()
  @test test_append_bond_type()
  @test test_append_ring_info()
  @test test_symbol_With_atom_map()
  @test test_tert_butylcoclopropane()
  @test test_not_sure_name1()
  @test test_quadricyclane()
  @test test_multiple_bonds()
  @test test_disconnected()
  @test test_cubane()
  @test test_observed_problem1()
  @test test_remove_atom_from_bond_topology()
  print("All tests complete\n")
end


#end  # module UtilitiesTest

# Given a Conformer proto, generate 3d query protos.

using ArgMacros
using Combinatorics
using DataStructures
using ProtoBuf

using dataset_pb2
using SmuUtilities

is_hydrogen = x->x == dataset_pb2.BondTopology_AtomType.ATOM_H

function compute_tolerance(d::Float32, tolerance::Float32)::Tuple{Float32, Float32}
  minval = max(d - tolerance * d, 0.0)
  maxval = d + tolerance * d
  (round(minval, digits=2), round(maxval, digits=2))
end

function make_3d_query(conformer_id,
                       bt::dataset_pb2.BondTopology,
                       geometry::dataset_pb2.Geometry,
                       tolerance::Number,
                       output::IO)
  
  dists = distances(geometry)
  smarts = Vector{String}()
  heavy_atoms = []
  ecount = DefaultDict{Integer, Integer}(0)
  for (i, atom) in enumerate(bt.atoms)
    atom == dataset_pb2.BondTopology_AtomType.ATOM_H && continue
    atnum = smu_atype_to_atomic_number(atom)
    push!(smarts, "[#$atnum]")
    push!(heavy_atoms, i)
    ecount[atnum] += 1
  end
  println(output, "query {")
  println(output, "  smarts: \"$(join(smarts, '.'))\"")
  println(output, "  comment: \"$(conformer_id)\"")
  for (k, v) in ecount
    println(output, "  elements_needed {")
    println(output, "    atomic_number: $k")
    println(output, "    hits_needed: $v")
    println(output, "  }")
  end
  println(output, "  geometric_constraints {")
  for (i, j) in Combinatorics.combinations(heavy_atoms, 2)
    d = dists[i, j]
    (minval,maxval) = compute_tolerance(d, tolerance)
    println(output, "     distances {")
    println(output, "       a1: $(i - 1)")
    println(output, "       a2: $(j - 1)")
    println(output, "       range {")
    println(output, "         min: $minval")
    println(output, "         max: $maxval")
    println(output, "         one_time_scaling_factor: 1.88973")
    println(output, "       }")
    println(output, "     }")
  end
  println(output, "  }")
  println(output, "}")
end

function main()
  @inlinearguments begin
    @argumentrequired String input_fname "-i" "--input"
    @argumentrequired String output_fname "-o" "--output"
    @argumentoptional Float32 tolerance "-t" "--tolerance"
  end
  input = open(input_fname, "r")
  output = open(output_fname, "w")
  if tolerance === nothing
    tolerance = 0.1
  end

  proto = readproto(input, dataset_pb2.Conformer())
  hasproperty(proto, :conformer_id) || return
  hasproperty(proto, :bond_topologies) || return
  hasproperty(proto, :optimized_geometry) || return
  make_3d_query(proto.conformer_id, proto.bond_topologies[1],
                proto.optimized_geometry, tolerance, output)
end

main()

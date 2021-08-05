# Discern connectivity from the coordinates in a Geometry

using ArgMacros
using Glob
using Logging
using ProtoBuf
using TFRecord

append!(LOAD_PATH, ["smu/jlout/", "smu/topology_from_geometry/", "smu/geometry"])

using BondLengthDistributions
using dataset_pb2
using TopologyFromGeometry

function get_topology_from_geometry(bond_length_distributions::AllBondLengthDistributions,
                conformer_id::Integer,
                bond_topology::BondTopology,
                geometry::Geometry,
                energy::Number,
                output_stream::IOStream)::Bool
  matching_parameters = MatchingParameters()
  @debug("in main $(length(bond_topology.atoms)) atoms and $(length(bond_topology.bonds)) bonds")
  result = bond_topologies_from_geom(bond_length_distributions, bond_topology,
                                     geometry, matching_parameters)
  hasproperty(result, :bond_topology) || return false
  number_bond_topologies = length(result.bond_topology)
  number_bond_topologies > 0 || return
  sep = ','
  energy = round(energy, digits=3)
  to_write = "$(conformer_id)$sep$(number_bond_topologies)$sep$(energy)"
  for bt in result.bond_topology
    score = round(bt.score, digits=3)
    to_write *= "$sep$(bt.smiles)$sep$(score)$sep"
    to_write *= bt.is_starting_topology ? "T" : "F"
  end
  println(output_stream, to_write)
  true
end

"""Return true if `conformer` has all the properties needed for calculations"""
function has_required_properties(conformer::dataset_pb2.Conformer)::Bool
  if ! hasproperty(conformer, :optimized_geometry)
    @warn("No optimized geometry, skipping")
    return false
  end
  if ! hasproperty(conformer, :bond_topologies)
    @warn("No bond_topology, skipping")
    return false
  end
  if !hasproperty(conformer, :properties)
    @warn("No properties in COnformer")
    return false
  end
  if !hasproperty(conformer.properties, :single_point_energy_pbe0d3_6_311gd)
    @warn("No single_point_energy_pbe0d3_6_311gd in propertoes")
    return false
  end

  true
end

"""Examine Conformer file(s) and determine the BondTopology of the
molecules, based on a set of observed bond length distributions.
Input is via either the --input_fname or --input_glob options.
"""
function main()
  @inlinearguments begin
    @argumentoptional String input_fname "-i" "--input"
    @argumentoptional String input_glob "-g" "--input_glob"
    @argumentrequired String output_fname "-o" "--output"
    @argumentrequired String bonds "-b" "--bonds"
    @argumentoptional Int64 nprocess "-N" "--nprocess"
    @argumentflag exclude_non_bonded "-x" "--xcldnonbond"
    @argumentflag ignore_hydrogens "-h" "--xhydrogen"
    @argumentflag debug "-debug" "--debug"
    @argumentflag verbose "-v" "--verbose"
  end
  if debug
    logger=Logging.SimpleLogger(stderr,Logging.Debug)
    global_logger(logger)
  end

  @info("input $(input_fname) bonds $(bonds) output $(output_fname) $nprocess")
  flush(stdout)
  nprocess === nothing && (nprocess = typemax(Int64))

  bond_length_distributions = AllBondLengthDistributions()
  exclude_non_bonded !== nothing && (bond_length_distributions.include_non_bonded = false)

  add_from_files!(bonds, bond_length_distributions) || @error("Cannot build bond length distribution $bonds")

  input_files = []
  if input_fname != nothing
    push!(input_files, input_fname)
  elseif input_glob != nothing
    input_files = filter(x->occursin(input_glob, x), readdir(dirname(input_glob)))
    if length(input_files) == 0
      @error("No files match $(input_glob)")
      return
    end
  else
    @error("Must specify an input source")
    return
  end

  molecules_read = 0
  molecules_processed = 0
  calculations_performed = 0

  open(output_fname, "w") do output_stream
    println(output_stream, "ConformerId,Energy,NBT,Score,Smiles")
    for conformer in TFRecord.read(input_files, record_type=Conformer)
      molecules_read += 1
      has_required_properties(conformer) || continue

      @info("processing $(conformer.conformer_id)")
      if get_topology_from_geometry(bond_length_distributions,
                                    conformer.conformer_id,
                                    conformer.bond_topologies[1],
                                    conformer.optimized_geometry,
                                    conformer.properties.single_point_energy_pbe0d3_6_311gd.value,
                                    output_stream)
        molecules_processed += 1
      end
      molecules_read > nprocess && break
    end
  end
  verbose && println("Read $(molecules_read) molecules, processed $(molecules_processed)")
  0
end


main()

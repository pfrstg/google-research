import sys

from absl import app
from absl import flags
from absl import logging

import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions
from apache_beam.io.tfrecordio import ReadFromTFRecord
import tensorflow as tf

from smu import dataset_pb2
from smu.geometry import bond_length_distribution
from smu.geometry import smu_molecule

import topology_from_geom

FLAGS = flags.FLAGS

flags.DEFINE_string("input", None, "TFDataRecord file containg Conformer protos")
flags.DEFINE_string("bonds", None, "File name stem for bond length distributions")
flags.DEFINE_string("output", None, "Output file")
flags.DEFINE_boolean("xnonbond", False, "Exclude non bonded interactions")
flags.DEFINE_boolean("beam", True, "Process with Beam")

def topology_matches_summary(topology_matches:dataset_pb2.TopologyMatches) -> str:
  """Create an output string containing summary data from `topology_matches`.

  Args:
    topology_matches:
  Returns:
    string suitable for writing.
  """
  energy = round(topology_matches.single_point_energy_pbe0d3_6_311gd, 2)
  result = f"{topology_matches.conformer_id},{len(topology_matches.bond_topology)},{energy}"

  for bt in topology_matches.bond_topology:
    result += f",{bt.score:.3f},{bt.smiles}"
    if bt.is_starting_topology:
      result += ",T"
    else:
      result += ",F"

  return result


def non_beam(bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions, input_fname: str,
             output_fname: str) -> None:
  """
  Args:
    bond_lengths:
    input_fname:
    output_fname:
  Returns:
  """
  items_read = 0
  calculations_performed = 0
  matching_parameters = smu_molecule.MatchingParameters()
  raw_dataset = tf.data.TFRecordDataset(input_fname)
  with open(output_fname, "w") as output:
    for raw_record in raw_dataset:
      conformer = dataset_pb2.Conformer()
      conformer.ParseFromString(raw_record.numpy())
      items_read += 1
      if conformer.fate != dataset_pb2.Conformer.FATE_SUCCESS:
        continue
      topology_matches = topology_from_geom.bond_topologies_from_geom(bond_lengths,
                     conformer.conformer_id,
                     conformer.bond_topologies[0],
                     conformer.optimized_geometry,
                     conformer.properties.single_point_energy_pbe0d3_6_311gd,
                     matching_parameters)
      print(topology_matches_summary(topology_matches), file=output)
      calculations_performed ++ 1

  print(f"Read {items_read} items, {calculations_performed} calculations performed", file=sys.stderr)

class SummaryData(beam.DoFn):
  """Given BondTopologies as input, yield summary data"""
  def process(self, topology_matches:dataset_pb2.TopologyMatches):
    yield topology_matches_summary(topology_matches)


def ReadConFormer(bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions, input: str,
                  output: str):
  """
  Args:
  Returns:
  """

  class GetAtoms(beam.DoFn):

    def process(self, item):
      yield item.optimized_geometry.atom_positions[0].x

# options = PipelineOptions(direct_num_workers=6, direct_running_mode='multi_processing')
  options = PipelineOptions()

  with beam.Pipeline(options=options) as p:
    protos = (p | beam.io.tfrecordio.ReadFromTFRecord(
        input, coder=beam.coders.ProtoCoder(dataset_pb2.Conformer)) |
              beam.ParDo(topology_from_geom.TopologyFromGeom(bond_lengths)) |
              beam.ParDo(SummaryData()) |
              beam.io.textio.WriteToText(output))

    return protos


def topology_from_geometry_main(unused_argv):
  del unused_argv

  bond_lengths = bond_length_distribution.AllAtomPairLengthDistributions()
  bond_lengths.add_from_files(FLAGS.bonds, 0.0, FLAGS.xnonbond)

  if FLAGS.beam:
    protos = ReadConFormer(bond_lengths, FLAGS.input, FLAGS.output)
    print(protos)
  else:
    non_beam(bond_lengths, FLAGS.input, FLAGS.output)
    


if __name__ == "__main__":
  flags.mark_flag_as_required("input")
  flags.mark_flag_as_required("bonds")
  flags.mark_flag_as_required("output")
  app.run(topology_from_geometry_main)

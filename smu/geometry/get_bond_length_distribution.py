"""Extracts bond length distributions from existing Conformer protos."""

import itertools

import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions
import numpy as np
import tensorflow as tf

from absl import app
from absl import flags
from absl import logging

from smu import dataset_pb2

import bond_lengths

FLAGS = flags.FLAGS

flags.DEFINE_string("input", None, "TFDataRecord file containg Conformer protos")
flags.DEFINE_string("output", None, "Output file")
flags.DEFINE_boolean("beam", True, "Use Beam")

def get_bond_length_distribution_nobeam(input_fname: str, output_fname: str):
  """Generate bond length distributions from the Conformer data in `input_fname`.
  Args:
    input_fname: An existing TFRecord file containing Conformer protos.
    output_fname: An output file that will be created that contains
      all bond length distributions - all bond types, all atom types.
      Requires post-processing to generate bond length distribution files.
  """
  data: Dict[Tuple[int, int, int], np.array] = {}
  atomic_numbers = [1, 6, 7, 8, 9]
  for (z1, z2) in itertools.combinations_with_replacement(atomic_numbers, 2):
    for btype in range(0, 4):
      key = (min(z1, z2), btype, max(z1, z2))
      data[key] = np.zeros(200001, dtype=np.int32)

  items_read = 0
  nprocessed = 0
  raw_dataset = tf.data.TFRecordDataset(input_fname)
  for raw_record in raw_dataset:
    conformer = dataset_pb2.Conformer()
    conformer.ParseFromString(raw_record.numpy())
    items_read += 1
    if conformer.fate != dataset_pb2.Conformer.FATE_SUCCESS:
      continue
    bond_lengths.bond_lengths(conformer, data)
    nprocessed += 1
#   if items_read > 20000:
#     break

  bond_lengths.write_bond_length_distributions(data, output_fname)
  print(f"Read {items_read} processed {nprocessed} FATE_SUCCESS conformers")


class BondDistToString(beam.DoFn):
  """Generate the dot separated form of a bond length distribution component """
  def process(self, bond_dist):
    key, value = bond_dist
    print(f"BondDistToString: key {key} value {value}")
    yield f"{key[0]}.{key[1]}.{key[2]}.{key[3]}.{value}"

class GroupBondTypes(beam.DoFn):
  def process(self, bond_dist):
    key, value = bond_dist
    print(f"GroupBondTypes: key #{key} value {value}")
    yield (key[0], key[1], key[2]), (key[3], value)

def get_bond_length_distribution_beam(input_fname: str, output_fname: str):
  """Generate bond length distibutions.

  Args:
    input_fname: An existing TFRecord file containing Conformer protos.
    output_fname: An output file that will be created that contains
      all bond length distributions - all bond types, all atom types.
      Requires post-processing to generate bond length distribution files.
  """
  print("Reading from {input_fname} output to {output_fname}")
  with beam.Pipeline(options=PipelineOptions()) as p:
    protos = (p
      | beam.io.tfrecordio.ReadFromTFRecord(input_fname, coder=beam.coders.ProtoCoder(dataset_pb2.Conformer().__class__))
      | beam.ParDo(bond_lengths.GetBondLengthDistribution())
      | beam.CombinePerKey(sum)
#     | beam.ParDo(GroupBondTypes())
#     | beam.GroupByKey()
      | beam.ParDo(BondDistToString())
      | beam.io.WriteToText(output_fname)
    )
    print(protos)

def get_bond_length_distribution(unused_argv):
  """Scan Conformer protos to extract bond length distributions."""
  del unused_argv
  print(f"Beam? {FLAGS.beam}")
  if FLAGS.beam:
    get_bond_length_distribution_beam(FLAGS.input, FLAGS.output)
  else:
    get_bond_length_distribution_nobeam(FLAGS.input, FLAGS.output)


if __name__ == "__main__":
  app.run(get_bond_length_distribution)

"""Create a subset of the Conformer data containing
  First BondTopology
  Optimized Geometry
  Single Point Energy
  and where the FATE is FATE_SUCCESS or FATE_DISASSOCIATED
"""

from typing import Optional

from absl import app
from absl import flags

import tensorflow as tf

from smu import dataset_pb2
from smu.geometry import utilities

FLAGS = flags.FLAGS

flags.DEFINE_string("input", None, "TFDataRecord file containg Conformer protos")
flags.DEFINE_string("output", None, "Riegeli output file")

def create_subset(conformer: dataset_pb2.Conformer) -> Optional[dataset_pb2.Conformer]:
  """Create a subset of `conformer`."""
  if not conformer.fate in [dataset_pb2.Conformer.FATE_SUCCESS, dataset_pb2.Conformer.FATE_DISASSOCIATED]:
    return None
  if conformer.duplicated_by > 0:
    return None

  subset = dataset_pb2.Conformer()
  subset.conformer_id = conformer.conformer_id
  subset.optimized_geometry.CopyFrom(utilities.geom_to_angstroms(conformer.optimized_geometry))
  subset.properties.single_point_energy_pbe0d3_6_311gd.CopyFrom(conformer.properties.single_point_energy_pbe0d3_6_311gd)
  subset.bond_topologies.extend(conformer.bond_topologies)
  subset.fate = conformer.fate
  return subset


def create_and_write_subset(conformer: dataset_pb2.Conformer,
              output) -> bool:
  """Create a subset of `conformer` and write to `output`.
  """
  subset = create_subset(conformer)
  if subset is None:
    return False

  output.write(subset.SerializeToString())
  return True

def make_subset_inner(input_fname: str, output_fname: str):
  """Reads DataSet protos from `input_fname` and writes a subset to `output_fname`.

  Args:
    input_fname: TFRecordDataset of dataset protos.
    output_fname: a TFrecord file we create.
  """
  items_read = 0
  items_written = 0
  raw_dataset = tf.data.TFRecordDataset(input_fname)
  with tf.io.TFRecordWriter(path=output_fname) as output:
    for raw_record in raw_dataset:
      conformer = dataset_pb2.Conformer()
      conformer.ParseFromString(raw_record.numpy())
      items_read += 1
      if create_and_write_subset(conformer, output):
        items_written += 1

  print(f"Read {items_read} write {items_written}")

def make_subset(unused_argv):
  """Creates a subset of the SMU data in TFDataRecord form"""
  make_subset_inner(FLAGS.input, FLAGS.output)

if __name__ == "__main__":
  flags.mark_flag_as_required("input")
  flags.mark_flag_as_required("output")
  app.run(make_subset)

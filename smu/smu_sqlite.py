# coding=utf-8
# Copyright 2022 The Google Research Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Copyright 2022 The Google Research Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Interface to a SQLite DB file for SMU data.

Provides a simpler interface than SQL to create and access the SMU data in an
SQLite database.

The majority of the data is stored as a blob, with just the bond topology id and
smiles string pulled out as fields.
"""
import datetime
import os
import sqlite3

from absl import logging
from rdkit import Chem

from smu import dataset_pb2
from smu.parser import smu_utils_lib
from smu.geometry import smu_molecule
from smu.geometry import topology_from_geom
import snappy

_CONFORMER_TABLE_NAME = 'conformer'
_BTID_TABLE_NAME = 'btid'
_SMILES_TABLE_NAME = 'smiles'


class ReadOnlyError(Exception):
  pass


class SMUSQLite:
  """Provides an interface for SMU data to a SQLite DB file.

  The class hides away all the SQL fun with just Conformer protobuf visible in
  the interface.

  Internal details about the tables:
  There are 3 separate tables
  * conformer: Is the primary table which has columns
      * cid: integer conformer id (unique)
      * conformer: blob wire format proto of a conformer proto
  * btid: Used for lookups by bond topology id which has columns
      * btid: integer bond topology id (not unique)
      * cid: integer conformer id (not unique)
  * smiles: Used to map smiles to bond topology ids with columns
      * smiles: text canonical smiles string (unique)
      * btid: integer bond topology id
    Note that if multiple smiles strings are associated with the same bond
    toplogy id, the first one provided will be silently kept.
  """

  def __init__(self, filename, mode='r'):
    """Creates SMUSQLite.

    Args:
      filename: database file, must be on local filesystem
      mode: 'c' (create, deletes existing), 'w' (writable), 'r' (read only)

    Raises:
      FileNotFoundError: if 'r' and file does not exist
    """
    if mode == 'c':
      if os.path.exists(filename):
        os.remove(filename)
      self._read_only = False
      self._conn = sqlite3.connect(filename)
      self._maybe_init_db()
    elif mode == 'w':
      self._read_only = False
      self._conn = sqlite3.connect(filename)
      self._maybe_init_db()
    elif mode == 'r':
      if not os.path.exists(filename):
        raise FileNotFoundError(filename)
      self._conn = sqlite3.connect(filename)
      self._read_only = True
    else:
      raise ValueError('Mode must be c, r, or w')

    self._conn = sqlite3.connect(filename)

  def _maybe_init_db(self):
    """Create the table and indices if they do not exist."""
    make_table = (f'CREATE TABLE IF NOT EXISTS {_CONFORMER_TABLE_NAME} '
                  '(cid INTEGER PRIMARY KEY, '
                  'exp_stoich STRING, '
                  'conformer BLOB)')
    self._conn.execute(make_table)
    self._conn.execute(f'CREATE UNIQUE INDEX IF NOT EXISTS '
                       f'idx_cid ON {_CONFORMER_TABLE_NAME} (cid)')
    self._conn.execute(f'CREATE INDEX IF NOT EXISTS '
                       f'idx_exp_stoich ON {_CONFORMER_TABLE_NAME} '
                       '(exp_stoich)')
    self._conn.execute(f'CREATE TABLE IF NOT EXISTS {_BTID_TABLE_NAME} '
                       '(btid INTEGER, cid INTEGER)')
    self._conn.execute(f'CREATE INDEX IF NOT EXISTS '
                       f'idx_btid ON {_BTID_TABLE_NAME} (btid)')
    self._conn.execute(f'CREATE TABLE IF NOT EXISTS {_SMILES_TABLE_NAME} '
                       '(smiles TEXT, btid INTEGER)')
    self._conn.execute(f'CREATE UNIQUE INDEX IF NOT EXISTS '
                       f'idx_smiles ON {_SMILES_TABLE_NAME} (smiles)')
    self._conn.execute('PRAGMA synchronous = OFF')
    self._conn.execute('PRAGMA journal_mode = MEMORY')
    self._conn.commit()

  def bulk_insert(self, encoded_conformers, batch_size=10000, limit=None):
    """Inserts conformers into the database.

    Args:
      encoded_conformers: iterable for encoded dataset_pb2.Conformer
      batch_size: insert performance is greatly improved by putting multiple
        insert into one transaction. 10k was a reasonable default from some
        early exploration.
      limit: maximum number of records to insert

    Raises:
      ReadOnlyError: if mode is 'r'
      ValueError: If encoded_conformers is empty.
    """
    if self._read_only:
      raise ReadOnlyError()
    if not encoded_conformers:
      raise ValueError()

    insert_conformer = (f'INSERT INTO {_CONFORMER_TABLE_NAME} '
                        'VALUES (?, ?, ?)')
    insert_btid = f'INSERT INTO {_BTID_TABLE_NAME} VALUES (?, ?)'
    insert_smiles = (
        f'INSERT OR IGNORE INTO {_SMILES_TABLE_NAME} VALUES (?, ?) ')

    cur = self._conn.cursor()

    start_time = datetime.datetime.now()

    pending_conformer_args = []
    pending_btid_args = []
    pending_smiles_args = []

    def commit_pending():
      cur.executemany(insert_conformer, pending_conformer_args)
      cur.executemany(insert_btid, pending_btid_args)
      cur.executemany(insert_smiles, pending_smiles_args)
      pending_conformer_args.clear()
      pending_btid_args.clear()
      pending_smiles_args.clear()
      self._conn.commit()

    idx = None
    for idx, encoded_conformer in enumerate(encoded_conformers, 1):
      conformer = dataset_pb2.Conformer.FromString(encoded_conformer)
      expanded_stoich = (
        smu_utils_lib.expanded_stoichiometry_from_topology(
          conformer.bond_topologies[0]))
      pending_conformer_args.append((conformer.conformer_id, expanded_stoich,
                                     snappy.compress(encoded_conformer)))
      for bond_topology in conformer.bond_topologies:
        pending_btid_args.append(
            (bond_topology.bond_topology_id, conformer.conformer_id))
        pending_smiles_args.append(
            (bond_topology.smiles, bond_topology.bond_topology_id))
      if batch_size and idx % batch_size == 0:
        commit_pending()
        elapsed = datetime.datetime.now() - start_time
        logging.info(
            'bulk_insert: committed at index %d, %f s total, %.6f s/record',
            idx, elapsed.total_seconds(),
            elapsed.total_seconds() / idx)

      if limit and idx >= limit:
        break

    # Commit a final time
    commit_pending()
    elapsed = datetime.datetime.now() - start_time
    logging.info('bulk_insert: Total records %d, %f s, %.6f s/record', idx,
                 elapsed.total_seconds(),
                 elapsed.total_seconds() / idx)

  def bulk_insert_smiles(self, smiles_btid_pairs, batch_size=10000):
    """Insert smiles to bond topology id mapping.

    Args:
      smiles_btid_pairs: iterable of pairs of (smiles, btid)
    """
    if self._read_only:
      raise ReadOnlyError()

    insert_smiles = (
        f'INSERT OR IGNORE INTO {_SMILES_TABLE_NAME} VALUES (?, ?) ')

    cur = self._conn.cursor()

    pending = []

    def commit_pending():
      cur.executemany(insert_smiles, pending)
      pending.clear()
      self._conn.commit()

    for idx, (smiles, btid) in enumerate(smiles_btid_pairs, 1):
      pending.append([smiles, btid])
      if batch_size and idx % batch_size == 0:
        commit_pending()

    # Commit a final time
    commit_pending()

  def vacuum(self):
    """Uses SQL VACUUM to clean up db.

    Args:
      filename to write to
    """
    if self._read_only:
      raise ReadOnlyError()
    cur = self._conn.cursor()
    cur.execute('VACUUM')
    self._conn.commit()

  def find_bond_topology_id_for_smiles(self, smiles):
    """Finds the bond_topology_id for the given smiles.

    Args:
      smiles: string to look up

    Returns:
      integer of bond_topology_id

    Raises:
      KeyError: if smiles not found
    """
    cur = self._conn.cursor()
    select = f'SELECT btid FROM {_SMILES_TABLE_NAME} WHERE smiles = ?'
    cur.execute(select, (smiles,))
    result = cur.fetchall()

    if not result:
      raise KeyError(f'SMILES {smiles} not found')

    # Since it's a unique index, there should only be one result and it's a
    # tuple with one value.
    assert len(result) == 1
    assert len(result[0]) == 1
    return result[0][0]

  def find_by_conformer_id(self, cid):
    """Finds the conformer associated with a conformer id.

    Args:
      cid: conformer id to look up.

    Returns:
      dataset_pb2.Conformer

    Raises:
      KeyError: if cid is not found
    """
    cur = self._conn.cursor()
    select = f'SELECT conformer FROM {_CONFORMER_TABLE_NAME} WHERE cid = ?'
    cur.execute(select, (cid,))
    result = cur.fetchall()

    if not result:
      raise KeyError(f'Conformer id {cid} not found')

    # Since it's a unique index, there should only be one result and it's a
    # tuple with one value.
    assert len(result) == 1
    assert len(result[0]) == 1
    return dataset_pb2.Conformer().FromString(snappy.uncompress(result[0][0]))

  def find_by_bond_topology_id_list(self, btids):
    """Finds all the conformer associated with a bond topology id.

    Args:
      btids: list of bond topology id to look up.

    Returns:
      iterable of dataset_pb2.Conformer
    """
    cur = self._conn.cursor()
    # DISTINCT is because the same cid can have the same btid multiple times.
    select = (''.join([
      f'SELECT DISTINCT cid, conformer '
      f'FROM {_CONFORMER_TABLE_NAME} '
      f'INNER JOIN {_BTID_TABLE_NAME} USING(cid) '
      f'WHERE {_BTID_TABLE_NAME}.btid IN (',
      ','.join('?' for _ in btids),
      ')']))
    cur.execute(select, btids)
    return (dataset_pb2.Conformer().FromString(snappy.uncompress(result[1]))
            for result in cur)

  def find_by_smiles_list(self, smiles):
    """Finds all conformer associated with a given smiles string.

    Args:
      smiles: list of string

    Returns:
      iterable for dataset_pb2.Conformer
    """
    canon_smiles = [smu_utils_lib.compute_smiles_for_molecule(
        Chem.MolFromSmiles(s, sanitize=False), include_hs=False)
                    for s in smiles]
    cur = self._conn.cursor()
    select = (''.join([
      f'SELECT btid FROM {_SMILES_TABLE_NAME} WHERE smiles IN (',
      ','.join('?' for _ in canon_smiles),
      ')']))
    cur.execute(select, canon_smiles)
    result = cur.fetchall()

    if not result:
      return []

    return self.find_by_bond_topology_id_list([r[0] for r in result])

  def find_by_expanded_stoichiometry_list(self, exp_stoichs):
    """Finds all of the conformers with a stoichiometry.

    The expanded stoichiometry includes hydrogens as part of the atom type.
    See smu_utils_lib.expanded_stoichiometry_from_topology for a
    description.

    Args:
      exp_stoichs: list of string

    Returns:
      iterable of dataset_pb2.Conformer
    """
    cur = self._conn.cursor()
    select = (''.join([
      f'SELECT conformer '
      f'FROM {_CONFORMER_TABLE_NAME} '
      f'WHERE exp_stoich IN (',
      ','.join('?' for _ in exp_stoichs),
      ')']))
    cur.execute(select, exp_stoichs)
    return (dataset_pb2.Conformer().FromString(snappy.uncompress(result[0]))
            for result in cur)

  def find_by_stoichiometry(self, stoich):
    """Finds all conformers with a given stoichiometry.

    The stoichiometry is like "C6H12".

    Internally, the stoichiometry is converted a set of expanded stoichiometries
    and the query is done to find all of those.

    Args:
      stoich: stoichiometry string like "C6H12", case doesn't matter
    Returns:
      Iterable of type dataset_pb2.Conformer.
    """
    exp_stoichs = list(
        smu_utils_lib.expanded_stoichiometries_from_stoichiometry(stoich))
    return self.find_by_expanded_stoichiometry_list(exp_stoichs)

  def find_by_topology(self, smiles, bond_lengths,
                       matching_parameters=smu_molecule.MatchingParameters()):
    """Find all conformers which have a detected bond topology.

    Note that this *redoes* the detection. If you want the default detected
    versions, you can just query by SMILES string. This is only useful if you
    adjust the distance thresholds for what a matching bond is.
    To adjust those, you probably want to use
    AllAtomPairLengthDistributions.add_from_string_spec

    Args:
      smiles: smiles string for the target bond topology
      bond_lengths: AllAtomPairLengthDistributions
      matching_parameters: controls the algorithm for matching topologies.
        Generally should not need to be modified.

    Yields:
      dataset_pb2.Conformer
    """
    query_bt = smu_utils_lib.molecule_to_bond_topology(
      smu_utils_lib.smiles_to_molecule(smiles))
    expanded_stoich = smu_utils_lib.expanded_stoichiometry_from_topology(
      query_bt)
    cnt_matched_conformer = 0
    cnt_conformer = 0
    logging.info('Starting query for %s with stoich %s',
                 smiles, expanded_stoich)
    for conformer in self.find_by_expanded_stoichiometry_list([expanded_stoich]):
      if not smu_utils_lib.conformer_eligible_for_topology_detection(conformer):
        continue
      cnt_conformer += 1
      matches = topology_from_geom.bond_topologies_from_geom(
          conformer,
          bond_lengths=bond_lengths,
          matching_parameters=matching_parameters)
      if smiles in [bt.smiles for bt in matches.bond_topology]:
        cnt_matched_conformer += 1
        del conformer.bond_topologies[:]
        conformer.bond_topologies.extend(matches.bond_topology)
        for bt in conformer.bond_topologies:
          try:
            bt.bond_topology_id = self.find_bond_topology_id_for_smiles(
              bt.smiles)
          except KeyError:
            logging.error('Did not find bond topology id for smiles %s',
                          bt.smiles)
        yield conformer
    logging.info('Topology query for %s matched %d / %d', smiles,
                 cnt_matched_conformer, cnt_conformer)

  def find_bond_topology_id_by_smarts(self, smarts):
    """Find all bond topology ids that match a smarts pattern.

    Args:
      smarts: SMARTS string

    Yields:
      int, bond topology id
    """
    pattern = Chem.MolFromSmarts(smarts)
    if not pattern:
      raise ValueError(f'Could not parse SMARTS {smarts}')

    for smiles, bt_id in self.smiles_iter():
      mol = smu_utils_lib.smiles_to_molecule(smiles)
      # This is not the prettiest thing in the world. In order for ring markings
      # in the SMARTS to work, RingInfo has to be added. The simplest way to get
      # RingInfo set is to call this function. We didn't put this into the
      # smu_utils_lib just in case it messes something else up.
      Chem.GetSymmSSSR(mol)

      if mol.GetSubstructMatches(pattern):
        yield bt_id

  def smiles_iter(self):
    """Iterates through all (smiles, btid) pairs in the DB.

    Yields:
      (smiles, bt_id)
    """
    cur = self._conn.cursor()
    cur.execute('SELECT smiles, btid FROM smiles')
    yield from cur

  def __iter__(self):
    """Iterates through all dataset_pb2.Conformer in the DB."""
    select = f'SELECT conformer FROM {_CONFORMER_TABLE_NAME} ORDER BY rowid'
    cur = self._conn.cursor()
    cur.execute(select)
    return (dataset_pb2.Conformer().FromString(snappy.uncompress(result[0]))
            for result in cur)

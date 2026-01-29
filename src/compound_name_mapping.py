"""
Compound Name Mapping Module
============================

Module for mapping BRD numbers to correct compound names.
Retrieves compound names from PubChem and ChEMBL, and uploads to Snowflake table.

Author: Claude Code
Date: 2026-01-20
"""

import pandas as pd
import numpy as np
import requests
import json
import time
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, List
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CompoundNameMapper:
    """Compound name mapping class"""

    def __init__(self, cache_dir: Path = None, rate_limit: float = 0.1, timeout: int = 5):
        """
        Parameters
        ----------
        cache_dir : Path
            Cache storage directory
        rate_limit : float
            API request interval (seconds)
        timeout : int
            API timeout (seconds)
        """
        self.cache_dir = cache_dir or Path('results/compound_name_mapping')
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.rate_limit = rate_limit
        self.timeout = timeout
        self._lock = threading.Lock()

        # Load cache
        self.pubchem_cache = self._load_cache('pubchem_cache.json')
        self.chembl_cache = self._load_cache('chembl_cache.json')

    def _load_cache(self, filename: str) -> Dict:
        """Load cache file"""
        cache_path = self.cache_dir / filename
        if cache_path.exists():
            with open(cache_path, 'r') as f:
                return json.load(f)
        return {}

    def _save_cache(self, cache: Dict, filename: str):
        """Save cache file (atomic write)"""
        cache_path = self.cache_dir / filename
        temp_path = cache_path.with_suffix('.tmp')
        with self._lock:
            with open(temp_path, 'w') as f:
                json.dump(cache, f, indent=2, ensure_ascii=False)
            temp_path.replace(cache_path)

    def query_pubchem(self, inchi_key: str) -> Dict:
        """
        Retrieve compound information from PubChem using InChIKey

        Returns
        -------
        Dict with keys: pubchem_cid, pubchem_name, iupac_name, molecular_formula
        """
        if not inchi_key or pd.isna(inchi_key):
            return {}

        # Check cache
        if inchi_key in self.pubchem_cache:
            return self.pubchem_cache[inchi_key]

        result = {}
        try:
            # Get CID and basic information
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/property/IUPACName,MolecularFormula,Title/JSON"
            resp = requests.get(url, timeout=self.timeout)

            if resp.status_code == 200:
                data = resp.json()
                props = data['PropertyTable']['Properties'][0]
                result = {
                    'pubchem_cid': props.get('CID'),
                    'iupac_name': props.get('IUPACName', ''),
                    'molecular_formula': props.get('MolecularFormula', ''),
                }

                # Get Title (from synonyms)
                syn_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{props['CID']}/synonyms/JSON"
                syn_resp = requests.get(syn_url, timeout=self.timeout)
                if syn_resp.status_code == 200:
                    syn_data = syn_resp.json()
                    synonyms = syn_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])

                    # Exclusion patterns (ID-format names)
                    exclude_prefixes = (
                        'BRD', 'MLS', 'CHEMBL', 'SCHEMBL', 'CHEBI', 'HMS',
                        'ZINC', 'AKOS', 'DTXSID', 'NCGC', 'SID', 'CID',
                        'MFCD', 'EN300', 'BBL', 'STK', 'ACM', 'AC1', 'AB0',
                        'SMR', 'UNII', 'DSS', 'LS-', 'NSC', 'SR-', 'CS-',
                        'CCG-', 'EINECS', 'HSDB', 'EPA', 'CCRIS'
                    )

                    # Find a good name
                    for name in synonyms[:30]:
                        name_upper = name.upper()
                        # Skip exclusion patterns
                        if any(name_upper.startswith(p) for p in exclude_prefixes):
                            continue
                        # Skip CAS number format
                        if len(name) < 15 and '-' in name and name.replace('-', '').isdigit():
                            continue
                        # Skip numbers only
                        if name.replace('-', '').replace(' ', '').isdigit():
                            continue
                        # Skip too long names
                        if len(name) > 50:
                            continue
                        # Adopt as valid name
                        result['pubchem_name'] = name
                        break

                    # If not found, use shortened IUPAC name
                    if 'pubchem_name' not in result:
                        result['pubchem_name'] = ''

        except Exception as e:
            logger.debug(f"PubChem query failed for {inchi_key}: {e}")

        # Save cache (thread-safe)
        with self._lock:
            self.pubchem_cache[inchi_key] = result
        time.sleep(self.rate_limit)

        return result

    def query_chembl(self, inchi_key: str) -> Dict:
        """
        Retrieve compound information from ChEMBL using InChIKey

        Returns
        -------
        Dict with keys: chembl_id, chembl_name
        """
        if not inchi_key or pd.isna(inchi_key):
            return {}

        # Check cache
        if inchi_key in self.chembl_cache:
            return self.chembl_cache[inchi_key]

        result = {}
        try:
            # Search using first part of InChIKey
            inchi_key_base = inchi_key.split('-')[0] if '-' in inchi_key else inchi_key
            url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key__startswith={inchi_key_base}"
            resp = requests.get(url, timeout=self.timeout)

            if resp.status_code == 200:
                data = resp.json()
                molecules = data.get('molecules', [])
                if molecules:
                    mol = molecules[0]
                    result = {
                        'chembl_id': mol.get('molecule_chembl_id', ''),
                        'chembl_name': mol.get('pref_name', ''),
                    }

        except Exception as e:
            logger.debug(f"ChEMBL query failed for {inchi_key}: {e}")

        # Save cache (thread-safe)
        with self._lock:
            self.chembl_cache[inchi_key] = result
        time.sleep(self.rate_limit)

        return result

    def determine_compound_name(self, row: pd.Series) -> Tuple[str, str]:
        """
        Determine the final compound name

        Priority:
        1. ChEMBL pref_name (if drug name exists)
        2. PubChem name (common name)
        3. IUPAC name (if short)
        4. Original pertname (fallback)

        Returns
        -------
        Tuple[str, str]: (compound_name, name_source)
        """
        pertname = row.get('pertname', '')
        chembl_name = row.get('chembl_name', '')
        pubchem_name = row.get('pubchem_name', '')
        iupac_name = row.get('iupac_name', '')

        # 1. ChEMBL pref_name (drug name)
        if chembl_name and not chembl_name.startswith('CHEMBL'):
            return chembl_name, 'chembl'

        # 2. PubChem name (common name)
        if pubchem_name and len(pubchem_name) < 100:
            return pubchem_name, 'pubchem'

        # 3. IUPAC name (only if short)
        if iupac_name and len(iupac_name) < 60:
            return iupac_name, 'iupac'

        # 4. Original pertname (fallback)
        return pertname, 'original'

    def process_compounds(self, df: pd.DataFrame,
                          save_interval: int = 100) -> pd.DataFrame:
        """
        Process compound list and add name information

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing pertname, pertid, inchi_key, canonical_smiles
        save_interval : int
            Cache save interval

        Returns
        -------
        pd.DataFrame
            DataFrame with name mapping information added
        """
        total = len(df)
        results = []

        for i, (_, row) in enumerate(df.iterrows()):
            inchi_key = row.get('inchi_key', '')

            # API retrieval
            pubchem_info = self.query_pubchem(inchi_key)
            chembl_info = self.query_chembl(inchi_key)

            # Create row data
            result = {
                'pertname': row['pertname'],
                'pertid': row.get('pertid', ''),
                'inchi_key': inchi_key,
                'canonical_smiles': row.get('canonical_smiles', ''),
                'pubchem_cid': pubchem_info.get('pubchem_cid'),
                'pubchem_name': pubchem_info.get('pubchem_name', ''),
                'iupac_name': pubchem_info.get('iupac_name', ''),
                'molecular_formula': pubchem_info.get('molecular_formula', ''),
                'chembl_id': chembl_info.get('chembl_id', ''),
                'chembl_name': chembl_info.get('chembl_name', ''),
            }

            # Determine final compound name
            compound_name, name_source = self.determine_compound_name(pd.Series(result))
            result['compound_name'] = compound_name
            result['name_source'] = name_source

            results.append(result)

            # Progress display
            if (i + 1) % 50 == 0:
                logger.info(f"Processed {i + 1}/{total} compounds")

            # Save cache
            if (i + 1) % save_interval == 0:
                self._save_cache(self.pubchem_cache, 'pubchem_cache.json')
                self._save_cache(self.chembl_cache, 'chembl_cache.json')

        # Final cache save
        self._save_cache(self.pubchem_cache, 'pubchem_cache.json')
        self._save_cache(self.chembl_cache, 'chembl_cache.json')

        df_result = pd.DataFrame(results)
        df_result['created_at'] = datetime.now()
        df_result['updated_at'] = datetime.now()

        return df_result

    def _process_single_compound(self, row: pd.Series) -> Dict:
        """Process single compound (for parallel processing)"""
        inchi_key = row.get('inchi_key', '')

        # API retrieval
        pubchem_info = self.query_pubchem(inchi_key)
        chembl_info = self.query_chembl(inchi_key)

        # Create row data
        result = {
            'pertname': row['pertname'],
            'pertid': row.get('pertid', ''),
            'inchi_key': inchi_key,
            'canonical_smiles': row.get('canonical_smiles', ''),
            'pubchem_cid': pubchem_info.get('pubchem_cid'),
            'pubchem_name': pubchem_info.get('pubchem_name', ''),
            'iupac_name': pubchem_info.get('iupac_name', ''),
            'molecular_formula': pubchem_info.get('molecular_formula', ''),
            'chembl_id': chembl_info.get('chembl_id', ''),
            'chembl_name': chembl_info.get('chembl_name', ''),
        }

        # Determine final compound name
        compound_name, name_source = self.determine_compound_name(pd.Series(result))
        result['compound_name'] = compound_name
        result['name_source'] = name_source

        return result

    def process_compounds_parallel(self, df: pd.DataFrame,
                                   max_workers: int = 4,
                                   save_interval: int = 500) -> pd.DataFrame:
        """
        Process compound list in parallel and add name information

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing pertname, pertid, inchi_key, canonical_smiles
        max_workers : int
            Number of parallel workers
        save_interval : int
            Cache save interval

        Returns
        -------
        pd.DataFrame
            DataFrame with name mapping information added
        """
        total = len(df)
        results = []
        processed = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(self._process_single_compound, row): i
                       for i, (_, row) in enumerate(df.iterrows())}

            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.warning(f"Error processing compound: {e}")

                processed += 1
                if processed % 100 == 0:
                    logger.info(f"Processed {processed}/{total} compounds")

                    # Periodically save cache
                    with self._lock:
                        self._save_cache(self.pubchem_cache, 'pubchem_cache.json')
                        self._save_cache(self.chembl_cache, 'chembl_cache.json')

        # Final cache save
        self._save_cache(self.pubchem_cache, 'pubchem_cache.json')
        self._save_cache(self.chembl_cache, 'chembl_cache.json')

        df_result = pd.DataFrame(results)
        df_result['created_at'] = datetime.now()
        df_result['updated_at'] = datetime.now()

        return df_result

    def save_to_parquet(self, df: pd.DataFrame, path: Path = None):
        """Save to local parquet"""
        if path is None:
            path = self.cache_dir / 'compound_mapping.parquet'
        df.to_parquet(path, index=False)
        logger.info(f"Saved to {path}")


def upload_to_snowflake(df: pd.DataFrame, conn, table_name: str = 'COMPOUND_NAME_MAPPING'):
    """
    Upload DataFrame to Snowflake

    Parameters
    ----------
    df : pd.DataFrame
        Data to upload
    conn : snowflake.connector.connection
        Snowflake connection
    table_name : str
        Table name
    """
    from snowflake.connector.pandas_tools import write_pandas

    # Convert column names to uppercase
    df_upload = df.copy()
    df_upload.columns = [c.upper() for c in df_upload.columns]

    # Upload to Snowflake
    success, nchunks, nrows, _ = write_pandas(
        conn,
        df_upload,
        table_name,
        database='BIOINFORMATICS',
        schema='LINCS',
        auto_create_table=True,
        overwrite=True
    )

    if success:
        logger.info(f"Uploaded {nrows} rows to {table_name}")
    else:
        logger.error(f"Failed to upload to {table_name}")

    return success

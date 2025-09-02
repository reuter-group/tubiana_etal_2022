#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2022 Reuter Group
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


""" PePrMInt dataset creation (previously on Notebook #2)

__author__ = ["Thibault Tubiana", "Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2022 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

import pandas as pd
import numpy as np
import os
import importlib
import pickle

from tqdm import trange, tqdm
from typing import Optional

from src.settings import Settings
from src.notebook_handle import NotebookHandle

# modules for objects in Notebook #3 and the auxiliary tools one (for IBS)
from src.alphafold_utils import AlphaFoldUtils
from src.ibs_tagging import IBSTagging

# module for object in Notebook #4
from src.figure_generator import FigureGenerator

class DatasetManager:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings

        if self.settings.XP_MODE:
            self._FULL_DATASET_FILENAME = "DATASET_peprmint_allatoms_XP"
            self._LIGHT_DATASET_FILENAME = "DATASET_peprmint_XP"
            self._IBS_FROM_CATH_DATASET_FILENAME = "TAGGED_MERGED_DS_FROM_CATH_XP"
            self._IBS_FROM_AF_DATASET_FILENAME = "TAGGED_MERGED_DS_FROM_AF_XP"
            self._IBS_FROM_CATHAF_DATASET_FILENAME = "TAGGED_MERGED_DS_FROM_CATHAF_XP"
        else:
            self._FULL_DATASET_FILENAME = "DATASET_peprmint_allatoms_d25"
            self._LIGHT_DATASET_FILENAME = "DATASET_peprmint_d25"
            self._IBS_FROM_CATH_DATASET_FILENAME = "TAGGED_MERGED_DS_FROM_CATH"
            self._IBS_FROM_AF_DATASET_FILENAME = "TAGGED_MERGED_DS_FROM_AF"
            self._IBS_FROM_CATHAF_DATASET_FILENAME = "TAGGED_MERGED_DS_FROM_CATHAF"

        self.DATASET = None
        self.alphafold_utils = None
        self.IBS_tagger = None
        self.figure_generator = None

        self.RECALCULATION = self.settings.config_file.getboolean(
            'DATASET_MANAGER', 'recalculate')

        self._libs_setup()

    def _libs_setup(self):
        # Pandas
        pd.options.mode.chained_assignment = (
            None  # remove warning when adding a new column; default='warn'
        )
        pd.set_option("display.max_columns", None)
        tqdm.pandas()   # activate tqdm progressbar for pandas

        # IPython
        if self.settings.USING_NOTEBOOK:
            self.settings.NOTEBOOK_HANDLE.dataset_manager_options()

    def load_full_dataset(self) -> bool:
        # loads the dataset built in a previous run (this full version is NOT 
        # used in the original notebooks)
        full_path = self.settings.WORKDIR + self._FULL_DATASET_FILENAME + ".pkl"
        if not os.path.isfile(full_path):
            print("Could not find the full dataset file ({0})".format(full_path))
            print("Use build() from DatasetManager to create it")
            return False
        else:
            self.DATASET = pd.read_pickle(full_path)
            print("Dataset (full version) loaded successfully")
            return True

    def load_light_dataset(self) -> bool:
        # loads the dataset built in a previous run (this light version is the 
        # one used throughout the original notebooks)
        full_path = self.settings.WORKDIR + self._LIGHT_DATASET_FILENAME + ".pkl"
        if not os.path.isfile(full_path):
            print("Could not find the light dataset file ({0})".format(full_path))
            print("Use build() from DatasetManager to create it")
            return False
        else:
            self.DATASET = pd.read_pickle(full_path)
            print("Dataset (light version) loaded successfully")
            return True

    def build(self, recalculate: Optional[bool] = None):
        # runs each step in the creation of the dataset (as of Notebook #2)
        if recalculate is not None:
            self.RECALCULATION = recalculate

        self._pepr2ds_setup()

        self.clean()
        self.compute_protusion()
        self.add_cluster_structural_info()
        self.add_uniprot_basic_info()
        self.add_prosite_info()
        self.add_sequences_without_structure()
        self.add_uniprot_protein_sheet_info()
        self.add_cluster_uniref_info()
        self.add_conservation()
        self.save_dataset()

        print("\nDataset built successfully")
        print("Dataset domains: ")
        print(list(self.DATASET.domain.unique()))
        print("Dataset 'data_type' in: ")
        print(list(self.DATASET.data_type.unique()))

    def fetch_alphafold_data(self,
                             EXCLUDE_SEQS: Optional[list] = None,
                             EXCLUDE_DOMAIN: Optional[list] = None):

        # Fetches AF data for domains in the current dataset
        # Originally on Notebook #3; ported to alphafold_utils.py
        self.alphafold_utils = AlphaFoldUtils(self.settings)
        self.alphafold_utils.run(self.DATASET,
                                 EXCLUDE_SEQS,
                                 EXCLUDE_DOMAIN)

        # TO DO: add preprocessing option here or leave it to the caller method?
        print("AlphaFold data fetched successfully; update dataset manager with build(recalculate=True)")

    def add_IBS_data(self, db="cath+af"):
        # Interfacial binding sites (IBS) tagging in the dataset
        # Originally on the "tools notebooks"; ported to ibs_tagging.py
        path_to_merged_file = self.settings.WORKDIR
        if db == "cath":
            path_to_merged_file += self._IBS_FROM_CATH_DATASET_FILENAME + ".pkl"
        elif db == "alphafold":
            path_to_merged_file += self._IBS_FROM_AF_DATASET_FILENAME + ".pkl"
        elif db == "cath+af":
            path_to_merged_file += self._IBS_FROM_CATHAF_DATASET_FILENAME + ".pkl"
        else:
            print(f"> Warning: did not recognize argument '{db}' to add_IBS_data")
            print(f"           - defaulting to 'cath+af'")
            db = "cath+af"
            path_to_merged_file += self._IBS_FROM_CATHAF_DATASET_FILENAME + ".pkl"

        self.IBS_tagger = IBSTagging(self.settings, data_type=db)
        self.IBS_tagger.run(self.DATASET)
        self.IBS_tagger.make_analysis_report()

        # TO DO: do we need to build again!? Figure out if the merged Dataset object in IBSTagging should be reflected here
        #print("Updating dataset with IBS data")
        #self.build(recalculate=False)

        # serialize tagger object e.g. to generate figures later on
        with open(path_to_merged_file, "wb") as outfile:
            pickle.dump(self.IBS_tagger, outfile)

    def load_IBS_data(self, db="cath+af") -> bool:
        path_to_merged_file = self.settings.WORKDIR
        if db == "cath":
            path_to_merged_file += self._IBS_FROM_CATH_DATASET_FILENAME + ".pkl"
        elif db == "alphafold":
            path_to_merged_file += self._IBS_FROM_AF_DATASET_FILENAME + ".pkl"
        elif db == "cath+af":
            path_to_merged_file += self._IBS_FROM_CATHAF_DATASET_FILENAME + ".pkl"
        else:
            print(f"> Warning: did not recognize argument '{db}' to load_IBS_data")
            print(f"           - defaulting to 'cath+af'")
            db = "cath+af"
            path_to_merged_file += self._IBS_FROM_CATHAF_DATASET_FILENAME + ".pkl"
        
        if not os.path.isfile(path_to_merged_file):
            print("Could not find the tagged dataset file ({0})".format(path_to_merged_file))
            print("Use add_IBS_data() from DatasetManager to compute it")
            return False
        else:
            with open(path_to_merged_file, "rb") as infile:
                self.IBS_tagger = pickle.load(infile)
            print("Tagged dataset loaded successfully")
            return True

    def get_figure_generator_after_IBS(self):
        if self.IBS_tagger is None:
            print("Error: (IBS) tagged dataset not computed/loaded")
        else:
            return FigureGenerator(self.settings, self.IBS_tagger.pepr2ds_dataset)

    def get_protusion_count_after_IBS(self, ibs_only=False):
        if self.IBS_tagger is None:
            print("Error: (IBS) tagged dataset not computed/loaded")
        else:
            return self.IBS_tagger._test_num_protrusions(ibs_only=ibs_only)

    """
    ### All methods below just encapsulate the steps in Notebook #2
    """

    def _pepr2ds_setup(self):
        import pepr2ds.builder.Builder as builderEngine
        importlib.reload(builderEngine)
        # TO DO: it seems that 'update' in pepr2ds only makes sense when 'recalculate' is True, right?
        self.builder = builderEngine.Builder(self.settings.SETUP, 
                                             recalculate = self.RECALCULATION,
                                             update = self.RECALCULATION,
                                             notebook = self.settings.USING_NOTEBOOK,
                                             core = self.settings.num_threads)

    def clean(self):
        print(">>> clean()...")
        self.builder.structure.clean_all_pdbs()
        self.DATASET = self.builder.structure.build_structural_dataset()

    def compute_protusion(self):
        print(">>> compute_protusion()...")
        self.DATASET = self.builder.structure.add_protrusions(self.DATASET)

    def add_cluster_structural_info(self):
        # TO DO: can we avoid recreating the builder object!?
        print(">>> _pepr2ds_setup()...")
        self._pepr2ds_setup()
        print(">>> add_cluster_structural_info()...")
        self.DATASET = self.builder.structure.add_structural_cluster_info(self.DATASET)

    def add_uniprot_basic_info(self):
        print(">>> add_uniprot_basic_info()...")
        try:
            self.DATASET = self.builder.sequence.add_uniprotId_Origin(self.DATASET)
        except Exception as e:
            print(f"Error adding UniProt basic info: {e}")

    def add_prosite_info(self):
        print(">>> add_prosite_info()...")
        self.DATASET = self.builder.sequence.match_residue_number_with_alignment_position(self.DATASET)

    def add_sequences_without_structure(self):
        print(">>> add_sequences_without_structure()...")
        self.DATASET = self.builder.sequence.add_sequence_in_dataset(self.DATASET)

    def add_uniprot_protein_sheet_info(self):
        print(">>> download_uniprot_data()...")
        # original method (individual HTTP requests) takes excessive time
        # self.builder.sequence.download_uniprot_data(self.DATASET)
        self.builder.sequence.download_uniprot_data_REST(self.DATASET)
        
        print(">>> add_info_from_uniprot()...")
        self.DATASET = self.builder.sequence.add_info_from_uniprot(self.DATASET)

    def add_cluster_uniref_info(self):
        print(">>> add_cluster_uniref_info()...")
        self.DATASET = self.builder.sequence.add_cluster_info(self.DATASET)

    def add_conservation(self):
        # TO DO: can we avoid recreating the builder object!?
        print(">>> _pepr2ds_setup()...")
        self._pepr2ds_setup()
        print(">>> add_conservation()...")
        self.DATASET = self.builder.sequence.add_conservation(self.DATASET,
                                                              gapcutoff=0.8)

    def save_dataset(self):
        self.DATASET = self.builder.optimize_size(self.DATASET)
        self.DATASET = self.DATASET.drop_duplicates(subset=['atom_number',
                                                            'atom_name',
                                                            'residue_name',
                                                            'residue_number',
                                                            'cathpdb',
                                                            'chain_id'])

        # save two versions of the dataset: a complete one (with all PDB atoms)
        # and a light one (only with `CA` and `CB` atoms).  
        self.builder.save_checkpoint_dataset(self.DATASET,
                                             self._FULL_DATASET_FILENAME)
        self.builder.save_checkpoint_dataset(self.DATASET.query("atom_name in ['CA','CB']"),
                                             self._LIGHT_DATASET_FILENAME)

    # TO DO: does not seem meaningful; consider removing in the future
    def test_alignment_for_C2DIS_domain(self):
        # generate alignment file list for C2DIS domain (S95, because everything is just too slow, too much structure)
        c2dis = self._selectUniquePerCluster(self.DATASET.query("domain== 'PH'"), 
                                             'S95', 
                                             'uniref90', 
                                             withAlignment=False, 
                                             pdbreference='2da0A00')
        #pdblist = c2dis.cathpdb.dropna().unique()
        #print(c2dis)
        print(self.DATASET.query("atom_name == 'CA' and domain =='PH'").columns)
        #print(self.DATASET.query("atom_name == 'CA' and domain =='PH' and data_type == 'cathpdb'")[["ASA_res_freesasa_florian","RSA_freesasa_florian","ASA_total_freesasa","ASA_mainchain_freesasa","ASA_sidechain_freesasa","RSA_sidechain_freesasa","RSA_total_freesasa_tien","RSA_sidechain_freesasa_tien"]])

    # TO DO: does not seem necessary; move to an ad-hoc module
    def _selectUniquePerCluster(self,
                                df, 
                                cathCluster, 
                                Uniref, 
                                withAlignment=True, 
                                pdbreference=None,
                                removeStrand=False):
        # return a dataset with only 1 datum per choosed cluster
        if cathCluster not in ["S35", "S60", "S95", "S100"]:
            raise ValueError('CathCluster given not in ["S35","S60","S95","S100"]')

        if Uniref not in ["uniref50", "uniref90", "uniref100"]:
            raise ValueError('CathCluster given not in ["uniref50","uniref90","uniref100"]')

        if withAlignment:
            df = df[~df.alignment_position.isnull()]

        cathdf = df.query("data_type == 'cathpdb'")
        seqdf = df.query("data_type == 'prosite'")

        def selectUniqueCath(group):
            uniqueNames = group.cathpdb.unique()
            if pdbreference:
                if pdbreference in uniqueNames:
                    select = pdbreference
                else:
                    select = uniqueNames[0]
            else:
                select = uniqueNames[0]

            # return group.query("cathpdb == @select")
            return select

        def selectUniqueUniref(group, exclusion):
            uniqueNames = group.uniprot_acc.unique()
            select = uniqueNames[0]
            # return group.query("uniprot_acc == @select")
            if select not in exclusion:
                return select

        # structures prior to sequences
        dfReprCathNames = cathdf.groupby(["domain", cathCluster]).apply(selectUniqueCath).to_numpy()
        print(dfReprCathNames)
        excludeUniref = df.query("cathpdb in @dfReprCathNames").uniprot_acc.unique()

        dfReprUnirefNames = seqdf.groupby(["domain", Uniref]).apply(selectUniqueUniref,exclusion=excludeUniref).to_numpy()
        dfReprCath = cathdf.query("cathpdb in @dfReprCathNames")
        dfReprUniref = seqdf.query("uniprot_acc in @dfReprUnirefNames")

        return (pd.concat([dfReprCath, dfReprUniref]))

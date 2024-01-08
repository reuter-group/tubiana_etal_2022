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


""" Data download (previously on Notebook #1)

The only difference here is that we moved notebook-related imports and initial 
settings to NotebookHandle. The Settings constructor checks whether we are
using a jupyter notebook and reflects the necessary configuration.

In this module, fetch() is the main (proxy) method for client classes, 
downloading most of the necessary data.

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
from collections import defaultdict

import os
import shutil
from pathlib import Path
import glob
import subprocess
import urllib
import tarfile
import gzip

from src.settings import Settings

class DataRetriever:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.UPDATE = False

    def fetch(self) -> bool:
        """ Fetches the main data needed to construct the dataset, from CATH,  
        EBI, and UniProt. Not included in this module are AlphaFold PDBs (that 
        were ported to a separate module, following Notebook #3) and the UniProt 
        xml sequence files (that are retrieved as part of the earlier  
        architecture of pepr2ds (called through DatasetManager).
        """

        if self.settings.FORMER_WORKING_DIR:
            print('WARNING! Working directory already exists! It is recommended to remove it if')
            print('         you want to fetch all data from scratch with DataRetriever.fetch()')

        print(f'> Data retriever (domains selected in .config file: {self.settings.active_superfamilies})')
        
        self._retrieve_cath_domains()
        self._retrieve_uniprot_to_pdb_correspondence()
        self._retrieve_prosite()
        self._retrieve_cath_pdb_files()
        self._copy_ref_pdb_files()
        
        return True

    def populate_from_local_data(self, path_to_data):
        """ Skip the CATH download step and copy files in the given directory
        instead. NB! Assuming that the files corresponding to superfamilies
        active in the peprmint_default.config file are available in the local
        folder!
        """
        if self.settings.FORMER_WORKING_DIR:
            print('WARNING! Working directory already exists! Please remove it before proceeding in this use case')

        print(f'> Data retriever **WITH LOCAL FILES** (domains selected in .config file: {self.settings.active_superfamilies})')
        
        self._copy_files_to_cath_folder(path_to_data)
        self._retrieve_cath_domains()
        self._retrieve_uniprot_to_pdb_correspondence()
        self._retrieve_prosite()
        
        return True

    def _copy_files_to_cath_folder(self, path_to_data):
        prefix = self.settings.CATHFOLDER
        aligned_suffix = self.settings.ALIGNED_SUBDIR

        for domain in self.settings.active_superfamilies:
            raw_folder = prefix + 'domains/' + domain + '/raw/'
            cleaned_folder = prefix + 'domains/' + domain + '/cleaned/'
            aligned_folder = prefix + 'domains/' + domain + '/' + aligned_suffix + '/'

            if not os.path.exists(raw_folder):
                os.makedirs(raw_folder)
            if not os.path.exists(cleaned_folder):
                os.makedirs(cleaned_folder)
            if not os.path.exists(aligned_folder):
                os.makedirs(aligned_folder)

            # copy each file into the expected folder
            filelist = Path(path_to_data).glob(f"{domain}/*.pdb")
            for f in filelist:
                out_name = aligned_folder + os.path.basename(f)
                exec_cmd = "cp " + str(f) + " " + out_name
                about = subprocess.run(exec_cmd, shell=True, capture_output=False)

            # save a copy of the reference structure in the expected folder
            ref_pdb = self.settings.config_file['PREPROCESSING']["ref_"+domain+"_pdb"]
            if ref_pdb is not None:
                original = aligned_folder + ref_pdb + '.pdb'
                out_name = self.settings.REF_FOLDER + ref_pdb + '.pdb'
                exec_cmd = "cp " + original + " " + out_name
                about = subprocess.run(exec_cmd, shell=True, capture_output=False)

    def _retrieve_cath_domains(self):
        # TO DO: release 4_2_0 or latest?
        self.dom_file = 'cath-domain-list.txt'
        url = self.settings.config_file['CATH']['domain_list_url']

        destination = self.settings.CATHFOLDER + self.dom_file
        if not os.path.isfile(destination) or self.UPDATE: 
            urllib.request.urlretrieve(url, destination)

        self.column_dom_file = [ 'Domain',
                                 'Class',
                                 'Architecture',
                                 'Topology',
                                 'Homologous',
                                 'S35',
                                 'S60',
                                 'S95',
                                 'S100',
                                 'S100Count',
                                 'DomSize',
                                 'resolution', ]

    def _retrieve_uniprot_to_pdb_correspondence(self):
        url = self.settings.config_file['UNIPROT']['url']

        destination = self.settings.CATHFOLDER + "pdb_chain_uniprot.csv"

        # NB! made a small fix here (open/write operations were outside the if)
        if not os.path.isfile(destination) or self.UPDATE:
            tmp_destination = destination + ".gz"
            
            urllib.request.urlretrieve(url, tmp_destination)
            
            with gzip.open(tmp_destination, 'rb') as f:
                file = f.read()
                with open(destination, 'wb') as output:
                    output.write(file)
            os.remove(tmp_destination)

    def _retrieve_prosite(self):
        url = self.settings.config_file['PROSITE']['url']

        destination = self.settings.PROSITEFOLDER + "prosite_alignments.tar.gz" 

        if self.UPDATE:
            if os.path.exists(self.settings.PROSITEFOLDER + "msa"):
                shutil.rmtree(self.settings.PROSITEFOLDER + "msa/")
            else:
                os.makedirs(self.settings.PROSITEFOLDER + "msa/")

        if not os.path.exists(self.settings.PROSITEFOLDER+"msa") or self.UPDATE:
            urllib.request.urlretrieve(url, destination)
            tf = tarfile.open(destination)
            tf.extractall(self.settings.PROSITEFOLDER)
            tf.close()
            os.rename(self.settings.PROSITEFOLDER + "prosite_alignments",
                      self.settings.PROSITEFOLDER + "msa")
            os.remove(destination)

    def _retrieve_cath_pdb_files(self):
        # cath_domains: the complete list from cath with all domains
        self.cath_domains = pd.read_csv(self.settings.CATHFOLDER + self.dom_file,
                                        comment = '#',
                                        sep = r"\s+",
                                        header = None)
        self.cath_domains.columns = self.column_dom_file

        if self.settings.PARALLEL:
            self.cath_domains['Superfamily'] = self.cath_domains.parallel_apply(
                lambda x : f"{x.Class}.{x.Architecture}.{x.Topology}.{x.Homologous}",
                axis=1)
        else:
            self.cath_domains['Superfamily'] = self.cath_domains.progress_apply(
                lambda x : f"{x.Class}.{x.Architecture}.{x.Topology}.{x.Homologous}",
                axis=1)

        # cath_superfamily: just superfamily-domain correspondence
        self.cath_superfamily = pd.DataFrame()
        self.cath_superfamily['Superfamily'] = self.cath_domains.Superfamily
        self.cath_superfamily['Domain'] = self.cath_domains.Domain

        # cath_domains_per_superfamily: a map superfamily -> list of domains (from cath_superfamily) 
        # NB! do not parallelize this one
        self.cath_domains_per_superfamily = defaultdict(list)
        _ = self.cath_superfamily.progress_apply(
                lambda x : self.cath_domains_per_superfamily[x.Superfamily].append(x.Domain),
                axis = 1)

        # gather data for domains enabled in config file
        for superfamily, domain in self.settings.SUPERFAMILY.items():
            if domain in self.settings.active_superfamilies:
                print(f"> Fetching domain {domain} files")
                self._fetch_dom_for_superfamily(superfamily, domain)

    def _fetch_dom_for_superfamily(self, superfamily, domName):
        prefix = self.settings.CATHFOLDER
        folder = prefix + 'domains/' + domName + '/raw/'
        if not os.path.exists(folder):
            os.makedirs(folder)
        if not os.path.exists(prefix + 'domains/' + domName + '/cleaned/'):
            os.makedirs(prefix + 'domains/' + domName + '/cleaned/')

        dom_list = self.cath_domains_per_superfamily[superfamily]

        # limit the dataset size when experimenting only
        if self.settings.XP_MODE and len(dom_list) > self.settings.xp_cath_limit:
                dom_list = dom_list[0:self.settings.xp_cath_limit]

                # make sure the reference pdb (to "align on z" later) is kept
                ref_pdb = self.settings.config_file['PREPROCESSING']["ref_"+domName+"_pdb"]
                if ref_pdb is not None and ref_pdb not in dom_list:
                    dom_list[0] = ref_pdb

        if self.settings.PARALLEL:
            pd.Series(dom_list).parallel_apply(
                lambda x : self._fetch_pdb_from_cath_dom(x, folder) )
        else:
            print(dom_list)
            pd.Series(dom_list).progress_apply(
                lambda x : self._fetch_pdb_from_cath_dom(x, folder) )

    def _fetch_pdb_from_cath_dom(self, dom, folder):
        url = self.settings.config_file['CATH']['fetch_pdb_url'] +  dom + ".pdb"
        destination = folder + dom + '.pdb'
        if not os.path.isfile(destination): 
            urllib.request.urlretrieve(url, destination)

    def _copy_ref_pdb_files(self):
        # keep a copy of the raw pdb for the reference protein in each superfamily
        for domain in self.settings.active_superfamilies:
            ref_pdb = self.settings.config_file['PREPROCESSING']["ref_"+domain+"_pdb"]
            if ref_pdb is not None:
                self._fetch_pdb_from_cath_dom(ref_pdb, self.settings.REF_FOLDER)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2023 Reuter Group
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


""" Rewriting the code to create and analyse data from the PePrMInt project.

As of now, we are ignoring the architecture of the future solution, and simply
porting the code previously developed in python notebooks by Thibault Tubiana 
towards a stand-alone script that we can test and maintain more easily.
In the next steps, we shall determine adequate modules and include options to
skip different phases of the workflow to allow researchers work with 
intermediate results.

__author__ = ["Thibault Tubiana", "Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2023 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

import os
import sys

from src.settings import Settings
from src.data_retriever import DataRetriever
from src.preprocessing import Preprocessing
from src.dataset_manager import DatasetManager
from src.figure_generator import FigureGenerator

def main():
    # notebook #0
    global_settings = Settings()   # setup reading standard configuration file
    #global_settings = Settings("/opt/cbu/my.config")  # different config file

    # support to specify superfamilies not implemented by default
    '''
    global_settings.add_new_superfamily(name = "C1",
                                        ref_pdb = "1ptrA00", # reference for tagging IBS
                                        ref_res1 = "243",
                                        ref_res2 = "257",
                                        ref_res3 = "237",
                                        cath_domain = "3.30.60.20",
                                        prosite_domain = "PS50081",
                                        interpro_domain = None,
                                        refine_AF_data_with_interpro = False)
    '''
 
    # notebook #1
    data_retriever = DataRetriever(global_settings)
    #data_retriever.fetch()
    data_retriever.populate_from_local_data(path_to_data='/scratch/cbu/tmp/')

    # superposition and reorientation of downloaded PDBs
    preprocess = Preprocessing(global_settings)
    #preprocess.run(database="cath", verbose=False, use_cath_superpose=False)

    # alternative, pairwise superimposition. NB! Might take too much time!
    #preprocess.run(database="cath", use_cath_superpose=True)

    # notebook #2
    dataset_manager = DatasetManager(global_settings)
    dataset_manager.build()   # build dataset from fetched data
    #dataset_manager.load_light_dataset()  # load dataset built on a previous run

    # notebook #3
    dataset_manager.fetch_alphafold_data()
    """
    dataset_manager.fetch_alphafold_data(EXCLUDE_SEQS = ["Q54C71",    # optional
                                                         "O94827",
                                                         "Q54C71",
                                                         "Q22070",
                                                         "P39960",
                                                         "Q62077",
                                                         "Q06839"],
                                         EXCLUDE_DOMAIN = ["PLA"])   # optional
    """
    preprocess.run(database="alphafold", verbose=True, use_cath_superpose=False)
    dataset_manager.build(recalculate=True)

    # auxiliary tools notebook for interfacial binding sites (IBS) tagging
    dataset_manager.add_IBS_data(db="cath+af")   # tag and save merged dataset
    #dataset_manager.load_IBS_data(db="cath+af")   # load from a previous run

    # notebook #4
    figure_gen = dataset_manager.get_figure_generator_after_IBS()
    figure_gen.make_figure_composition_of_exposed_IBS()
    figure_gen.make_figure_protrusions()
    figure_gen.make_figure_composition_for_proteins_with_HP_at_IBS()
    figure_gen.make_figure_neighbourhood_composition()
    figure_gen.make_figure_number_of_structures_w_and_wo_HP_at_IBS()
    figure_gen.make_figure_composition_for_proteins_without_HP_at_IBS()
    figure_gen.make_figure_superfamily_decomposition_exposed_env_HP()

if __name__ == '__main__':
    print('Running under Python {0[0]}.{0[1]}.{0[2]}'.format(sys.version_info),
        file=sys.stderr)
    main()

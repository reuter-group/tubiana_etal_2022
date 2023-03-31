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


""" Data transformation (previously in ad-hoc scripts e.g. align_on_z.py)

This module handles the superposition and the reorientation (along the z axis)
of downloaded PDBs. This ensures that the data is ready for the later tagging
of IBS (interfacial binding sites).

__author__ = ["Thibault Tubiana", "Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2022 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

import numpy as np
import vg
from tqdm import trange, tqdm

from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import os
import stat
from pathlib import Path
import glob
import subprocess
import warnings

import shutil
import urllib

from src.settings import Settings

class Preprocessing:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.superpose_tool = None

    def run(self) -> bool:
        """
        Proxy method for running the superposition over downloaded PDBs 
        (cath-superpose external tool, with built-in SSAP for all-vs-all 
        pairwise structure alignment) and the reorientation along the z axis.
        """
        
        print(f"> Superposing downloaded PDBs with cath-superpose and SSAP alignments")

        if not self._fetch_cath_superpose_binary():
            return False

        # TO DO: consider adding an option to avoid overwriting downloaded PDBs with superimposed ones
        if not self._run_cath_superpose():
            return False

        print(f"> Orienting superposed PDBs along z axis with respect to the domain representative")
        return self._orient_on_z_axis()


    ############################################################################
    # superpose-related methods below

    def _fetch_cath_superpose_binary(self) -> bool:
        if self.settings.OS == "linux":
            tool_url = self.settings.config_file['CATH']['superpose_url_linux']
        elif self.settings.OS == "macos":
            tool_url = self.settings.config_file['CATH']['superpose_url_macos']
        else:
            print('Error: current superposition method (cath-superpose) only available on GNU/Linux and MacOS')
            return False

        self.superpose_tool = self.settings.CATHFOLDER + 'cath-superpose'
        if not os.path.isfile(self.superpose_tool): 
            urllib.request.urlretrieve(tool_url, self.superpose_tool)

        # make executable: 'current permission bits' BINARY OR 'executable bit'
        st = os.stat(self.superpose_tool)
        os.chmod(self.superpose_tool, st.st_mode | stat.S_IEXEC)

        return True

    def _run_cath_superpose(self, verbose=False) -> bool:
        success_only = True
        error_message = ""

        for domain in self.settings.active_superfamilies:
            # temporary output folders
            # TO DO: if the alignments are not needed, remove the --do-the-ssaps option
            output_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/super/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            tmp_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/cath_tools_tmp/'
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)

            # domain files as arguments to cath-superpose
            filelist = Path(self.settings.CATHFOLDER).glob(f"domains/{domain}/raw/*.pdb")
            run_input = ""
            for f in filelist:
                run_input += "--pdb-infile " + str(f) + " "

            # need to export environment var for cath-ssap
            raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
            env_setting = f"export CATH_TOOLS_PDB_PATH={raw_dir} ; "

            # run cath-superpose
            run_options = "--do-the-ssaps " + tmp_dir
            run_output = "--sup-to-pdb-files-dir " + output_dir
            trigger = f"{self.superpose_tool} {run_options} {run_output} {run_input}"
            about = subprocess.run(env_setting+trigger, shell=True, capture_output=True)
            if verbose:
                print("done executing:")
                x = str(about.args)
                print(x)
                print("return code is:")
                x = int(about.returncode)
                print(str(x))
                print("stdout is:")
                x = str(about.stdout, "utf-8")
                print(x)
                print("stderr is")
                x = str(about.stderr, "utf-8")
                print(x)

            # TO DO: maybe call about.check_returncode() instead (raise exception)
            if int(about.returncode) != 0:
                success_only = False
                out = str(about.stderr, "utf-8")
                error_message += out

            # delete extra folder with alignments, replace raw with superposed
            shutil.rmtree(tmp_dir)
            raw_files = len([name for name in os.listdir(raw_dir) if os.path.isfile(os.path.join(raw_dir, name))])
            super_files = len([name for name in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, name))])
            if raw_files == super_files:
                shutil.rmtree(raw_dir)
                os.rename(output_dir, raw_dir)
            else:
                success_only = False
                out = f"Error: superposed file count for {domain} does not match raw file count\n"
                out += "Check directory: " + output_dir + "\n"
                error_message += out
        
        if not success_only:
            print("Errors found while executing cath-superpose:")
            print(error_message)

        return success_only
    
    ############################################################################
    # reorient-related methods below

    def _orient_on_z_axis(self, verbose=False) -> bool:
        success_only = True
        error_message = ""

        for domain in self.settings.active_superfamilies:
            ref_pdb = self.settings.config_file['REORIENT_ALONG_Z']["ref_"+domain+"_pdb"]
            if ref_pdb is None:
                print(f"  Warning: no reference protein for '{domain}' defined on .config file - skipping reorientation along z")
                continue
            
            # read pdbs in raw, save in temporary output folder, then replace them
            raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
            output_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/z_oriented/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            rotation, translation = self._get_transformation_from_reference(domain)
            self._transform_pdbs(domain, rotation, translation)

            if verbose:
                print(f"  '{domain}' orientation:")
                print(f"  rotation matrix: {rotation}")
                print(f"  translation vector: {translation}")

            # replace raw with reoriented
            raw_files = len([name for name in os.listdir(raw_dir) if os.path.isfile(os.path.join(raw_dir, name))])
            oriented_files = len([name for name in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, name))])
            if raw_files == oriented_files:
                shutil.rmtree(raw_dir)
                os.rename(output_dir, raw_dir)
            else:
                success_only = False
                out = f"Error: reoriented file count for '{domain}' does not match raw file count\n"
                out += "Check directory: " + output_dir + "\n"
                error_message += out
        
        if not success_only:
            print("Errors found while orienting along z axis:")
            print(error_message)

        return success_only

    def _get_transformation_from_reference(self, domain):
        raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
        ref_pdb = self.settings.config_file['REORIENT_ALONG_Z']["ref_"+domain+"_pdb"]
        structure =  self._get_structure(f"{raw_dir}/{ref_pdb}.pdb")

        res1 = self.settings.config_file.getint('REORIENT_ALONG_Z', 'ref_'+domain+'_res1')
        res2 = self.settings.config_file.getint('REORIENT_ALONG_Z', 'ref_'+domain+'_res2')
        res3 = self.settings.config_file.getint('REORIENT_ALONG_Z', 'ref_'+domain+'_res3')

        translation = self._get_translation_vector(structure, res1, res2, res3)
        
        structure = self._apply_translation_vector(structure, translation)
        rotation = self._get_rotation_matrix(structure, res1, res2, res3, domain)

        return (rotation, translation)

    def _transform_pdbs(self,
                        domain,
                        rotation: np.array,       # rotation matrix (3x3)
                        translation: np.array):   # translation vector (1x3)

        raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
        output_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/z_oriented/'
        pdb_list = glob.glob(f"{raw_dir}/*.pdb")

        for pdb_path in tqdm(pdb_list):
            structure = self._get_structure(pdb_path)
            structure = self._apply_translation_vector(structure,translation)
            structure = self._apply_rotation_matrix(structure, rotation)

            # write new file
            pdb_io = PDBIO()
            pdb_io.set_structure(structure)
            pdb_io.save(f"{output_dir}/{structure.id}.pdb")

    def _get_structure(self, pdb_path):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', PDBConstructionWarning)

                parser = PDB.PDBParser()
                pdb_code = Path(pdb_path).stem   # filename without extension
                return parser.get_structure(id = pdb_code, file = pdb_path)
        except:
            print(f"Reading Error: {pdb_path}")

    def _get_translation_vector(self,
                                structure: PDB.Structure.Structure,   # Biopython PDB Structure
                                res1: int,
                                res2: int,
                                res3: int) -> np.array:
        """
        Get the translation vector between the origin and the centroid of the
        triangle formed by the CA of three aminoacid residues.
        """
        chain = structure[0].child_list[0].id
        p1 = structure[0][chain][res1]['CA'].get_coord()
        p2 = structure[0][chain][res2]['CA'].get_coord()
        p3 = structure[0][chain][res3]['CA'].get_coord()

        # translation vector (1x3)
        return list(-np.mean(np.array([p1, p2, p3]), axis=0))

    def _apply_translation_vector(self,
                                  structure: PDB.Structure.Structure,  # Biopython PDB Structure
                                  translation: np.array):              # (1x3) translation vector
        """
        Apply translation without rotation.
        """
        null_rotation = np.identity(3).tolist()
        structure[0].transform(null_rotation, translation)
        return structure

    def _apply_rotation_matrix(self,
                               structure: PDB.Structure.Structure,  # Biopython PDB Structure
                               rotation:np.array):                  # (3x3) rotation matrix
        """
        Apply rotation without translation
        """
        structure[0].transform(rotation, [0, 0, 0])
        return structure

    def _get_rotation_matrix(self,
                             structure: PDB.Structure.Structure,   # Biopython PDB Structure
                             res1: int,
                             res2: int,
                             res3: int,
                             domain: str,
                             orientation = 'z') -> np.array:       # axis used for alignment
        """
        Get the rotation matrix between the normal of a triangle formed by the 
        CA of three aminoacid residues and a plane (x,y,z)
        """

        N, COM = self._get_normal_COM(structure, res1, res2, res3)

        # This norm will be our translation vector to all our atoms
        axis = {'x':[-1,0,0],
                'y':[0,-1,0],
                'z':[0,0,-1]}
        #Create the reference vector, per default we want to align on the z axis so it will be [0,0,1]
        refVector = PDB.vectors.Vector(axis[orientation])
        normal = PDB.vectors.Vector(N) #The normal should be a biopython object

        # Transformation 1
        temp = structure.copy()

        # Case 1 : normal and rotation
        rotation = PDB.vectors.rotmat(normal, refVector)
        temp.transform(rotation, [0, 0, 0])

        # If it doesn't work, case 2: -normal and rotation
        if not self._test_rotation(temp, res1, res2, res3):
            temp = structure.copy()
            rotation = PDB.vectors.rotmat(-normal, refVector)
            temp.transform(rotation, [0, 0, 0])
            # If it doesn't work, case3: normal and rotation.T
            if not self._test_rotation(temp, res1, res2, res3):
                temp = structure.copy()
                rotation = PDB.vectors.rotmat(normal, refVector).T
                temp.transform(rotation, [0, 0, 0])
                # If it doesn't work, case2: -normal and rotation.T
                if not self._test_rotation(temp, res1, res2, res3):
                    temp = structure.copy()
                    rotation = PDB.vectors.rotmat(-normal, refVector).T

        # (3x3) rotation matrix
        return rotation

    def _get_normal_COM(self, structure, res1, res2, res3):
        """
        Auxiliary method, to get the normal vector and and the geom. center of 
        a structure
        """

        # compute vectors (get new coordinates - rotation with null translation)
        chain = structure[0].child_list[0].id
        p1 = structure[0][chain][res1]['CA'].get_coord()
        p2 = structure[0][chain][res2]['CA'].get_coord()
        p3 = structure[0][chain][res3]['CA'].get_coord()

        # compute normal to the triangle: the cross product of two vectors
        A = p2 - p1
        B = p3 - p1
        N = np.cross(A, B)

        coords = np.array([x.get_coord() for x in structure.get_atoms()])
        COM = coords.mean(axis=0)

        return N, COM

    def _test_rotation(self, structure, res1, res2, res3):
        """
        Auxiliary method to recalculate angle, etc.
        """

        N, COM = self._get_normal_COM(structure, res1, res2, res3)

        # Recalculate angle
        angle = vg.angle(N, np.array([0, 0, -1]))
        return (angle == 0 and COM[2] > 0) or (angle == 180 and COM[2] > 0)
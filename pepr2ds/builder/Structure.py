from .Attributes import Attributes

import pandas as pd
import os,sys
import numpy as np
import requests
import pbxplore as pbx
import freesasa
from Bio import Entrez
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

from biopandas.pdb import PandasPdb
from scipy.spatial import ConvexHull
import scipy.spatial.distance as scidist

import shutil
from pathlib import Path


Entrez.email = "thibault.tubiana@uib.no"
Entrez.tool = "fetch sequence from UniprotID"

import warnings
warnings.filterwarnings("ignore")



class Structure(Attributes):
    def __init__(self, SETUP: dict, recalculate=True, update=True, notebook=True, core=4):
        super().__init__(SETUP, recalculate, update, notebook, core)


        self.ERRORS = []
        if notebook:
            from tqdm.notebook import tnrange, tqdm
            self.tqdm = tqdm
            self.tnrange = tnrange
            self.tqdm.pandas()  # activate tqdm progressbar for pandas apply
        else:
            print("notebook = False")
            from tqdm import tnrange, tqdm
            self.tqdm.pandas()  # activate tqdm progressbar for pandas apply


    def simplify_SS(self, ss):
        """
        Secondary Structure simplification from DSSP.
        All helix types are just helix
        All sheet type are just sheet
        turn, bend and coil are 'C' (coil)
        Args:
            ss (str): secondary structure letter
        Returns:
            str: simplified secondary structure letter.
        """
        if ss in ["H", "G", "I"]:
            return "H"
        elif ss in ["E", "B"]:
            return "E"
        elif ss in ["T", "S", "-"]:
            return "C"

    def smooth_ss(self, ss, window=3):
        """
        Smooth secondary structure by removing regions smaller than 3 amino acids
        and replacing it by the previous SS.
        for a helix, the minimum size is 3.6 to have a complet turn, so it will remove small helices
        for sheet, it will remove small flat sheet
        for coil, the minimum size appears to be 3 also (https://doi.org/10.1002/pro.5560051223)
        Example:
        original: CCCCCCSSHHHHHCCCCCCCCHHCCCCCSSSSSSHHH
        smoothed: CCCCCCCCHHHHHCCCCCCCCCCCCCCCSSSSSSHHH
        Args:
            ss <list>: secondary structure of the protein as a list (1 aa per element)
            windows <int>: smoothing window
        return:
            smoothedSS <list>: "simplified"/clean secondary structure list
        """

        smoothedSS = np.asarray(['C'] * len(ss))
        ssSegment = []

        # Definition, otherwise these variables will be local only
        newSegment = False
        stopSegment = None
        currentSS = None
        start = None
        segmentLength = 0

        i = 0
        while i < len(ss) - window:
            # If we have a non coil SS.
            if ss[i] in ['H', 'E']:
                # and if it's the begining of the segment
                if segmentLength == 0:
                    # setup the start,segmentLength and currentSS.
                    start = i
                    segmentLength = 1
                    currentSS = ss[i]
                    continueSegment = True
                # If we are in a in segment.
                while continueSegment and i < len(ss) - 2:
                    # get the next element
                    i += 1
                    if ss[i] == currentSS:  # if the next element is still the same SS, continue!
                        segmentLength += 1
                    else:  # Otherwise...
                        # Check the next+1 element (to detect patern like 'EECE')
                        if ss[i + 1] == currentSS:  # i t's a patern like 'EECE'
                            segmentLength += 1  # so we continue
                        else:  # end the segment
                            if segmentLength >= window:  # If the segment size is >= the minimum segment size
                                ssSegment.append(
                                    [currentSS, start, segmentLength])  # We add the segment in our segment list
                            # reset variables
                            stopSegment = None
                            currentSS = None
                            start = None
                            segmentLength = 0
                            continueSegment = False
                            i -= 1  # "cancel the last iteraction"

                # For the last SS elements.
                if len(ss) - 2 <= i < len(ss) - 1:
                    if segmentLength >= window:
                        ssSegment.append([currentSS, start, segmentLength])
            # General Coil iterator
            i += 1

        # Let's replace the secondary structure element.
        for segment in ssSegment:
            sselem = segment[0]
            start = segment[1]
            length = segment[2]
            smoothedSS[start:start + length] = sselem

        return smoothedSS

    def define_SS_segment(self, ssAtom):
        """
        Define continuous secondary structure segment segment
        Idea:
        CCCCHHHHCCEEEEHHHHHEEECCCCC
         C1  H1 C2 E1  H2   E2  C3

        NOTE: this can be modified to improve the coding of the SS segment
        Args:
            ssAtom (str): string of SS sequence (for every amino acids)
        Returns:
            SSegments (list): list of secondary structure segment.
        """
        SSelements = list(set(ssAtom))
        SScount = {x: 0 for x in SSelements}  # Counter for every amino acid type

        # instanciate first element
        SSsegments = []
        for i in range(len(ssAtom)):
            # check if the segment change
            currentSS = ssAtom[i]
            if i == 0:  # increment the first element since we start with 0
                SScount[currentSS] += 1
            else:
                if currentSS != ssAtom[i - 1]:  # segment change
                    SScount[currentSS] += 1  # increment segment index
            SSsegments.append(f"{currentSS}{SScount[currentSS]}")

        return SSsegments

    def pdb_to_datafram_with_dssp(self, pdbpath):
        """
        Convert a PDB into a dataframe and add secondary structure info with DSSP
        Args:
            pdbpath (str): path to the PDB file
        Returns:
            df (pd.DataFrame): pdbfile as a dataframe with DSSP and SASA data.
        """
        REL_SASA = {"TRP": 285,
                    "PHE": 240,
                    "GLY": 104,
                    "ALA": 129,
                    "VAL": 174,
                    "ILE": 197,
                    "LEU": 201,
                    "MET": 224,
                    "PRO": 159,
                    "TYR": 263,
                    "SER": 155,
                    "THR": 172,
                    "ASN": 195,
                    "GLN": 255,
                    "CYS": 167,
                    "LYS": 236,
                    "ARG": 274,
                    "HIS": 224,
                    "ASP": 193,
                    "GLU": 223}

        pdbdf = PandasPdb().read_pdb(pdbpath).df["ATOM"]
        #pdbdf.drop_duplicates
        
        """
        
        #remove duplicated residues.
        pdbdf = pdbdf.drop_duplicates(subset=["atom_name","chain_id","residue_number"])

        #Remove the last residue when sometimes you just have the "N"
        lastresid = pdbdf.residue_number.values[-1]
        if len(pdbdf.query('residue_number == @lastresid')) == 1:
            pdbdf = pdbdf.query("residue_number != @lastresid")
            
        """

        ########################
        # calculate dssp and SASA
        ########################

        # pdbmd = md.load_pdb(pdbpath)
        # ssRes = md.compute_dssp(pdbmd)[0] #DSSP
        # sasa_res_mdtraj = md.shrake_rupley(pdbmd, mode="residue")[0] #sasa per atoms
        # sasa_atom_mdtraj = md.shrake_rupley(pdbmd, mode="atom")[0] #sasa per atoms

        p = PDBParser(QUIET=False)
        structure = p.get_structure("struc", pdbpath)
        # Calcul FreeSasa
        try:
            struc_freesasa = freesasa.Structure(pdbpath)
            sasaFreesasa = freesasa.calc(struc_freesasa)
        except Exception as e:
            print("FreeSASA problem", pdbpath, e)
        

        model = structure[0]
        try:
            dssp = DSSP(model, pdbpath, dssp="mkdssp")
        except Exception as e:
            print(pdbpath, e)

    
        ssResSmooth, sasa_res = zip(*[(self.simplify_SS(dssp[x][2]),
                                       dssp[x][3] * 100) for x in dssp.keys()])


        ssRes = [dssp[x][2] for x in dssp.keys()]
        ids = [r_id[1] for c_id,r_id in dssp.keys()]
        pdbdf = pdbdf[pdbdf["residue_number"].isin(ids)]


        # Structure for protein block (structural alphabet)
        # with io.capture_output(stderr=True, display = False) as capture:
        #    structureForPB = pbx.chains_from_files([pdbpath])

        structureForPB = pbx.structure.PDB.PDB(pdbpath).get_chains()
        # dihedrals = structureForPB.get_phi_psi_angles()
        for chain in structureForPB:
            dihedrals = chain.get_phi_psi_angles()
            proteinBlocks = list(pbx.assign(dihedrals))

        # SMOOTH SECONDARY STRUCTURE (helix has to be at least made of 3.6 amino acids... stuff like this)
        ssResSmooth = self.smooth_ss(ssResSmooth)
        # Since DSSP is per amino acids, this is to "copy" the SS to an atom level
        # But I also use this loop to asign the sasa per amino acids.
        ssAtom = []
        ssFullAtom = []
        pbAtom = []
        sasa_res_atom = []
        ASA_res_freesasa_florian = []  # residue sasa per atoms
        freesasa_atom_atom = []
        RSA_sidechain_freesasa = []
        ASA_mainchain_freesasa = []
        ASA_sidechain_freesasa = []
        RSA_total_freesasa_tien = []
        RSA_sidechain_freesasa_tien = []
        RSA_freesasa_florian = []
        ASA_total_freesasa = []
        #atoms_list = list(structure.get_atoms())
        atom_list_biopython =  list(structure.get_atoms())
        atoms_list = pdbdf.atom_number.values
        #residues_list = list(structure.get_residues())
        residues_list = pdbdf.residue_number.values
        residues_insertion_list = pdbdf.insertion.values

        residues_index = {res:index for index,res in enumerate(pdbdf.residue_number.unique())}
        # residues_index = {
        #     res.get_full_id()[3][1]: index for index, res in enumerate(residues_list)
        # }  # looks like 308: 0, 309: 1, 310: 2, 311: 3.....

        for i in range(len(atoms_list)):
            #residue_number = atom_list_biopython[i].get_parent().get_full_id()[3][1]
            #residue_number = residues_index[i]
            residue_number = residues_list[i]
            resIndex = residues_index[residue_number]
            resName = pdbdf.residue_name.values[i]
            chainID = pdbdf.chain_id.values[i]
            ssAtom.append(ssResSmooth[resIndex])
            ssFullAtom.append(ssRes[resIndex])
            pbAtom.append(proteinBlocks[resIndex])
            sasa_res_atom.append(sasa_res[resIndex])
                                   
            ## Note that FreeSaSa needs the insertion information of residues                       
            residue_sasa_florian = freesasa.selectArea([f"x, resi {residue_number}" + residues_insertion_list[i]], struc_freesasa, sasaFreesasa)["x"]
            ASA_res_freesasa_florian.append(residue_sasa_florian)
            # freesasa_atom_atom.append(sasaFreesasa.atomArea(i))
            if resName == "GLY": #GLY don't have sidechains. We will set the value at 0 (exposition = 0)
                RSA_sidechain_freesasa.append(0)
                RSA_sidechain_freesasa_tien.append(0)
                ASA_sidechain_freesasa.append(0)
            else:
                relativeSideChainFlorian = sasaFreesasa.residueAreas()[chainID][str(residue_number)].relativeSideChain * 100
                sideChain = sasaFreesasa.residueAreas()[chainID][str(residue_number)].sideChain
                relativeSideChainTien = sideChain/ REL_SASA[resName] * 100

                RSA_sidechain_freesasa.append(relativeSideChainFlorian)
                ASA_sidechain_freesasa.append(sideChain)
                RSA_sidechain_freesasa_tien.append(relativeSideChainTien)

            mainChain = sasaFreesasa.residueAreas()[chainID][str(residue_number)].mainChain
            totalsasa = sasaFreesasa.residueAreas()[chainID][str(residue_number)].total
            relative_total_tien = totalsasa / REL_SASA[resName] * 100

            ASA_total_freesasa.append(totalsasa)
            ASA_mainchain_freesasa.append(mainChain)
            RSA_total_freesasa_tien.append(relative_total_tien)

            RSA_freesasa_florian.append(residue_sasa_florian / REL_SASA[resName] * 100)
        # Creating the dataframe from PandasPDB
        df = pdbdf[
            [
                "atom_number",
                "atom_name",
                "residue_name",
                "chain_id",
                "residue_number",
                "x_coord",
                "y_coord",
                "z_coord",
                "occupancy",
                "b_factor",
                # 'segment_id',
                # 'element_symbol',
                # 'charge'
            ]
        ]

        # Adding calculated values (SS & SASA)
        df["sec_struc"] = ssAtom
        df["sec_struc_full"] = ssFullAtom
        df["prot_block"] = pbAtom
        df["sasa_rel_dssp"] = sasa_res_atom
        # df["atom_sasa_freesasa"] = freesasa_atom_atom
        df["ASA_res_freesasa_florian"] = ASA_res_freesasa_florian #ASA fetch florian's style
        df["RSA_freesasa_florian"] = RSA_freesasa_florian #REL ASA florian's Style

        df["ASA_total_freesasa"] = ASA_total_freesasa
        df["ASA_mainchain_freesasa"] = ASA_mainchain_freesasa
        df["ASA_sidechain_freesasa"] = ASA_sidechain_freesasa
        df["RSA_sidechain_freesasa"] = RSA_sidechain_freesasa
        df["RSA_total_freesasa_tien"] = RSA_total_freesasa_tien
        df["RSA_sidechain_freesasa_tien"] = RSA_sidechain_freesasa_tien

        # Adding the SS Segment
        df["sec_struc_segment"] = self.define_SS_segment(ssAtom)

        df = df.round({"sasa_rel_dssp": 2})
        return df

    # df = pdb_to_datafram_with_dssp(pdbtest)
    # instanciate the progressbar
    def clean_pdb(self, pdb_path):
        """
        Clean a PDB file (remove alternative residues, rename terminus atoms...)
        This if for DSSP.
        TODO: use this regex to only keep PDB from CATH "^\S{4}\w\d{2}"
        Args:
            pdb_path (str): path to the PDB Path
        Returns:
            None
        """
        f = open(pdb_path, "r")
        outpath = pdb_path.replace("raw","cleaned")
        change = False

        try:
            lines = f.readlines()
        except UnicodeDecodeError as e:
            print(pdb_path)
            print(e)
        for i in range(len(lines)):
            if lines[i].startswith("ATOM"):
                #REMOVE ALTERNATIV COLUMNS
                if lines[i][16:17] != " ":
                    lines[i] = lines[i][:16] + " " + lines[i][17:]
                    change = True
                #if lines[i][26:27] != " ":
                    #lines[i] = lines[i][:26] + " " + lines[i][27:]
                    #change = True
                #Change 'OXT' terminus oxigen into a regular O, otherwise dssp will not work
                if lines[i][12:16] == "OXT":
                    lines[i] = lines[i][:12] + "O  " + lines[i][16:]
                    change = True
        f.close()
        # If the file needs to be modified, we overwrite it.
        if change:
            f = open(outpath, "w")
            for line in lines:
                f.write(line)
            f.close()
        else:
            if pdb_path != outpath:
                shutil.copy2(pdb_path, outpath)





    def process_pdb_list(self, pdblist, datatype='cath', domname = None, mappingFile = None):
        """
        Main function that will create the dataset
         1. Read the "pdb_chain_uniprot.csv"
         format:
             PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
             101m	 A	      P02185	  1	      154	  0	         153	1	  154
          --> Give the correspondance between PDB (in a chain/residue level) and UniprotID

         2. For each pdb in pdblist (pdblist is the list of the pdb from CATH)
            Get all info (pdbname, chain...)
            Transform the PDB into a dataframe and add secondary structure info
         3. Try to get the UniprotID and add it into the dataframe
         4. return all the results, with a ERRORS list as well (list of PDBs were an error has occured)

        Args:
            pdblist (list(str)): a list of PDB path
        Returns:
            dataset (pd.DataFrame): dataset of output results
            ERRORS (list): list of PDBs were an error has occured.
        """

        if datatype == 'custom':
            pdb_uniprot_mapping = pd.read_csv(f"{mappingFile}", comment="#") 
        else:
            pdb_uniprot_mapping = pd.read_csv(f"{self.CATHFOLDER}/pdb_chain_uniprot.csv",
                                          comment="#")  # taken from https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html



        def process_pdb(pdb_path, pdb_uniprot_mapping, datatype='cath', domname = None):

            # FIRST : REMOVE ALTERNATIV POSITIONS
            #self.clean_pdb(pdb_path)

            splittedPath = pdb_path.split("/")

            if datatype=='cath':
                pdbname = splittedPath[-1]
                domain = splittedPath[-3]
                cathpdbcode = pdbname.split(".")[0]
                pdbcode = cathpdbcode[:4].upper()
                chain = cathpdbcode[4]
                data_type = "cathpdb"
                try:
                    uniprot_acc = \
                        pdb_uniprot_mapping.query("PDB == @pdbcode.lower() and CHAIN == @chain")["SP_PRIMARY"].values[0]
                except:
                    uniprot_acc = None

            elif datatype == 'alfafold':
                pdbname = splittedPath[-1]
                domain = splittedPath[-3]
                cathpdbcode = pdbname.split(".")[0]
                pdbcode = pdbname
                chain = "A"
                uniprot_acc = cathpdbcode
                data_type = 'alfafold'

            elif datatype == 'custom':
                pdbname = splittedPath[-1]
                domain = domname
                cathpdbcode = pdbname.split(".")[0]
                pdbcode = cathpdbcode[:4].upper()
                data_type = 'custom'



            try:
                # Error when MDtraj read some PDB with missing segments. they are cut. So the number of residues are not the same
                # It would be nice to find an alternativ to calculate secondary structures.
                df = self.pdb_to_datafram_with_dssp(pdb_path)

            except (ValueError, IndexError, Exception) as e:
                print(f"{domain} - {cathpdbcode} - {e}")
                self.ERRORS.append(f"{domain} - {cathpdbcode} - {e}")
                return None

            if datatype == 'custom':
                chain = df["chain_id"].unique()[0]
                try:
                    uniprot_acc = pdb_uniprot_mapping.query("PDB == @pdbcode.lower() and CHAIN == @chain")["SP_PRIMARY"].values[0]
                except:
                    uniprot_acc = None
            
            df["pdb"] = pdbcode
            df["domain"] = domain
            df["cathpdb"] = cathpdbcode
            df["chain"] = chain
            # Add accesssible surface aera
            df["uniprot_acc"] = uniprot_acc
            df["data_type"] = data_type

            return df

        pdb_list_series = pd.Series(pdblist)

        if self.PARALLEL and not pdb_list_series.empty:
            print("parralel processing", pdb_list_series)
            dataset_list_dataframe = pdb_list_series.parallel_apply(lambda x: process_pdb(x, pdb_uniprot_mapping, datatype, domname=domname))
        else:
            dataset_list_dataframe = pdb_list_series.progress_apply(lambda x: process_pdb(x,pdb_uniprot_mapping, datatype, domname=domname))


        dataset_list_dataframe = dataset_list_dataframe.tolist()

        print(dataset_list_dataframe)
        try:
            returndf = pd.concat(dataset_list_dataframe)
            return returndf
        except ValueError as e:
            print(e)
            print("no update to make, returning None")
            return None


    def add_experimental_method(self, DATASET):
        print("fetching experimental method from PDBe")

        pdb_list = [x.lower() for x in DATASET.pdb.unique()]

        SERVER = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/"
        exp_method = {}

        for pdb in self.tqdm(pdb_list):
            if len(pdb) == 4:
                try:  # Try to get the PDB
                    json_results = requests.get(SERVER + pdb).json()
                except:
                    print("No pdb called", pdb, file=sys.stderr)
                    json_results = None

                try:
                    exp = json_results[pdb][0]['experimental_method']
                except:
                    exp = "unkown"
            else:
                exp="AFmodel"

            exp_method[pdb] = exp



        print("Adding it into dataset")
        if self.PARALLEL:
            DATASET["Experimental Method"] = DATASET.parallel_apply(lambda x: exp_method[x.pdb.lower()], axis=1)
        else:
            DATASET["Experimental Method"] = DATASET.progress_apply(lambda x: exp_method[x.pdb.lower()], axis=1)

        return DATASET

    def clean_all_pdbs(self, custom_mode = False):
        """
        custom_mode is a boolean to specify if we are in a custom mode. 
        """
        if not self.RECALCULATE:
            return None

        pdblist = [str(x) for x in Path(self.CATHFOLDER).glob("domains/**/raw/*.pdb")]
        pdblist_alfafold = ([str(x) for x in Path(self.ALFAFOLDFOLDER).glob("**/extracted/*.pdb")])  # impl:AFS

        allpdbs = pdblist + pdblist_alfafold

        for pdb in self.tqdm(allpdbs, desc='cleaning'):
            self.clean_pdb(pdb)




    def build_structural_dataset(self, checkpoint_name = "checkpoint_structure"):
        updatemode = False

        if not self.RECALCULATE:
            if os.path.exists(f"{self.WORKDIR}/{checkpoint_name}.pkl"):
                print("> Reading checkpoint 1")
                DATASET = pd.read_pickle(f"{self.WORKDIR}/checkpoint_structure.pkl")
            else:
                print("> no checkpoint file, recalculating")

        else:



            pdblist = [str(x) for x in Path(self.CATHFOLDER).glob("domains/**/cleaned/*.pdb")]
            pdblist_alfafold = ([str(x) for x in Path(self.ALFAFOLDFOLDER).glob("**/extracted/*.pdb")]) #impl:AFS

            if self.UPDATE and os.path.isfile(f"{self.WORKDIR}/checkpoint_structure.pkl"):
                print("update mode, loading checkpoint backup file")
                DATASET_backup = pd.read_pickle(f"{self.WORKDIR}/checkpoint_structure.pkl")
                pdb_already_in_db =DATASET_backup["cathpdb"].unique()
                print(f"{len(pdb_already_in_db)} cathpdb found in backup")
                #Removing pdbs already present in the backup file

                #TODO, remove element in common in two lists
                def remove_duplicated_element(pdb_already_processed, pdbpaths):
                    """
                    Remove pdbs that were already processed in the backup file.
                    Args:
                        pdb_already_processed: list of PDBIDS ('2da0A00') already processed.
                        pdbpaths: list of new PDBPATH ('/users/XXX/Y/Z/2da0A00.pdb')
                    Returns:
                        list of new pathpaths

                    """
                    returnlist = []
                    for pdbpath in pdbpaths:
                        pdbname = Path(pdbpath).stem
                        if pdbname not in pdb_already_processed:
                            returnlist.append(pdbpath)
                    return returnlist

                print("removing duplicated pdbs (CATH)")
                beforesize = len(pdblist)
                pdblist = remove_duplicated_element(pdb_already_in_db, pdblist)
                print (f"-----{beforesize - (beforesize - len(pdblist))} new pdbs")



                print("removing duplicated pdbs (ALFAFOLD)")
                beforesize = len(pdblist_alfafold)
                pdblist_alfafold = remove_duplicated_element(pdb_already_in_db, pdblist_alfafold)
                print(f"        {beforesize - (len(pdblist_alfafold)-beforesize)} new pdbs")
                updatemode = True


            #Process alfafold structures
            print(".. Processing alfafold structures ..")
            DATASET_AF = self.process_pdb_list(pdblist_alfafold, datatype='alfafold')

            #process cath structures
            print(".. Processing cath structures ..")
            DATASET_cath = self.process_pdb_list(pdblist, datatype='cath')





            print(f"{len(self.ERRORS)} errors")
            DATASET = pd.concat([DATASET_cath,DATASET_AF])
            if updatemode == True:
                DATASET = pd.concat([DATASET_backup, DATASET])

            DATASET = DATASET.reset_index(drop=True)


            #DATASET["data_type"] = "cathpdb"
            if len(self.ERRORS) > 0:
                print(f"{len(self.ERRORS)} pdbs wasn't processed. Check 'ERRORS' list.")
                print(f"in {self.PEPRMINT_FOLDER}/pdb_processing_errors.log")
                with open(f"{self.PEPRMINT_FOLDER}/pdb_processing_errors.log",'w') as errorlog:
                    for pdb in self.ERRORS:
                        errorlog.write(pdb+'\n')

            DATASET = self.add_experimental_method(DATASET)
            self.save_checkpoint_dataset(DATASET, 'checkpoint_structure')


        return(DATASET)


    def add_protrusions(self,
                        DATASET,
                        atom_selection=["CA", "CB"],
                        sel_lowdens=["CA", "CB"],
                        DISTANCECUTOFF=10,
                        LOWDENSITYTHRESHOLD=25,
                        CI_DISTANCECUTOFF=15,
                        ):

        """
        Get all atoms from a PDB groupby object, filter by [atom_selection], compute the comvexhull based on
        DISTANCECUTOFF, protrusions based on LOWDENSITYTHRESHOLD and CI_DISTANCECUTOFF.

        Args:
            DATASET (pandas.DataFrame): Full Dataset object
            atom_selection (int): Atoms for convexhull calculation
            sel_lowdens (list(str)): Atoms selection for low_density calculation
            DISTANCECUTOFF (int): Distance cutoff for low_density calculation
            LOWDENSITYTHRESHOLD (int): neighborood threshold for low_dens_calculation
            CI_DISTANCECUTOFF (int): cutoff distance for co-insertable.

        Returns:

            DATASET (pandas.DataFrame): Dataset news properties as columns

        """
        print("... Conputing protrusions .. ")

        def calc_protrusions_on_group(
                pdbdata,
                atom_selection=["CA", "CB"],
                sel_lowdens=["CA", "CB"],
                DISTANCECUTOFF=10,
                LOWDENSITYTHRESHOLD=25,
                CI_DISTANCECUTOFF=15,
        ):

            HYDROPHOBICS = ['LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'CYS', 'MET']

            #####################################
            # PROTRUSION CALCULATION
            #####################################
            # instantiate the convexhull flag value
            N = len(pdbdata)

            pdbdata["convhull_vertex"] = 0
            # Taking only atoms selected
            subsel = pdbdata.query("atom_name in @atom_selection")

            # keeping 3D coordinate like [[x,y,z],[x,y,z],....] of the SUBSELECTION
            coords = subsel[["x_coord", "y_coord", "z_coord"]].values

            # Calculating convexhull on the SUBSELECTION
            hull = ConvexHull(coords)

            # Changing convexhull flag (0 the atom is not a vertex, 1 it is) IN THE SUBSELECTION
            # this step is needed because the convexhull indexes match the subselection indexes (with iloc).
            # BUT subsel.index match the original subdataset index
            subsel.iloc[hull.vertices, subsel.columns.get_loc("convhull_vertex")] = 1
            # now we can change the values in the ORIGINAL GROUPED DATASET
            pdbdata["convhull_vertex"][subsel.query("convhull_vertex == 1").index] = 1
            # pdbdata["convhull_vertex"] = convhull_vertex

            #####################################
            # LOW DENSITY PROTRUSION CALCULATION (AND NEIGHBOURS)
            #####################################

            # original implementation
            subset = pdbdata.query(
                "atom_name in @sel_lowdens")  # carefull, of matches if we want to change this. TODO : work only with indexes.
            subset["density"] = 0
            subset["protrusion"] = 0

            # to speed up calculation, save the data in a list to access them faster with numpy ;)
            Nsubset = len(subset)
            columnNeighbours = [np.nan] * Nsubset
            columnNeighboursFull = [np.nan] * Nsubset
            columnDensity = np.repeat(0, Nsubset)
            columnLowDens = np.repeat(0, Nsubset)
            columnHydroProt = np.repeat(0, Nsubset)
            columnIdx = np.repeat(-1, Nsubset)
            columnCoInsertable = np.repeat(0, Nsubset)
            columnsCoInsertableNeighbors = [[]] * Nsubset

            convhullinfo = subset["convhull_vertex"].values
            resnameinfo = subset["residue_name"].values
            resnumber = subset["residue_number"].values
            atomnameinfo = subset["atom_name"].values

            distVect = scidist.pdist(subset[["x_coord", "y_coord", "z_coord"]].values)
            distmat = scidist.squareform(distVect)

            result = []  # to adapt the co_insertable protrusion calculation

            # low density atoms on the convex hull
            lowdens_index = []
            for i, row in enumerate(distmat):
                ######################################
                # NEIGHBOURS
                ######################################
                if atomnameinfo[i] == "CB":
                    density = row[row < DISTANCECUTOFF].shape[0] - 1
                    # calculate density for convexhull verteces only
                    if convhullinfo[i] == 1:
                        # close = distance of atoms below a certain cutoff
                        close = np.where(row < DISTANCECUTOFF)[0]

                        # HERE WE WILL SEARCH FOR CLOSE RESIDUE, TO MAKE STAT ON THEM
                        neighbours = []
                        neighbours_FULL = []
                        current_resnumber = resnumber[i]
                        for neighbor in close:  # TODO OPTIMIZE
                            if not neighbor == i:  # ignore current atoms
                                neighborID = resnumber[neighbor]
                                neighborName = resnameinfo[neighbor]
                                if neighborID not in neighbours: #If the neighbor not in the list
                                    if current_resnumber != neighborID:
                                        neighbours.append(neighborID)  # add the resname in the list
                                        neighbours_FULL.append(f"{neighborName}-{neighborID}")  # add the resname in the list


                        neighboursIDAsString = ';'.join(list(map(str, neighbours)))
                        neighboursFULLAsString = ';'.join(neighbours_FULL)
                        columnNeighbours[i] = neighboursIDAsString
                        columnNeighboursFull[i] = neighboursFULLAsString
                        columnDensity[i] = density

                        # if the number of close atoms are bellow a certain threshold
                        if density <= LOWDENSITYTHRESHOLD:
                            # we consider vertex as a low density protrusion.
                            columnLowDens[i] = 1
                            if resnameinfo[i] in HYDROPHOBICS:
                                columnHydroProt[i] = 1
                    columnIdx[i] = i

            subset["neighboursID"] = columnNeighbours
            subset["neighboursList"] = columnNeighboursFull
            subset["density"] = columnDensity
            subset["protrusion"] = columnLowDens
            subset["is_hydrophobic_protrusion"] = columnHydroProt
            subset["idx"] = columnIdx  # To keep the idx object for the distance matrix
            subset.iloc[lowdens_index, subset.columns.get_loc("protrusion")] = 1

            # COINSERTABLE
            simplices = pd.DataFrame(hull.simplices, columns=["v1", "v2", "v3"])  # Vertex 1, vertex2, vertex3
            protrusions = subset.query("is_hydrophobic_protrusion == 1")
            co_insertable = []
            for residue in protrusions.iterrows():
                # Get 0-based residue index (same used for convexhull)
                idx = residue[1]["idx"]

                # Get all triangles where the current residue is found
                dd = simplices.query(" (v1 == " + str(idx) + " or v2 == " + str(idx) + " or v3 == " + str(idx) + ")")
                # remove duplicates, we just want which vertex our residue is attached to.
                uniques = np.unique(dd.values)
                # remove itself
                uniques = uniques[uniques != idx]
                # Look on protrusions

                # add it to the coinsertable list.
                uniques = [x for x in uniques if x in protrusions['idx'].to_numpy()]

                # Last check, for long edges, it has to be bellow our cutoff distance
                uniques = [x for x in uniques if distmat[idx][x] < CI_DISTANCECUTOFF]

                # Add the co-insertable neighbors
                columnsCoInsertableNeighbors[idx] = [resnumber[x] for x in uniques]

                co_insertable.extend(uniques)
            # remove duplicates
            co_insertable = np.unique(co_insertable)

            # Collect data
            if len(co_insertable) > 0:
                columnCoInsertable[co_insertable] = 1

            subset["is_co_insertable"] = columnCoInsertable
            subset["co_insertable_neighbors"] = columnsCoInsertableNeighbors
            # now we can change the values in the original dataset

            newcolumns = subset.columns.difference(pdbdata.columns).drop("idx")

            # t = subset[["protrusion", "neighbours","density"]]
            pdbdata = pd.concat([pdbdata, subset[newcolumns]], axis=1)


            #Calc lowest density co-insertable
            pdbdata['LDCI'] = False
            try:
                minIndex = pdbdata.query("is_co_insertable == 1").density.idxmin()
                pdbdata.at[minIndex, "LDCI"] = True
            except:
                pass


            return pdbdata

        # test = DATASET.query("cathpdb == '1rlwA00'")
        # %time p = calc_protrusions_on_group(test)

        if self.PARALLEL:  # PARRALEL:
            DATASET = DATASET.groupby(["domain","cathpdb"]).parallel_apply(lambda x: calc_protrusions_on_group(x,
                                                                                                    atom_selection,
                                                                                                    sel_lowdens,
                                                                                                    DISTANCECUTOFF,
                                                                                                    LOWDENSITYTHRESHOLD,
                                                                                                    CI_DISTANCECUTOFF
                                                                                                    )
                                                                )
        else:
            DATASET = DATASET.groupby("cathpdb").progress_apply(lambda x: calc_protrusions_on_group(x,
                                                                                                    atom_selection,
                                                                                                    sel_lowdens,
                                                                                                    DISTANCECUTOFF,
                                                                                                    LOWDENSITYTHRESHOLD,
                                                                                                    CI_DISTANCECUTOFF
                                                                                                    )
                                                                )
        DATASET = DATASET.reset_index(drop=True)
        return DATASET



    def add_structural_cluster_info(self, DATASET):
        domfile = "cath-domain-list.txt"
        colomnDomFile = [
            "cathpdb",
            "Class",
            "Architecture",
            "Topology",
            "Homologous",
            "S35",
            "S60",
            "S95",
            "S100",
            "S100Count",
            "DomSize",
            "resolution",
        ]
        cathDomains = pd.read_csv(self.CATHFOLDER + domfile, comment="#", sep=r"\s+", header=None)
        cathDomains.columns = colomnDomFile
        # Renumber the clusters
        cathDomains["S35"] = cathDomains["S35"].astype(str)
        cathDomains["S60"] = cathDomains["S35"] + "." + cathDomains["S60"].astype(str)
        cathDomains["S95"] = cathDomains["S60"] + "." + cathDomains["S95"].astype(str)
        cathDomains["S100"] = cathDomains["S95"] + "." + cathDomains["S100"].astype(str)

        cathDomains = cathDomains[
            ["cathpdb", "S35", "S60", "S95", "S100", "S100Count", "resolution"]
        ]

        # Mergin with the previous dataset, on cathPDB.
        DATASET_cath= pd.merge(DATASET, cathDomains, on="cathpdb")
        DATASET_af = DATASET.query("data_type == 'alfafold'")
        DATASET = pd.concat([DATASET_cath, DATASET_af])
        return DATASET
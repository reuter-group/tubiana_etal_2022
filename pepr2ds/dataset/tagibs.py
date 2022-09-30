# Imports and settings
import pandas as pd
import numpy as np
import math
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.spatial import ConvexHull
from biopandas.pdb import PandasPdb

from . import widget
from . import analysis


plt.rcParams["figure.dpi"] = 200

# %matplotlib inline
sns.set_style("darkgrid")
import ipywidgets as widgets
from IPython.display import display, Markdown, clear_output
# from tqdm.auto import tqdm
from tqdm.notebook import tnrange, tqdm

from termcolor import colored
from IPython.display import display, HTML
import weasyprint

tqdm.pandas()  # activate tqdm progressbar for pandas apply
pd.options.mode.chained_assignment = (
    None  # default='warn', remove pandas warning when adding a new column
)
pd.set_option("display.max_columns", None)
from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"
# %config InlineBackend.figure_format ='svg' #better quality figure figure
np.seterr(divide='ignore', invalid='ignore')

import MDAnalysis as mda
import nglview as nv

from Bio import PDB
from Bio.PDB import PDBParser
import numpy as np


from pandarallel import pandarallel


class Dataset():
    def __init__(self, dataset, PEPRMINT_FOLDER):
        self.dataset = dataset
        self.haspdb = False
        self.CATHFOLDER = f"{PEPRMINT_FOLDER}/databases/cath/"
        self.FIGURESFOLDER = f"{PEPRMINT_FOLDER}/figures/"
        self.PEPRMINT_FOLDER = PEPRMINT_FOLDER
        #self.ui = widget.Widgets(self)
        self.analysis = analysis.Analysis(self)
        self.ibs = None
        self.nonibs = None
        self.domainDf = None
        self.noAlignment = False
        self.domainLabel = ''

    def get_df_objects(self):
        return self.df_objects




    def load_dataset(self, name, path=None):
        if path == None:
            path = f"{self.PEPRMINT_FOLDER}/dataset"

        self.domainDf = pd.read_pickle(f"{path}/{name}.pkl")
        self.domainLabel = "+".join(self.domainDf.domain.unique())


    def save_dataset(self, name, path=None):
        if path == None:
            path = f"{self.PEPRMINT_FOLDER}/dataset"
        self.domainDf.to_pickle(f"{path}/{name}.pkl")



    #############
    ## DATASET FUNCTIONS
    #############
    def tag_ibs(self, dataset,
                domain,
                pdbreference,
                includeResidueRange=[],
                excludeResidueRange=[],
                extendSS=True,
                withAlignment=False,
                onlyC=False,
                cathCluster=None,
                Uniref=None,
                addSequence=True,
                extendAlign=True,
                excludeStrand=False,
                overide_axis_mode=False,
                zaxis=0,
                extendCoilOnly=True,
                coordinates_folder_name = None, #If coordinate folder is given, all the X,Y,W coordinates will be updated from the PDB inside the folder.
                filter_uniprot_acc = None,
                data_type = None,
                base_folder = 'cath',
                silent=True,
                ):
        """
        TODO
        """

        AATYPE = {
            "LEU": "Hydrophobic,H-non-aromatic",
            "ILE": "Hydrophobic,H-non-aromatic",
            "CYS": "Hydrophobic,H-non-aromatic",
            "MET": "Hydrophobic,H-non-aromatic",
            "TYR": "Hydrophobic,H-aromatic",
            "TRP": "Hydrophobic,H-aromatic",
            "PHE": "Hydrophobic,H-aromatic",
            "HIS": "Positive",
            "LYS": "Positive",
            "ARG": "Positive",
            "ASP": "Negative",
            "GLU": "Negative",
            "VAL": "Non-polar",
            "ALA": "Non-polar",
            "SER": "Polar",
            "ASN": "Polar",
            "GLY": "Non-polar",
            "PRO": "Non-polar",
            "GLN": "Polar",
            "THR": "Polar",
            "UNK": "none"
        }


        # Reset AATYpe if needed... #TODO -> Remove this for final prod, it's only a trick to avoid recalculating the full dataset if we want to change a definition
        if not silent:
            print("Domain=",domain)
        dataset["type"] = dataset.residue_name.apply(lambda x: AATYPE[x])
        dataset.uniprot_acc = dataset.uniprot_acc.astype(str)


        dataset["exposition"] = np.where(dataset['RSA_freesasa_florian'] >= 20,
                                         "exposed",
                                         "buried")
        #Same to "exposed" condition
        dataset["exposed"] = dataset["RSA_freesasa_florian"].apply(lambda x: True if x >= 20 else False)

        #with self.ui.out:
        if not silent:
            print("selecting amino acids")
        #######################
        # Defining borders and alignment position
        #######################
        # Checking
        if not domain in dataset.domain.unique():
            raise ValueError("domain not recognized")

        # Get domain
        df = dataset.query("domain == @domain")
        self.domainLabel = domain


        if data_type == 'cath':
            print("data_type is cathpdb")
            df = df.query("data_type == 'cathpdb'")
        elif data_type == 'alfafold':
            print("data_type is alfafold")
            df = df.query("data_type == 'alfafold'")
        elif data_type == 'cath+af':
            print("data_type is cath+af")
            df = df.query("data_type in ['cathpdb','alfafold']")
        else:
            raise("Problem with datatype")


        df["matchIndex"] = list(range(len(df)))

        #If SH2, clean with CHO data

        from sys import platform
        if platform == "linux" or platform == "linux2":
            sh2_cho = "/mnt/g/clouds/OneDrive - University of Bergen/projects/peprmint/data/Cho_SH2_transformed.xlsx"
        else:
            sh2_cho = "/Users/thibault/OneDrive - University of Bergen/projects/peprmint/data/Cho_SH2_transformed.xlsx"

        if domain == 'SH2':
            cho = pd.read_excel(
                sh2_cho,
            engine='openpyxl').dropna(subset=["Range"])
            uniprot_cho = list(cho.uniprot_acc)
            uniprot_dataset = list(
                df.query("domain == 'SH2' and atom_name == 'CA' and alignment_position == 0").uniprot_acc.unique())
            common = list(set(uniprot_cho).intersection(uniprot_dataset))
            df = df.query("uniprot_acc in @common")


        # CLUSTER REDUNDANCY.
        if cathCluster and Uniref:
            df = self.selectUniquePerCluster(df, cathCluster, Uniref, withAlignment, pdbreference)

        # KEEP ONLY MATCH WITH UNIPROT_ACC
        if filter_uniprot_acc:
            def select_only_one(group):
                keep = group.cathpdb.unique()[0]
                return group.query('cathpdb == @keep')

            df = df.query("uniprot_acc in @filter_uniprot_acc")
            df = df.groupby('uniprot_acc', as_index=False).apply(lambda group: select_only_one(group))


        # check if several borders given, the format should be a list
        if not any(isinstance(i, list) for i in includeResidueRange):
            includeResidueRange = [includeResidueRange]

        if not overide_axis_mode:
            if extendAlign:
                includeAliRange = []
                for s, e in includeResidueRange:
                    start_ali = df.query(
                        "cathpdb == @pdbreference and atom_name == 'CA' and residue_number == @s").alignment_position.values[
                        0]
                    end_ali = df.query(
                        "cathpdb == @pdbreference and atom_name == 'CA' and residue_number == @e").alignment_position.values[
                        0]
                    includeAliRange.append([start_ali, end_ali])
                SelectionString = ' or '.join(
                    ["{} <= alignment_position <= {}".format(s, e) for s, e in includeAliRange])
            else:
                SelectionString = ' or '.join(
                    ["{} <= residue_number <= {}".format(s, e) for s, e in includeResidueRange])

            # Exclusion
            if len(excludeResidueRange) > 0:
                if not any(isinstance(i, list) for i in excludeResidueRange):
                    excludeResidueRange = [excludeResidueRange]
                SelectionString = '(' + SelectionString + ') and not ( ' + ' or '.join(
                    ["{} <= residue_number <= {}".format(s, e) for s, e in excludeResidueRange]) + ')'

            # if extendSS:
            #    ssSegments = df.query("cathpdb == @pdbreference and atom_name == 'CA' and ({0})".format(SelectionString)).sec_struc_segment.unique()
            #    alignmentPos = df.query("cathpdb == @pdbreference and atom_name == 'CA' and sec_struc_segment in @ssSegments").alignment_position
            # else:
            #    alignmentPos = df.query("cathpdb ==  @pdbreference and atom_name == 'CA' and ({0})".format(SelectionString)).alignment_position

            # Change IBS status to True
            df.loc[df.eval(f"{SelectionString}", engine='python'), "IBS"] = True
            df.loc[~df.eval(f"{SelectionString}", engine='python'), "IBS"] = False
        elif overide_axis_mode == True:
            self.noAlignment = True

            if not coordinates_folder_name is None:
                if base_folder == 'cath':
                    base_folder = f"{self.CATHFOLDER}/domains"
                else:
                    base_folder = f"{self.PEPRMINT_FOLDER}/databases/{base_folder}"

                if data_type == 'custom':
                    coordinates_folder = coordinates_folder_name
                else:
                    coordinates_folder = f"{base_folder}/{domain}/{coordinates_folder_name}"
                    
                if not silent:
                    print(coordinates_folder)
                if os.path.exists(coordinates_folder):
                    if not silent:
                        print("UPDATING COORDINATES")

                    def update_coords(group, coordinates_folder):
                        pdb = group.cathpdb.unique()[0]

                        if not os.path.isfile(f"{coordinates_folder}/{pdb}.pdb"):
                            return None
                        pl1 = PandasPdb().read_pdb(f"{coordinates_folder}/{pdb}.pdb").df["ATOM"]
                        # When sometimes we remove duplicated residues we have to be sure that we update
                        # the residues list we take from the new coordinates files
                        residues_number_list = group.residue_name.unique()
                        pl1 = pl1.query("residue_name in @residues_number_list")

                        pdbdf = pl1[["atom_name", "residue_number", "x_coord", "y_coord", "z_coord"]]


                        _merged = group.merge(pdbdf, on=["atom_name", "residue_number"], how="left")

                        # Sometimes there is duplicated atoms, we just keep the first one by removing duplicates.
                        if len(_merged["x_coord_y"].values) != len(group["x_coord"]):
                            _merged = _merged.drop_duplicates(subset=["atom_name", "residue_number"])

                        try:
                            group["x_coord"] = _merged["x_coord_y"].values
                            group["y_coord"] = _merged["y_coord_y"].values
                            group["z_coord"] = _merged["z_coord_y"].values
                        except:
                            print(group)
                            1/0

                        return (group)
                if not silent:
                    df = df.groupby("cathpdb", as_index=False).progress_apply(
                        lambda x: update_coords(x, coordinates_folder))
                else:
                    df = df.groupby("cathpdb", as_index=False).apply(
                        lambda x: update_coords(x, coordinates_folder))

            df.loc[df.eval("z_coord <= @zaxis ", engine='python'), "IBS"] = True
            df.loc[df.eval("z_coord > @zaxis", engine='python'), "IBS"] = False






        def tag_MLIP(df):
            # TODO: Improve this and put it directly in the dataset generation and calculate the real MLIP
            df['LDCI'] = False
            try:
                minIndex = df.query('is_co_insertable == 1').density.idxmin()
                df.at[minIndex, "LDCI"] = True
            except:
                return df

            return df

        if "LDCI" not in df.columns:
            if not silent:
                print("taggin MLIP")
            if not silent:
                df = df.groupby("cathpdb").progress_apply(tag_MLIP)
            else:
                df = df.groupby("cathpdb").apply(tag_MLIP)

        self.domainDf = df


        if not silent:
            print("taggin IBS")
            ibs_nonibs = df.groupby('cathpdb').progress_apply(
                lambda x: self.get_ibs_and_non_ibs(x, extendSS, onlyC, excludeStrand,extendCoilOnly))
        else:
            ibs_nonibs = df.groupby('cathpdb').apply(
                lambda x: self.get_ibs_and_non_ibs(x, extendSS, onlyC, excludeStrand,extendCoilOnly))

        ibs = pd.concat([x[0] for x in ibs_nonibs]).reset_index(drop=True)
        nonibs = pd.concat([x[1] for x in ibs_nonibs]).reset_index(drop=True)

        if not silent:
            print(f"len IBS {len(ibs.cathpdb.unique())}")
            print(f"len nonIBS {len(nonibs.cathpdb.unique())}")

        if addSequence:
            if not silent:
                print("adding sequences")
            ibsSeq = df.query(f"data_type == 'prosite' and ({SelectionString})", engine='python')
            nonIbsSeq = df.query(f"data_type == 'prosite' and not ({SelectionString})", engine='python')
            ibs = pd.concat([ibs, ibsSeq])
            nonibs = pd.concat([nonibs, nonIbsSeq])

        self.ibs = ibs
        self.nonibs = nonibs

        #Update "domainDf" with new tag

        matchIndexIBS = ibs.matchIndex
        matchIndexnonIBS = nonibs.matchIndex

        self.domainDf.loc[self.domainDf.matchIndex.isin(matchIndexIBS), "IBS"] = True
        self.domainDf.loc[self.domainDf.matchIndex.isin(matchIndexnonIBS), "IBS"] = False


        #remove NaN Values
        self.domainDf = self.domainDf.dropna(subset=["residue_name"])
        #return (ibs, nonibs)

    def selectUniquePerCluster(self, df, cathCluster, Uniref, withAlignment=True, pdbreference=None,
                               removeStrand=False):
        """
        Return a datasert with only 1 data per choosed clusters.
        """
        
        # Prints for dubugging
        
        #print("Uniref: ", Uniref)
        #print("Cathcluster: ", cathCluster)
        #print(pdbreference)
        
        if cathCluster not in ["S35", "S60", "S95", "S100"]:
            raise ValueError('CathCluster given not in ["S35","S60","S95","S100"]')

        if Uniref not in ["uniref50", "uniref90", "uniref100"]:
            raise ValueError('CathCluster given not in ["uniref50","uniref90","uniref100"]')

        if withAlignment:
            df = df[~df.alignment_position.isnull()]

        cathdf = df.query("data_type == 'cathpdb'")
        seqdf = df.query("data_type == 'prosite' or data_type == 'alfafold'")

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

        dfReprCathNames = list(cathdf.groupby(["domain", cathCluster]).apply(selectUniqueCath).to_numpy())
        
        
        # df.to_csv('tmp_file.csv') Only for debugging purpose
        if len(dfReprCathNames) > 0:
            excludeUniref = df.query(
                "cathpdb in @dfReprCathNames").uniprot_acc.unique()  # Structures are prior to sequences.
            dfReprUnirefNames = list(seqdf.groupby(["domain", Uniref]).apply(selectUniqueUniref,
                                                                        exclusion=excludeUniref).to_numpy())


        else:
            dfReprUnirefNames = list(seqdf.groupby(["domain", Uniref]).apply(selectUniqueUniref,
                                                                        exclusion = []).to_numpy())



        dfReprCath = cathdf.query("cathpdb in @dfReprCathNames")
        uniproc_acc_cath = dfReprCath.uniprot_acc.unique()
        dfReprUniref = seqdf.query("uniprot_acc in @dfReprUnirefNames")

        return (pd.concat([dfReprCath, dfReprUniref]))

    def get_ibs_and_non_ibs(self, cathpdb, extendSS=True, onlyC=False, excludeStrand=False, extendCoilOnly=False):
        # get secondary structure loops
        # _ssSegments = cathpdb.query("@start <= alignment_position <= @end").sec_struc_segment.unique()
        _ssSegments = list(map(str,cathpdb.query("IBS == True").sec_struc_segment.unique()))
        # check if segment is a loop

        if onlyC:
            ssSegments = []
            for segment in _ssSegments:
                if segment.startswith('C'):
                    ssSegments.append(segment)

            ibs = cathpdb.query("sec_struc_segment in @ssSegments")
            nonibs = cathpdb.query("sec_struc_segment not in @ssSegments")

        else:
            # ibs = cathpdb.loc[cathpdb.eval("@start <= alignment_position <= @end")]
            # nonibs = cathpdb.loc[~cathpdb.eval("@start <= alignment_position <= @end")]
            if extendSS:
                ssSegs = cathpdb.query("IBS == True").sec_struc_segment.unique()

                #TODO: EXTEND COIL ONLY

                if extendCoilOnly:
                    ssSegs = [x for x in ssSegs if x.startswith('C')]
                    ibs = cathpdb.query("IBS == True or sec_struc_segment in @ssSegs")
                    nonibs = cathpdb.query("IBS == False and sec_struc_segment not in @ssSegs")
                    ibs.IBS = True
                    nonibs.IBS = False
                    return(ibs,nonibs)



                if not excludeStrand:
                    ibs = cathpdb.query("sec_struc_segment in @ssSegs")
                else:
                    ibs = cathpdb.query("sec_struc_segment in @ssSegs and not sec_struc == 'E'")



                nonibs = cathpdb.query("sec_struc_segment not in @ssSegs")
                # Update IBS
                ibs.IBS == True

            else:
                ibs = cathpdb.query("IBS == True")
                nonibs = cathpdb.query("IBS == False")




        return (ibs, nonibs)

    def generate_picutre_of_IBS(self, subfolder='raw'):
        import MDAnalysis as mda

        cathfolder = f"{self.PEPRMINT_FOLDER}/databases/cath"
        domain = self.ibs.domain.unique()[0]
        pdbfolder = f"{cathfolder}/domains/{domain}/{subfolder}"
        outputPDB = f"{cathfolder}/domains/{domain}/IBS/pdb"
        outputPNG = f"{cathfolder}/domains/{domain}/IBS/png"

        view = {
            "PH": "set_view (0.579240322,    0.505582690,   -0.639426887,    -0.760016978,    0.618558824,   -0.199398085,     0.294711024,    0.601475298,    0.742546558,    -0.000000262,    0.000000060, -192.327224731,     1.087764502,   -0.489078194,    0.073249102,   155.025390625,  229.629058838,  -20.000000000 )",
            'C2': "set_view (0.732501030,    0.379622459,   -0.565089643, -0.677367389,    0.323625565,   -0.660633683, -0.067914240,    0.866690934,    0.494200438, -0.000011362,   -0.000007764, -146.056518555,  4.814969063,   -0.570566535,   10.014616013,115.152069092,  176.960968018,  -20.000000000 )",
            'START':"set_view (     0.997152150,   -0.059061568,    0.046864305,     0.044335004,   -0.043418523,   -0.998070359,     0.060982693,    0.997307599,   -0.040675215,     0.000018209,   -0.000058129, -155.363037109,   -10.252207756,    6.747013569,   12.991518974,   115.271072388,  195.459533691,  -20.000000000 )",
            'C1':"set_view (     0.313648969,    0.817849040,    0.482437909,    -0.183394849,    0.550686538,   -0.814315557,    -0.931659758,    0.166932240,    0.322710901,     0.000000011,    0.000000708, -117.470191956,    -2.323507547,    1.068632245,    1.322908878,    99.383796692,  135.557312012,  -20.000000000 )",
            'SH2':"set_view (     0.552586615,   -0.076772571,    0.829912364,     0.832381070,    0.101391248,   -0.544850767,    -0.042316202,    0.991879404,    0.119931221,     0.000000000,    0.000000000, -132.395614624,    -2.113309860,    3.920848846,    8.741382599,   109.715965271,  155.075241089,  -20.000000000 )",
            "C2DIS":"set_view (     0.143325791,    0.803093255,   -0.578357637,     0.467322886,    0.460218340,    0.754856706,     0.872391641,   -0.378470331,   -0.309342563,     0.000000000,    0.000000000, -154.193801880,     0.543474197,   -0.320308685,    1.579742432,   121.567565918,  186.820037842,  -20.000000000 )",
            "FYVE":"set_view (     0.886591494,    0.351747036,   -0.300380319,    -0.339166492,    0.052801486,   -0.939242423,    -0.314515173,    0.934602559,    0.166114911,    -0.000005476,    0.000001019, -188.998291016,    -1.669404864,   -6.018759727,    0.102616847,   112.505050659,  265.490844727,  -20.000000000 )",
            "PX":"set_view (     0.807477295,    0.112026922,    0.579163909,     0.589899063,   -0.153209910,   -0.792809069,    -0.000082285,    0.981823087,   -0.189797714,     0.000000000,    0.000000000, -159.604339600,     4.754449844,    6.864978790,    9.598649979,   125.833282471,  193.375396729,  -20.000000000 )",
            "ENTH":"set_view (     0.679595053,    0.116319597,    0.724306464,     0.711486697,   -0.345033497,   -0.612157106,     0.178704709,    0.931353927,   -0.317244768,     0.000000000,    0.000000000, -156.046157837,     0.794780731,   -7.394954681,   11.364852905,   102.870040894,  209.222259521,  -20.000000000 )",
            "PLD":"set_view (    -0.942126930,    0.016740968,    0.334838033,     0.324299693,    0.298773795,    0.897533059,    -0.085015431,    0.954176843,   -0.286911786,     0.000000000,    0.000000000, -155.298629761,    -3.902690887,   -0.973564148,   13.115238190,   124.827445984,  185.769866943,  -20.000000000 )",
            "ANNEXIN":"set_view (    -0.082639754,    0.010267707,    0.996525764,     0.996567369,    0.005868928,    0.082582556,    -0.005001415,    0.999930799,   -0.010716723,     0.000000000,    0.000000000, -107.264595032,    -3.597750664,   -2.025382519,    5.951478481,    86.362419128,  128.166778564,  -20.000000000 )",
            "PLA":"set_view (    -0.986954212,   -0.063351221,   -0.148008540,    -0.151886493,    0.061519083,    0.986482680,    -0.053388733,    0.996094048,   -0.070337638,    -0.000003427,   -0.000011377, -137.629989624,     4.810346603,   -4.336011410,    8.157093048,   108.358093262,  166.902557373,  -20.000000000 )",
        }
        if not os.path.isdir(outputPDB):
            os.makedirs(outputPDB)
        if not os.path.isdir(outputPNG):
            os.makedirs(outputPNG)

        pdblist = self.ibs.cathpdb.unique()
        print(len(pdblist))
        for pdb in tqdm(pdblist):
            # read pdb and tag IBS in the beta-factor column
            file = f"{pdbfolder}/{pdb}.pdb"
            ibsResidue = list(map(int, self.ibs.query("cathpdb == @pdb").residue_number.unique()))
            selectionString = " or ".join([f"resnum {x}" for x in ibsResidue])
            U = mda.Universe(file)
            U.atoms.tempfactors = 0
            U.select_atoms(selectionString).tempfactors = 1
            U.atoms.write(f"{outputPDB}/{pdb}.pdb")

            # generate picture.

            pymolCmd = f"load {outputPDB}/{pdb}.pdb; " + \
                       f"{view[domain]};" + \
                       "as cartoon; select ibs, b > 0.5; color red, ibs;" + \
                       f"bg white; png {outputPNG}/{pdb}.png, ray=0"
            # "select ax, z < 0.0001; show spheres, ax and name CA; set sphere_color, blue; set sphere_scale, 0.4;" + \

            _ = os.system(f'pymol -Q -c -d "{pymolCmd}"')



    # START
    def show_structure_and_plane(self, idpdb, folder='raw'):
        parser = PDB.PDBParser()
        from Bio.PDB.PDBExceptions import PDBConstructionWarning
        import warnings
        import nglview as nv

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = parser.get_structure(id=idpdb,
                                             file=f"/Users/thibault/Documents/WORK/peprmint/databases/cath/domains/{self.domainLabel}/{folder}/{idpdb}.pdb")


        view = nv.show_biopython(structure)
        shape = view.shape
        chain = structure[0].child_list[0].id


        shape.add_sphere([0, 0, 0], [1, 0, 0], 1)
        shape.add_sphere([10, 10, 0], [0, 1, 0], 0.5)
        shape.add_sphere([-10, -10, 0], [0, 1, 0], 0.5)
        shape.add_sphere([10, -10, 0], [0, 1, 0], 0.5)
        shape.add_sphere([-10, 10, 0], [0, 1, 0], 0.5)

        shape.add_arrow([0, 0, 0], [0, 0, -10], [1, 0, 0], 1.0)

        mesh = [20, 20, 0,
                20, -20, 0,
                -20, -20, 0,
                -20, -20, 0,
                20, 20, 0,
                -20, 20, 0]

        color = [[0, 1, 0]] * len(mesh)  # RGB, for now let's fix it to blue
        color = np.asarray(color).flatten().tolist()
        shape.add_mesh(mesh, color)

        view.update_representation(component=7, repr_index=0, opacity=0.5, side="double")
        return view


  

    def export_dataset_PePrMInt(self, outputfile=None):
        """
        Clean and prepare the dataset for PePrMInt database (PrPr2DS)
        Returns: df (pd.DataFrame)
        """
        long2short = {'domain': 'dm',
            'cathpdb': 'cath',
            'pdb': 'pdb',
            'uniprot_acc': 'uacc',
            'uniprot_id': 'uid',
            'residue_name': 'rna',
            'IBS': 'ibs',
            'chain_id': 'chain',
            'residue_number': 'rnu',
            'b_factor': 'bf',
            'sec_struc': 'ss',
            'sec_struc_full': 'ssf',
            'prot_block': 'pb',
            'data_type': 'dt',
            'Experimental Method': 'em',
            'resolution': 'rsl',
            'RSA_total_freesasa_tien': 'rsa',
            'convhull_vertex': 'cv',
            'protrusion': 'pro',
            'is_hydrophobic_protrusion': 'hypro',
            'is_co_insertable': 'coin',
            'neighboursList': 'nbl',
            'density': 'den',
            'exposed': 'expo',
            'S35': 's35',
            'S60': 's60',
            'S95': 's95',
            'S100': 's100',
            'uniref50': 'u50',
            'uniref90': 'u90',
            'uniref100': 'u100',
            'origin': 'origin',
            'location': 'loc',
            'taxon': 'taxon'}



        #tips to save only 1 atom per PDB. take ONLY CB for residue except for GLY and get rid of the atom_name 
        df = self.domainDf.query("(atom_name == 'CB') or (residue_name == 'GLY')")        
        #Correct alphafold typo...
        df['data_type'].replace({'alfafold':'alphafold'}, inplace=True)
        #Actually, remove alphafold for pepr2web
        df = df.query("data_type != 'alphafold'")



        keepColNames = list(long2short.keys())
        #Add the columns if does not exist (example, from custom CSV generation).
        for col in keepColNames:
            if col not in df.columns:
                df[col] = np.nan

        df = df[
            ["domain", "cathpdb", 'pdb',
            'uniprot_acc', 
            'uniprot_id', 
            # 'atom_number', 
            #'atom_name', 
            'residue_name', 
            'IBS', #BINDING SITE 
             'chain_id', 
             'residue_number',
             'b_factor', 'sec_struc', 'sec_struc_full', 'prot_block', # Structure stuff
             'data_type', 'Experimental Method', 'resolution', 'RSA_total_freesasa_tien', # Experimental stuff
             'convhull_vertex', 'protrusion', 'is_hydrophobic_protrusion', 'is_co_insertable', 'neighboursList', 'density', # Hydrophobic protrusions 
             'exposed',  #Exposision
              'S35', 'S60', 'S95', 'S100', 'uniref50', 'uniref90', 'uniref100',  #Clusters
             'origin', 'location', 'taxon',
             ]]



        #Renaming the structures if needed
        df = df.rename({
            'cathpdb': 'structure_name',
            'RSA_total_freesasa_tien':'RSA',
            }
        )




        df.to_csv(outputfile, index=False,) # Full dataset, Full name

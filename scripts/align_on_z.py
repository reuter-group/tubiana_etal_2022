from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from pathlib import Path
import warnings
import numpy as np
import glob
import vg

import sys
import os
from tqdm import tqdm
import argparse



def parseArg():
    """
    This fonction will the list of pdb files and the distance
    @return: dictionnary of arguments
    """
    arguments = argparse.ArgumentParser(prog="CalcAngle Through Trajectory")
    arguments.description = "\
            This program is made to calcul the angle between 2 protein domains\
            through a trajectory. MSF is calculed before spliting domains to \
            remove amino acids with high MSF (amino acids which high \
            fluctuation). MSF/boxplot/angle is produced, as well as a logfile.\
            "

    arguments.add_argument('-d', "--domain", help="domain", required=True)
    arguments.add_argument('-ref', '--ref', help="reference file", required=True)
    arguments.add_argument('-res1', '--res1', help="residue 1", required=True, type=int)
    arguments.add_argument('-res2', '--res2', help="residue 2", required=True, type=int)
    arguments.add_argument('-res3', '--res3', help="residue 3", required=True, type=int)
    arguments.add_argument('-dir', '--directory', help="Directory where Domains folder are", default="cath")
    arguments.add_argument('-i', '--inputfolder', help="Folder where the structures are", default="raw")
    arguments.add_argument('-o', '--outputfolder', help="where to save the structures",default="zaligned")
    args = vars(arguments.parse_args())
    return (args)

def get_structure(path:str):
    """
    Read a PDB file with PDB parser
    Args:
        path: String. Path to the variable

    Returns:

    """
    parser = PDB.PDBParser()
    pdbCode = Path(path).stem #Get PDB code
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = parser.get_structure(id=pdbCode,
                                             file=path)
    except:
        print(f"Reading Error - {path}")
    return structure


def write_structure(output_folder:str, structure:PDB.Structure.Structure):
    """
    Write the Biopython structure
    Note that the 'structure.id' will be used as pdbname.
    Args:
        output_folder: (output pdb path)
        structure (PDB.Structure.Structure): Biopython structure

    Returns:
        None
    """
    pdbCode = structure.id
    outputPath = f"{output_folder}/{pdbCode}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(outputPath)


def get_translation_vector(structure:PDB.Structure.Structure, res1:int, res2:int, res3:int):
    """
    Get the translation vector between the centroid of a triangle formed by 3 amino acids's CA and the origin (0,0,0)
    Args:
        structure (PDB.Structure.Structure): Biopython PDB Structure
        res1 (int): Residue number forming Triangle vertex 1
        res2 (int): Residue number forming Triangle vertex 2
        res3 (int): Residue number forming Triangle vertex 3

    Returns:
        translation (np.array):  translation vector (1x3)

    """
    chain = structure[0].child_list[0].id
    p1 = structure[0][chain][res1]['CA'].get_coord()
    p2 = structure[0][chain][res2]['CA'].get_coord()
    p3 = structure[0][chain][res3]['CA'].get_coord()

    translation = list(-np.mean(np.array([p1, p2, p3]), axis=0))

    return translation


def get_rotation_matrix(structure:PDB.Structure.Structure, res1:int, res2:int, res3:int, domain:str, orientation='z'):
    """
    Get the rotation matrix between the normal of a triangle formed by 3 amino acids's CA and a plane (x,y,z)
    Args:
        structure (PDB.Structure.Structure): Biopython PDB Structure
        res1 (int): Residue number forming Triangle vertex 1
        res2 (int): Residue number forming Triangle vertex 2
        res3 (int): Residue number forming Triangle vertex 3
        orientation (str): Axis used for alignment (default = 'z'

    Returns:
        rotation (np.array):  Rotation matrix (3x3)

    """

    def get_normal_COM(structure, res1, res2, res3):
        """
        Calcul the normal and the geom center of a structure
        Args:
            structure:
            res1:
            res2:
            res3:

        Returns:

        """
        # ROTATION
        # 1. compute vectors
        # Get new coordinates
        chain = structure[0].child_list[0].id
        p1 = structure[0][chain][res1]['CA'].get_coord()
        p2 = structure[0][chain][res2]['CA'].get_coord()
        p3 = structure[0][chain][res3]['CA'].get_coord()

        # Translation = Nul
        A = p2 - p1
        B = p3 - p1
        # 2. compute triangle NORM which is the cross product of vector A/B
        N = np.cross(A, B)

        coords = np.array([x.get_coord() for x in structure.get_atoms()])
        COM = coords.mean(axis=0)

        return N, COM

    def test_rotation(structure, res1, res2, res3):
        # Recalculate angle etc....
        N, COM = get_normal_COM(structure, res1, res2, res3)
        # Recalculate angle
        angle = vg.angle(N, np.array([0, 0, -1]))
        return (angle == 0 and COM[2] > 0) or (angle == 180 and COM[2] > 0)

    N,COM = get_normal_COM(structure, res1, res2, res3)

    # This norm will be our translation vector to all our atoms
    axis = {'x':[-1,0,0],
            'y':[0,-1,0],
            'z':[0,0,-1]}
    #Create the reference vector, per default we want to align on the z axis so it will be [0,0,1]
    refVector = PDB.vectors.Vector(axis[orientation])
    normal = PDB.vectors.Vector(N) #The normal should be a biopython object

    # Transformation 1
    temp = structure.copy()

    #case1 : normal and rotation
    rotation = PDB.vectors.rotmat(normal, refVector)
    temp.transform(rotation, [0, 0, 0])

    #If it doesn't work, case2: -normal and rotation
    if not test_rotation(temp, res1, res2, res3):
        temp = structure.copy()
        rotation = PDB.vectors.rotmat(-normal, refVector)
        temp.transform(rotation, [0, 0, 0])
        # If it doesn't work, case3: normal and rotation.T
        if not test_rotation(temp, res1, res2, res3):
            temp = structure.copy()
            rotation = PDB.vectors.rotmat(normal, refVector).T
            temp.transform(rotation, [0, 0, 0])
            # If it doesn't work, case2: -normal and rotation.T
            if not test_rotation(temp, res1, res2, res3):
                temp = structure.copy()
                rotation = PDB.vectors.rotmat(-normal, refVector).T

    return rotation



def apply_rotation_matrix(structure:PDB.Structure.Structure, rotation:np.array):
    """
    Apply a translation on the structure based on the vector Translation
    Args:
        structure (PDB.Structure.Structure): Biopython PDB Structure
        rotation (np.array):  Rotation matrix (3x3)

    Returns:
        structure (PDB.Structure.Structure): rotated Biopython PDB Structure
    """
    #Rotation without translation
    structure[0].transform(rotation, [0, 0, 0])
    return structure


def apply_translation_vector(structure:PDB.Structure.Structure, translation:np.array):
    """
    Apply a translation on the structure based on the vector Translation
    Args:
        structure (PDB.Structure.Structure): Biopython PDB Structure
        translation (np.array):  Translation vector (1x3)

    Returns:
        structure (PDB.Structure.Structure): Translated Biopython PDB Structure
    """
    rotation = np.identity(3).tolist()
    structure[0].transform(rotation, translation)
    return structure



def get_transformation_from_reference(templateCode:str, pdbFolder:str, res1:int, res2:int, res3:int, domain:str):

    pdbPath = f"{pdbFolder}/{templateCode}.pdb"

    structure =  get_structure(pdbPath)
    #1. Get translation
    translation = get_translation_vector(structure, res1, res2, res3)
    #2. Apply translation
    structure = apply_translation_vector(structure, translation)
    #3. get rotation
    rotation = get_rotation_matrix(structure, res1, res2, res3, domain)

    #It's enough, no need to rotate the reference since we will transform all structures in the next step

    return (rotation,translation)


def create_output_dir(directory_path:str, domain:str, outputFolderName):
    """
    Create the "zaligned" output folder if it not extists yet.
    Args:
        directory_path (str): CATH folder
        domain (str): Domain name

    Returns:
        None
    """
    output_folder = f"{directory_path}/{domain}/{outputFolderName}"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


def transform_pdbs(directory_path:str, domain:str, rotation:np.array, translation:np.array,inputFolderName,outputFolderName):
    """
    Search for all PDB in the domain PDB folder, apply a translation and a rotation and then save the PDB in a new folder (zaligned)
    Args:
        directory_path (str): CATH folder
        domain (str): Domain name
        rotation (np.array): Rotation matrix (3x3)
        translation (np.array): Translation vector (1x3)

    Returns:
        None
    """

    pdbFolder = f"{directory_path}/{domain}/{inputFolderName}"
    outputFolder = f"{directory_path}/{domain}/{outputFolderName}"



    pdb_list = glob.glob(f"{pdbFolder}/*.pdb")

    for pdb in tqdm(pdb_list):
        structure = get_structure(pdb)
        structure = apply_translation_vector(structure,translation)
        structure = apply_rotation_matrix(structure, rotation)
        write_structure(outputFolder, structure)



def set_directory_path(folder):
    """
    Set the directory path, where all domains folders are located.
    Args:
        folder: string. Folder name ("cath" for example).

    Returns:
        directory_path.
    """
    import platform
    system = platform.system()
    if folder == "cath":
        folder = "cath/domains"
    if system == "Linux":
        print(os.getcwd())
        directory_path = os.path.join(os.getcwd(),f"peprmint/databases/{folder}/")
    elif system == "Darwin":
        directory_path = f"/Users/thibault/Documents/WORK/peprmint/databases/{folder}/"
    else:
        print("Path is not configure for windows. exiting")
        sys.exit(1)

    if os.path.exists(directory_path):
        return directory_path
    else:
        print(f"directory {directory_path} doesn't exist. Exiting")
        sys.exit(1)



if __name__ == "__main__":
    #Retrieving Arguments


    args = parseArg()
    domain = args["domain"]
    pdbFileCode = args["ref"]
    res1 = args["res1"]
    res2 = args["res2"]
    res3 = args["res3"]
    directory_path = set_directory_path(args["directory"])


    inputFolderName = args["inputfolder"]
    outputFolderName = args["outputfolder"]


    #Create the output dir that will contain the aligned pdbs
    create_output_dir(directory_path, domain, outputFolderName)


    #Set PDB FOlder
    pdbFolder = f"{directory_path}/{domain}/{inputFolderName}"

    #Get Rotation matrix and translation
    rotation, translation = get_transformation_from_reference(pdbFileCode, pdbFolder,
                                                                      res1,res2,res3,
                                                              domain)

    print(f"rotation matrix: {rotation}")
    print(f"translation vector: {translation}")

    #Now LET'S TRANSFORM ALL PDDDBBBBBB
    transform_pdbs(directory_path, domain, rotation, translation, inputFolderName, outputFolderName)

    print("Done")









from . import Attributes as Attributes
from . import Structure as Structure
from . import Sequence as Sequence
import pandas as pd
import numpy as np
import importlib
importlib.reload(Structure)
importlib.reload(Attributes)
importlib.reload(Sequence)



class Builder(Attributes.Attributes):
    def __init__(self, SETUP: dict, recalculate=True, update=True, notebook=True, core=4):
        super().__init__(SETUP=SETUP, recalculate=recalculate, update=update, notebook=notebook, core=core)
        self.structure = Structure.Structure(SETUP=SETUP, recalculate=recalculate, update=update, notebook=notebook,core=core)
        self.sequence = Sequence.Sequence(SETUP=SETUP, recalculate=recalculate,update=update, notebook=notebook, core=core)


    def optimize_size(self, DATASET):
        print("> datatypes optimisation")
        # Calculate the aatype

        typedf = pd.DataFrame(list(set(self.AATYPE.items())), columns=["residue_name", "type"])
        DATASET = pd.merge(DATASET, typedf, on="residue_name")

        beforeOptimisation = DATASET.memory_usage().sum() / 1024 ** 2

        print("Size BEFORE optimization {0:.2f} MB".format(beforeOptimisation))
        
        DATASET.head()
        DATASET["type"] = DATASET["type"].astype("category")
        DATASET["sec_struc"] = DATASET["sec_struc"].astype("category")
        DATASET["sec_struc_full"] = DATASET["sec_struc_full"].astype("category")

        DATASET["prot_block"] = DATASET["prot_block"].astype("category")
        DATASET["sec_struc_segment"] = DATASET["sec_struc_segment"].astype("category")
        DATASET["protrusion"] = (
            DATASET["protrusion"].fillna(False).astype(bool)
        )
        DATASET["convhull_vertex"] = (
            DATASET["convhull_vertex"].fillna(False).astype(bool)
        )
        DATASET["density"] = DATASET["density"].fillna(0).astype(np.uint16)
        DATASET["S100Count"] = DATASET["S100Count"].fillna(0).astype(np.uint16)
        DATASET["atom_number"] = DATASET["atom_number"].fillna(0).astype(np.uint16)
        DATASET["domain"] = DATASET["domain"].astype("category")
        DATASET["residue_number"] = DATASET["residue_number"].astype(np.int16)
        # DATASET["alignment_position"] = DATASET['alignment_position'].astype(np.int16)
        DATASET["data_type"] = DATASET["data_type"].astype("category")
        DATASET["atom_name"] = DATASET["atom_name"].astype("category")
        DATASET["residue_name"] = DATASET["residue_name"].astype("category")
        DATASET["chain_id"] = DATASET["chain_id"].astype("category")
        DATASET["pdb"] = DATASET["pdb"].astype("category")
        # DATASET["cathpdb"] = DATASET['cathpdb'].astype('category')
        DATASET["S35"] = DATASET["S35"].astype("category")
        DATASET["S60"] = DATASET["S60"].astype("category")
        DATASET["S95"] = DATASET["S95"].astype("category")
        DATASET["S100"] = DATASET["S100"].astype("category")
        DATASET["uniprot_acc"] = DATASET["uniprot_acc"].astype("category")
        DATASET["uniprot_id"] = DATASET["uniprot_id"].astype("category")
        DATASET["origin"] = DATASET["origin"].astype("category")
        DATASET["taxon"] = DATASET["taxon"].astype("category")
        DATASET["uniref50"] = DATASET["uniref50"].astype("category")
        DATASET["uniref90"] = DATASET["uniref90"].astype("category")
        DATASET["uniref100"] = DATASET["uniref100"].astype("category")
        DATASET["prositeName"] = DATASET["prositeName"].astype("category")
        DATASET["prositeID"] = DATASET["prositeID"].astype("category")
        DATASET["residue_index"] = DATASET["residue_index"].astype(np.int16)
        #DATASET["alignment_position"] = DATASET["alignment_position"].astype(pd.Int16Dtype())
        DATASET["is_hydrophobic_protrusion"] = DATASET["is_hydrophobic_protrusion"].fillna(0).astype(bool)
        DATASET["convhull_vertex"] = DATASET["convhull_vertex"].fillna(0).astype(bool)
        DATASET["is_co_insertable"] = DATASET["is_co_insertable"].fillna(0).astype(bool)

        afterOptimisation = DATASET.memory_usage().sum() / 1024 ** 2

        print("Size AFTER optimization {0:.2f} MB".format(afterOptimisation))
        print(
            "{0:.2f}% of the original size".format(afterOptimisation / beforeOptimisation * 100)
        )
        return DATASET





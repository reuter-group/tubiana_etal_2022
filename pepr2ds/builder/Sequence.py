from .Attributes import Attributes
from .Uniprot import *

import pandas as pd
import numpy as np
import requests
import os
import urllib.request, urllib.parse
from collections import defaultdict
import difflib
import re
from Bio import Entrez
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import glob


class Sequence(Attributes):
    def __init__(self, SETUP: dict, recalculate=True, update=True, notebook=True, core=4):
        super().__init__(SETUP, recalculate, update, notebook, core)

        if notebook:
            from tqdm.notebook import tnrange, tqdm
            self.tqdm = tqdm
            self.tnrange = tnrange
            self.tqdm.pandas()  # activate tqdm progressbar for pandas apply
        else:
            print("notebook = False")
            from tqdm import tnrange, tqdm
            self.tqdm.pandas()  # activate tqdm progressbar for pandas apply

    def add_uniprotId_Origin(self, DATASET):
        def get_uniprotID_from_ACC(uniprot_accs):
            """
            Search in uniprot for GENEID from uniprot_acc :
            args:
                pdbs <list> : list of all pdb names
            return:
                mapping <dict> : mapping of {pdb:uniprot}
            """
            
            #### DEPRECATED CODE ####
            """
            url = "https://www.uniprot.org/uploadlists/"

            params = {
                "from": "ACC",
                "to": "ID",
                "format": "tsv",
                "query": ",".join(uniprot_accs[1:5]),
            }
            data = urllib.parse.urlencode(params)
            data = data.encode("utf-8")
            req = urllib.request.Request(url, data)
            print("doodle", data)
            with urllib.request.urlopen(req) as f:
                response = f.read()
                
            """
            
            job_id = submit_id_mapping(from_db = "UniProtKB_AC-ID", to_db = "UniProtKB", ids= uniprot_accs)
            
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                results = get_id_mapping_results_search(link)
                print(results["results"][1:2])
                mapping= {x["from"] : x["to"]["uniProtkbId"] for x in results["results"]}
                return mapping

            # Yes I know... not understandable, but since "response" is a binary text of the results, this is just
            # to convert the result in a dictionnary... :-)
            """
            mapping = {}

            twoByTwo = list(zip(*[iter(response.decode("utf-8").split()[2:])] * 2))

            for uniprot, gene in twoByTwo:
                mapping[uniprot] = gene

            return mapping
            """

        DATASET.fillna(np.NaN, inplace=True)  # replace nan properly
        uniprot_accs = DATASET.uniprot_acc.unique()
        uniprot_accs = [x for x in uniprot_accs if type(x) == type("string")]
        mapDict = get_uniprotID_from_ACC(uniprot_accs)
        # add NAN
        mapDict[np.NaN] = np.NaN
        DATASET["uniprot_id"] = DATASET.uniprot_acc.progress_apply(lambda x: mapDict[x])
        DATASET["origin"] = DATASET["uniprot_id"].str.split("_").str[1]

        return DATASET

    def add_cluster_info(self, DATASET):
        def get_unirefID(uniprotAccs, to="NF90"):
            """
            Search in uniprot for GENEID from uniprot_acc :
            database identification
            https://www.uniprot.org/help/api_idmapping
            args:
                pdbs <list> : list of all pdb names
            return:
                mapping <dict> : mapping of {pdb:uniprot}
               
            """
            #### DEPRACATED CODE ####
            """
            url = "https://www.uniprot.org/uploadlists/"

            params = {
                "from": "ACC+ID",
                "to": to,
                "format": "tab",
                "query": ",".join(uniprotAccs),
            }
            data = urllib.parse.urlencode(params)
            data = data.encode("utf-8")
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read()

            # Yes I know... not understandable, but since "response" is a binary text of the results, this is just
            # to convert the result in a dictionnary... :-)
            mapping = {}

            twoByTwo = list(zip(*[iter(response.decode("utf-8").split()[2:])] * 2))

            for uniprot, uniref in twoByTwo:
                mapping[uniprot] = uniref.split('_')[1]



            return mapping
            """
            job_id = submit_id_mapping(from_db = "UniProtKB_AC-ID", to_db = to,ids=uniprot_accs)
            
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                results = get_id_mapping_results_search(link)
                print(results["results"][1:2])
                mapping= {x["from"] : x["to"]["uniProtkbId"] for x in results["results"]}
                return mapping
            

        def try_fetch_database(uniprotaccs, db, iteration=0):
            import time
            try:
                print("trying to fetch " + db)
                r = get_unirefID(uniprotaccs, db)
                print("  >done")
                return r
            except urllib.error.HTTPError:
                if iteration < 5:
                    print(">HTTP error (propably internal server error). Waiting 5 second and trying again.")
                    time.sleep(5)
                    return try_fetch_database(uniprotaccs, db, iteration + 1)  # Let's try again!
                else:
                    print("Too many tries... Try again later... Maybe use a VPN (UiB or another one)")
                    raise

        uniprotaccs = DATASET.uniprot_acc.dropna().unique()
        uniprotaccs = uniprotaccs[np.logical_not(uniprotaccs == np.nan)]
        # print("> Database queries, if fail, start again later")

        uniref50 = try_fetch_database(uniprotaccs, "NF50")
        uniref90 = try_fetch_database(uniprotaccs, "NF90")
        uniref100 = try_fetch_database(uniprotaccs, "NF100")

        print("> mapping with dataset")
        DATASET["uniref50"] = DATASET["uniprot_acc"].apply(
            lambda acc: uniref50[acc] if acc in uniref50 else "obsolete")
        DATASET["uniref90"] = DATASET["uniprot_acc"].apply(
            lambda acc: uniref90[acc] if acc in uniref90 else "obsolete")
        DATASET["uniref100"] = DATASET["uniprot_acc"].apply(
            lambda acc: uniref100[acc] if acc in uniref100 else "obsolete")
        print("  >ok<  ")
        return (DATASET)

    def return_aligned_seq_as_dict(self, alignmentfiles):
        """
        Read an alignment file and return the alignement as a Dictionnary
        - Keys = uniprot_id
        - Values = alignement (FASTA format)
        Args:
            alignementfile (string): path to the alignment file
        Returns:
            seqDict (dict): Dictionnary with all aligned sequences.
        """

        from Bio import AlignIO

        # Transformation in a list if only 1 unique file... it will be easier to deal with
        seqDictFull = dict()
        # entropyDictFull = dict()
        if type(alignmentfiles) == type(''):
            alignmentfiles = [alignmentfiles]

        for alignmentfile in alignmentfiles:
            alignment = AlignIO.read(alignmentfile, "fasta")

            # Calculate Shannon entropy
            msa = [str(x.seq).replace('.', '-').upper() for x in alignment]
            # aliSize = len(msa[0])
            # entropy = []
            # for i in range(aliSize):
            #    entropy.append(calculate_shannon_entropy([x[i] for x in msa]))

            seqDict = defaultdict(dict)

            regexFull = re.compile("^(\S+_\S+)\|(\S+)\/(\d+-\d+): (.+)\|(\w+)(\/.+)?")
            #   regex =re.compile("^(\S+_\S+)\|(\S+)\/(\d+-\d+):")

            prositeDomName = None
            for record in alignment:
                header = record.description
                match = regexFull.match(header)
                if match:
                    uniprot_id = match.group(1)
                    uniprot_acc = match.group(2)
                    seqrange = match.group(3)
                    prositeName = match.group(4)
                    prositeID = match.group(5)
                    # unkownparameter = match.group(6)
                else:
                    print(header)
                    uniprot_id = "error"
                    seqrange = "error"
                    uniprot_acc = "error"
                    domainName = "error"
                # uniprot_id = record.id.split("/")[0]
                # seqrange = record.id.split("/")[1]
                # Some uniprot_id sequences are splitted in two, resulting multiple sequences with the same "uniprot_id"

                seqDict[uniprot_id][seqrange] = record
            # entropyDictFull[prositeName + "|" + prositeID] = entropy
            seqDictFull[prositeName + "|" + prositeID] = seqDict
        return seqDictFull

    def match_residue_number_with_alignment_position(self, DATASET):
        """

        Args:
            DATASET:

        Returns:

        """
        self.CATHPDB_PROSITESEQ = {}

        def get_seq_from_pdbdf(pdbdf):
            """
            Convert a PDB into a sequence
            Just take the Amino-acid 3-letter code and convert it into a 1letter code.
            Return the sequence as a string
            Args:
                pdbfg (pd.DataFrame): the pdb as a DataFrame
            Returns:
                seq (string): sequence of the PDB.
            """
            # First : convert the sequences.
            letters = {
                "ALA": "A",
                "ARG": "R",
                "ASN": "N",
                "ASP": "D",
                "CYS": "C",
                "GLU": "E",
                "GLN": "Q",
                "GLY": "G",
                "HIS": "H",
                "ILE": "I",
                "LEU": "L",
                "LYS": "K",
                "MET": "M",
                "PHE": "F",
                "PRO": "P",
                "SER": "S",
                "THR": "T",
                "TRP": "W",
                "TYR": "Y",
                "VAL": "V",
            }

            CA = pdbdf.query("atom_name == 'CA'")

            resnames = CA.residue_name.values
            resid = CA.residue_number.values
            seq3 = []
            seq3.append(resnames[0])
            for i in range(1, len(resid)):
                if resid[i] != resid[i - 1]:
                    seq3.append(resnames[i])

            seq = "".join([letters[x] for x in seq3])
            return seq



        def return_diff_with_score(sequence, alignment):
            """
            Compute the differences between an alignement (in FASTA) and a sequence.
            The calculation use the "difflib package" (work just like the diff command)
            Also Calcul a very basic score between the alignment and the Sequence.
            +1 with a match, 0 otherwise...
            Args:
                sequence (str): Sequence of the PDB
                alignment (str): alignement of the corresponding uniprotID (in prosite)
            Returns:
                diff (list): list of differences in every sequence
                score (int): matching score.

            """
            # matcher = difflib.SequenceMatcher(a=sequence, b=alignment)
            # match = matcher.find_longest_match(0, len(matcher.a), 0, len(matcher.b))

            d = difflib.Differ()
            diff = list(d.compare(alignment, sequence))
            # diff = list(d.compare(sequence, alignment))

            score = 0
            for element in diff:
                if element[0] == " ":
                    score += 1
            return (diff, score)

        def replace_in_string(string, position, letter):
            stringList = list(string)
            stringList[position] = letter
            return "".join(stringList)

        def is_mutation(i, diff):
            nbMinus = 0
            nbPlus = 0
            for j in range(i + 1, len(diff) - 1):
                if nbMinus != 0 and nbMinus == nbPlus:
                    return True

                if diff[j][0] == "-":
                    if diff[j][2] not in ["-", "."]:
                        nbMinus += 1
                    else:
                        return False

                elif diff[j][0] == "+":
                    if diff[j][2] not in ["-", "."]:
                        nbPlus += 1
                    else:
                        return False
                else:
                    return False

        def find_position_in_sequence(sequence: str, recordDict: dict):
            """
            This COMPLEEEEEXE function is made to match the alignment position and the sequence.
            To be honest, I don't even remember how it ended up like this.... It works that's all T_T
            Args:
                sequence (str): Sequence object
                recordDict (dict): TODO

            Returns:

            """
            if len(recordDict) <= 1:
                alignment = str(list(recordDict.values())[0].seq).upper()
                diff, score = return_diff_with_score(sequence, alignment)
            else:  # if there the uniprotID has multiples entry in the alignment
                results = {}
                for r in recordDict.values():
                    alignment = str(r.seq).upper()
                    diff, score = return_diff_with_score(sequence, alignment)
                    results[score] = (diff, alignment)
                best = max(results.keys())
                diff = results[best][0]
                alignment = results[best][1]
                score = best

            # Since a uniprot sequence can contains several of the same domain,
            # And only on is referenced in prosite, you can have mismatches
            # So here we calculate the "alignment cover percentage" and if it is bellow
            # than 95%, we remove it !
            seqInAliSize = len(alignment.replace(
                " ", "").replace(".", "").replace("-", ""))
            if score / seqInAliSize < 0.75:
                return None, None

            position = []
            # first search first match. it will be index 0!
            indexbase = 0
            StartIndexInDiff = 0
            for i in range(len(diff)):
                if diff[i][0] == "-":
                    if diff[i + 1][0] in ["-"]:
                        indexbase += 0  # '-' at start means that the alignment is bigger
                if diff[i][0] == "+":
                    indexbase -= 1
                if diff[i][0] == " ":
                    # let's say that we need at least 3 residue to match to avoid false positiv like SSQH with SQH
                    if diff[i + 1][0] == " ":  # and diff[i+3][0] == ' ':
                        # indexbase = -i +minIncr
                        StartIndexInDiff = i
                        break

            mutationFound = False
            nbMutation = 0
            position = []
            missing = 0
            numberOfGoodMatch = 0
            insertion = 0
            deletion = 0
            add = False
            retain = 0
            startingGaps = 0

            for i in range(len(diff)):
                match = diff[i]
                if match[0] == " ":
                    add = True

                elif match[0] == "+":
                    if i > abs(StartIndexInDiff):  # look after start
                        # if i > abs(indexbase)+retain:
                        if i < len(diff) - 1:
                            if diff[i - 1][0] == " " and diff[i + 1][0] == " ":
                                retain += 1
                                position.append(np.NaN)
                                continue
                            elif diff[i - 1][0] == "-" and diff[i - 1][2] not in [".", "-"]:
                                retain += 1
                                continue
                            elif nbMutation != 0:
                                nbMutation -= 1
                                retain += 1
                                continue

                        insertion += 0
                    add = True

                elif match[0] == "-":
                    if i < len(diff) - 1:
                        if i > StartIndexInDiff:
                            if diff[i + 1][0] == "+":
                                add = True
                            elif is_mutation(i, diff):
                                add = True
                                nbMutation += 1
                    if i < abs(StartIndexInDiff) and diff[i][2] == "-":
                        startingGaps += 1
                    if match[2] == "." and diff[i - 1][0] == "+":
                        retain += 1

                if add == True:
                    pos = indexbase + i + insertion - deletion - retain  # - startingGaps

                    position.append(pos)
                    add = False

                elif add == "NaN":
                    position.append(np.NaN)
                    add = False
            return alignment, position

        #######

        def add_position_in_PDB(pdbdf, alignmentDicts):
            # First we need to find the write alignement dict
            residueIndex = pdbdf.groupby(["residue_number", "chain_id"]).ngroup()
            pdbdf["residue_index"] = residueIndex
            # get sequence
            sequence = get_seq_from_pdbdf(pdbdf)  # 1.5ms
            # get alignment
            uniprot_id = pdbdf.uniprot_id.iloc[0]  # take first uniprot_ID as reference
            # Biopython record alignment (as a dict, the key is the range)

            # Search for our uniprotID in the a alignent File.
            uniprotInAlignment = False
            for name, dico in alignmentDicts.items():
                if uniprot_id in dico.keys():
                    prositeName = name.split("|")[0]
                    prositeID = name.split("|")[1]
                    alignmentDict = dico
                    uniprotInAlignment = True
                    #entropy = entropyDict[name]
                    break

            if uniprotInAlignment == False:
                # print(f"no alignement found for {pdbdf.cathpdb.iloc[0]}")
                pdbdf["prositeName"] = np.NaN
                pdbdf["prositeID"] = np.NaN
                pdbdf["alignment_position"] = np.NaN
                return pdbdf

            recordDict = alignmentDict[uniprot_id]
            # alignement = str(record.seq).upper()
            # get position INDEX based on the alignement
            alignment, positions = find_position_in_sequence(sequence, recordDict)  # <1ms
            # CATHPDB_PROSITESEQ[]
            # if pdbdf.cathpdb.iloc[0] in CATHPDB_PROSITESEQ:
            #    # print(f"DUPLICATE {pdbdf.cathpdb.iloc[0]}")
            #    1 + 1
            self.CATHPDB_PROSITESEQ[pdbdf.cathpdb.iloc[0]] = alignment

            if positions == None:  # if there is any problem to find the position (see comments in find_position_in_sequence function)
                pdbdf["alignment_position"] = np.NaN
                return pdbdf  # just return pdbdf.
            for x in pdbdf.residue_index:
                try:
                    positions[x]
                except:
                    print(x)
            # this list comprehension is made to extend the PER AA numbering to PER ATOM.
            pdbdf["alignment_position"] = [positions[x] for x in pdbdf.residue_index]
            pdbdf["prositeName"] = prositeName
            pdbdf["prositeID"] = prositeID
            # add entropy
            #entropyPosition = defaultdict(lambda: np.NaN)
            #for i in range(len(entropy)):
            #    entropyPosition[i] = entropy[i]
            #pdbdf["entropy"] = [entropyPosition[x] for x in pdbdf.alignment_position]
            return pdbdf

        def check_seq_ali(cathpdb, DATASET):
            CA = DATASET.query("cathpdb == @cathpdb and atom_name == 'CA'")
            pos = CA.alignment_position
            resname = CA.residue_name
            resPosInAli = dict(zip(pos, resname))

            alignment = self.CATHPDB_PROSITESEQ[cathpdb]
            alisize = len(alignment)

            match = 0
            total = 0
            unmatch = []
            for p in pos:
                try:
                    p = int(p)
                except:
                    continue  # Value is NaN which means that it is an insertion
                if p >= 0 and p < alisize:
                    total += 1
                    # resname = CA.query("alignment_position == @p").residue_name.values[0]
                    resname = resPosInAli[p]
                    resname1letter = self.map3to1[resname]

                    if resname1letter == alignment[p]:
                        match += 1
                    else:
                        if alignment[p] in ["-", "."]:
                            match += 1
                        else:
                            unmatch.append(f"{resname1letter}-{alignment[p]}={p}")
            matching = match / total * 100
            return (matching, unmatch)

        def iterate_on_domain_add_sequence_index_to_pdb(domdf, dom):
            """
            This is the master loop function to iterate over domains in the dataset (with a groupby("domain"))
            This fetch the fasta alignement file and the iterate over cathpdb to update the sequence index.
            """


            paths = self.DOMAIN_SEQ[dom]
            isPathsOK = True
            # check if every paths are okay
            if type(paths) == type(''):
                paths = [paths]

            for path in paths:
                if not os.path.exists(path):
                    isPathsOK = False

            if isPathsOK == True:
                # read the sequence as a dict (the key is uniprotID)
                alignmentDicts = self.return_aligned_seq_as_dict(paths)
            else:
                print("no alignement file found. please check DOMAIN_SEQ object")
                return domdf

            # update the sequence position index in alignment for PDB residue ID.
            domdf = domdf.groupby("cathpdb").progress_apply(add_position_in_PDB, alignmentDicts=alignmentDicts)

            return domdf

        # Work on a copy for now.... just in case.
        # DATASET_SEQ = DATASET.copy(deep=True)
        print("> Matching between the structure and the sequence aligned...")
        DATASET = DATASET.groupby("domain").progress_apply(
            lambda x: iterate_on_domain_add_sequence_index_to_pdb(x, x.name))

        return (DATASET)




    def add_sequence_in_dataset(self, DATASET):
        """
        This is the master loop function to iterate over domains in the dataset (with a groupby("domain"))
        This fetch the fasta alignement file and the iterate over cathpdb to update the sequence index.
        """
        regex = re.compile("^(\S+_\S+)\|(\S+)\/(\d+-\d+): (.+)\|(.+)")

        SequencesDataframeList = []
        for dom in self.tqdm(self.DOMAIN_SEQ):
            print(dom)
            paths = self.DOMAIN_SEQ[dom]

            isPathsOK = True
            # check if every paths are okay
            if type(paths) == type(''):
                paths = [paths]

            for path in paths:
                if not os.path.exists(path):
                    print("no alignement file found. please check DOMAIN_SEQ object")
                    return domdf

            alignmentDicts = self.return_aligned_seq_as_dict(paths)

            # Instantiate the list that will contain all the sequences
            seqDatasetList = []
            # read the sequence as a dict (the key is uniprotID)

            # replace PI3K_C2 by C2 (prosite has 2 C2 domain, I will just keep PI3K_C2)
            # if dom == "PI3K_C2":
            #    dom = "C2"

            for name, dico in alignmentDicts.items():
                prositeName = name.split("|")[0]
                prositeID = name.split("|")[1]
                alignmentDict = dico
                #entropy = entropyDicts[name]

                #entropyPosition = defaultdict(lambda: np.NaN)
                #for i in range(len(entropy)):
                #    entropyPosition[i] = entropy[i]

                # for each sequence
                for uniprot_id in alignmentDict:

                    records = alignmentDict[uniprot_id]

                    for alirange, ali in records.items():
                        header = ali.description
                        match = regex.match(header)
                        uniprot_acc = match.group(2)
                        # this (complex) comprehension list is just to get the index number of each amino acids while ignoring gaps.
                        position_in_alignement = [i for i in range(len(ali.seq)) if ali.seq[i] not in [".", "-"]]
                        # keep only the sequence
                        seq = [x.upper() for x in ali.seq if x not in [".", "-"]]
                        start = int(match.group(3).split("-")[0])
                        residue_numbers = [i + start for i in range(len(seq))]
                        # NOTE : Secondary structure prediction is possible with raptorX-SS8. No jpred since Jpred is only on webserver.
                        N = len(seq)
                        data = {
                            "atom_name": np.repeat("CA", N),
                            "ali_range": np.repeat(alirange, N),
                            "residue_name": [self.map1to3[x] for x in seq],
                            "domain": np.repeat(dom, N),
                            "uniprot_id": np.repeat(uniprot_id, N),
                            "uniprot_acc": np.repeat(uniprot_acc, N),
                            "data_type": np.repeat("prosite", N),
                            "residue_index": position_in_alignement,
                            "residue_number": residue_numbers,
                            "alignment_position": position_in_alignement,
                            "prositeName": np.repeat(prositeName, N),
                            "prositeID": np.repeat(prositeID, N),
                            #"entropy": [entropyPosition[x] for x in position_in_alignement],
                        }
                        seqDatasetList.append(pd.DataFrame.from_dict(data))


            seqDatasetList = pd.concat(seqDatasetList)
            seqDatasetList.reset_index(drop=True)
            SequencesDataframeList.append(seqDatasetList)
        SequencesDataframe = pd.concat(SequencesDataframeList)
        SequencesDataframe.reset_index(drop=True)


        print("> Adding the Sequence data in the dataset...")
        DATASET = pd.concat([DATASET, SequencesDataframe]).reset_index(drop=True)

        # Adding the ORIGIN
        DATASET["origin"] = DATASET["uniprot_id"].str.split("_").str[1]
        return DATASET


    def download_uniprot_data(self, DATASET):
        """
        Download UNIPROT page for every entries we have.
        On uniprot xml page you can have a lot of information that we can add on our dataset later on.

        Returns:
            None, it's just to download datas
        """
        #Get all differents uniprot ids
        uniprot_ids = DATASET.uniprot_id.dropna().unique()
        import requests

        link = "https://www.uniprot.org/uniprot/"

        exists = 0
        downloaded = 0

        for uniprot_id in self.tqdm(uniprot_ids):
            file = self.UNIPROTFOLDER + uniprot_id + ".xml"
            if os.path.isfile(self.UNIPROTFOLDER + uniprot_id + ".xml"):
                exists += 1
                continue
            url = link + uniprot_id + ".xml"
            r = requests.get(url)
            open(file, "wb").write(r.content)
            downloaded += 1

        # note: it takes 30mn for 4200 sequences...

        print(f"> {downloaded} new files downloaded and {exists} will be reused.")



    def add_info_from_uniprot(self, DATASET):
        from bs4 import BeautifulSoup

        def read_xml_file(path="/Users/thibault/Documents/WORK/dataset_tests/Q08945.xml"):
            """
            Read a XML file and return the XML object
            Args:
                infile (str): path to the xml file
            returns:
                soupe (beautifulSoup object): xml object.
            """
            infile = open(path, "r")
            contents = infile.read()
            soup = BeautifulSoup(contents, "xml")
            infile.close()
            return soup

        def parse_xml(soup):
            """
            Parse XML file
            """
            try:
                uniprot_id = soup.uniprot.entry.find("name").text
            except:
                uniprot_id = np.NaN

            try:
                uniprot_acc = soup.find_all("accession")[
                    0
                ].text  # take only the citable uniprot_acc
            except:
                uniprot_acc = np.NaN

            try:
                location = [
                    x.text
                    for x in soup.find_all("comment", {"type": "subcellular location"})[
                        0
                    ].find_all("location")
                ]
            except:
                location = ["unkown"]

            taxonsXML = soup.find_all("taxon")
            taxon = '/'.join([x.text for x in taxonsXML[:2]])

            try:
                prosite = [x["id"] for x in soup.find_all("dbReference", {"type": "PROSITE"})]
                prosite_name = [
                    x.find_all("property", {"type": "entry name"})[0]["value"]
                    for x in soup.find_all("dbReference", {"type": "PROSITE"})
                ]
            except:
                prosite = ["unknown"]
                prosite_name = ["unknown"]

            dataDict = {
                "uniprot_id": [uniprot_id],
                "uniprot_acc": [uniprot_acc],
                "location": [location],
                "CR:prositeID": [prosite],
                "taxon": [taxon],
                "CR:prositeName": [prosite_name],
            }
            df = pd.DataFrame.from_dict(dataDict)
            # df = df.explode("location").explode("CR:prosite").reset_index(drop=True)

            return df

        def create_DF_XML(repo="orinoco"):
            if repo == "orinoco":
                XMLfolder = "/net/orinoco/thibault/projects/peprmint/database/uniprot"
            elif repo == "local":
                XMLfolder = self.UNIPROTFOLDER
            else:
                print("unknown folder type. try again.")
                return None

            XMLfiles = glob.glob(XMLfolder + "/*.xml")
            XMLlist = []
            for xmlPath in self.tqdm(XMLfiles):
                soup = read_xml_file(xmlPath)
                df = parse_xml(soup)
                XMLlist.append(df)

            parsedXML = pd.concat(XMLlist).reset_index(drop=True)
            return parsedXML


        print("> Parsing uniprot files")
        UNIPROTDATA = create_DF_XML("local")
        DATASET = DATASET.merge(how="left", on="uniprot_id", right=UNIPROTDATA.drop("uniprot_acc", axis=1))
        return DATASET


    def add_conservation(self, DATASET, uniref='uniref90', normalized=True, gapcutoff=0.8):
        """
        Add conservation and entropy calculation

        Args:
            DATASET:

        Returns:

        """

        def norm_0_1(array, minvalue=0, maxvalue=None):
            array = np.array(array)

            gaps_indexes = np.where(array == -1)
            if not maxvalue:
                maxvalue = array.max()

            array = 1 - ((array - minvalue) / (maxvalue - minvalue))
            #array[array > 1] = 1
            array[gaps_indexes] = -1

            return array

        def get_representative_uniprot_accs(unirefCluster):
            """
            Return only representative uniprot_acc.
            Args:
                unirefCluster (dict): dictionnary of acc:cluster_id
            """

            df = pd.DataFrame([list(unirefCluster.keys()),
                               list(unirefCluster.values())],
                              index=["uniprot_acc", "cluster"]).T

            def select_seq(cluster):
                clusterid = cluster.cluster.unique()[0]
                accs = cluster.uniprot_acc.dropna().unique()
                accs.sort()
                if clusterid in accs:
                    return clusterid
                else:
                    return accs[0]

            results = df.groupby("cluster").apply(select_seq)

            return results.unique()

        def filter_alignment(alignment, repr_accs):
            """
            return only alignments which uniprot_acc is in "repr_accs"
            """
            from Bio.Align import MultipleSeqAlignment
            uniprot_acc_RE = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
            keep_ali = []
            for seq in alignment:
                seqid = seq.id
                search = uniprot_acc_RE.search(seqid)
                if search:
                    uniprot_acc = search.group(0)
                    if uniprot_acc in repr_accs:
                        keep_ali.append(seq)

            return MultipleSeqAlignment(keep_ali)

        def convert_H10_as_array(alignment):
            """
            Convert all sequences in alignment into a H10 index following the convention
            of Mirny and Shakhnovich (1999) (SEE http://thegrantlab.org/bio3d/reference/entropy.html)

            Return the alignment as
            """
            # Define Conversion Dict following the convention of Mirny and Shakhnovich (1999) (SEE http://thegrantlab.org/bio3d/reference/entropy.html)
            H10conversion = {'V': 'B',  # HYDROPHOBIC / ALIPHATIC
                             'I': 'B',  # HYDROPHOBIC / ALIPHATIC
                             'L': 'B',  # HYDROPHOBIC / ALIPHATIC
                             'M': 'B',  # HYDROPHOBIC / ALIPHATIC
                             'F': 'A',  # AROMATIC
                             'W': 'A',  # AROMATIC
                             'Y': 'A',  # AROMATIC
                             'S': 'S',  # SER/THR
                             'T': 'S',  # SER/THR
                             'N': 'O',  # POLAR
                             'Q': 'O',  # POLAR
                             'H': 'H',  # POSITIVE
                             'K': 'P',  # POSITIVE
                             'R': 'P',  # POSITIVE
                             'D': 'N',  # NEGATIVE
                             'E': 'N',  # NEGATIVE
                             'A': 'T',  # TINY
                             'G': 'T',  # TINY
                             'P': 'P',  # PROLINE
                             'C': 'C',  # CYSTEINE
                             '-': '-',  # GAPS OR NONSTANDARD AMINO ACIDS
                             'X': '-',  # GAPS OR NONSTANDARD AMINO ACIDS
                             '.': '-',  # GAPS OR NONSTANDARD AMINO ACIDS
                             }
            alignment_array = []
            for seq in alignment:
                alignment_array.append(np.array([H10conversion[x.upper()] for x in seq.seq]))
            return np.array(alignment_array)

        def calculate_shannon_entropy(column, gapcutoff=0.1):
            """
            Calculate the shannon entropy on a column of a multiple alignment sequence
            based on https://github.com/ffrancis/Multiple-sequence-alignment-Shannon-s-entropy/blob/master/msa_shannon_entropy012915.py
            Args:
                column <list>: all the amino acids aligned in a column
            Returns:
                sh_entropy <float>: shannon entropy of this specific alignment position
            """
            import math
            unique_aa = set(column)  # get all the unique amino acids in a column
            # gaps = ['-','.','X']
            # unique_aa_nogap = [x for x in unique_aa if x not in gaps]

            gap_counts = column.count('-')
            gap_counts += column.count('.')
            gap_counts += column.count('X')
            if (gap_counts / len(column)) > gapcutoff:
                return -1

            M = len(column)
            # M = len(unique_aa_nogap)
            entropy_list = []
            for amino in unique_aa:
                if not amino == '-':  # do not calculate H for gaps, but keep it as diversity
                    n_i = column.count(amino)  # How many times we have this amino acids
                    P_i = n_i / M  # The propability to get it
                    entropy_i = P_i * (math.log(P_i, 2))  # multiply by a log2
                    entropy_list.append(entropy_i)
            sh_entropy = -(sum(entropy_list))
            return sh_entropy

        def calc_conservation_per_domain(alignmentfiles, repr_accs, normalize = True,gapcutoff=0.8):
            """
            Read an alignment file and return the alignement as a Dictionnary
            - Keys = uniprot_id
            - Values = alignement (FASTA format)
            Args:
                alignementfile (string): path to the alignment file
                repr_accs (list(str)): list of representative structures.
            Returns:
                seqDict (dict): Dictionnary with all aligned sequences.
            """

            from Bio import AlignIO
            regexFull = re.compile("^(\S+_\S+)\|(\S+)\/(\d+-\d+): (.+)\|(\w+)(\/.+)?")

            entropyDictFull = defaultdict(lambda: np.NaN)
            entropyH10DictFull = defaultdict(lambda: np.NaN)
            if type(alignmentfiles) == type(''):
                alignmentfiles = [alignmentfiles]

            for alignmentfile in alignmentfiles:
                alignment = AlignIO.read(alignmentfile, "fasta")
                alignment_filtered = filter_alignment(alignment, repr_accs)

                # Calculate Shannon entropy
                msa = [str(x.seq).replace('.', '-').upper() for x in alignment_filtered]
                msaH10 = convert_H10_as_array(alignment_filtered)
                aliSize = len(msa[0])
                entropy = []
                entropyH10 = []



                for record in alignment:
                    header = record.description
                    match = regexFull.match(header)
                    if match:
                        prositeName = match.group(4)
                        prositeID = match.group(5)
                    else:
                        print(header)

                shannon = []
                shannonH10 = []
                for i in range(aliSize):
                    shannon.append(calculate_shannon_entropy([x[i] for x in msa], gapcutoff))
                    shannonH10.append(calculate_shannon_entropy([x[i] for x in msaH10],gapcutoff))

            #We will use a dataframe to store our data because it's faster to merge with our original dataset later on
            if normalize:
                shannon = norm_0_1(shannon, maxvalue=4.34)
                shannonH10 = norm_0_1(shannonH10, maxvalue=3.15)

            df = pd.DataFrame(list(zip(shannon, shannonH10)),
                              columns=["shannon", 'shannonH10'])
            df["prositeName"] = prositeName
            df["prositeID"] = prositeID
            df = df.reset_index()
            #transform the index (0-based) into "alignment_position" which is also 0-based
            df = df.rename(columns={'index':'alignment_position'})

            return df


        ############
        # Main function
        ############
        conservationList = []
        for dom in self.tqdm(DATASET.domain.unique()):
            domain = DATASET.query("domain == @dom")

            if uniref not in ['uniref50','uniref90','uniref100']:
                print("wrong Uniref given, choosing uniref90")
                uniref = 'uniref90'

            unirefCluster = {k: v for k, v in [tuple(x) for x in domain[["uniprot_acc", uniref]].drop_duplicates().dropna().values]}
            alignment_files = self.DOMAIN_SEQ[dom]
            repr_accs = get_representative_uniprot_accs(unirefCluster)
            conservationList.append(calc_conservation_per_domain(alignment_files, repr_accs, normalized, gapcutoff))


        #Concat all dataframe list into an original dataset
        conservationdf = pd.concat(conservationList)
        conservationdf = conservationdf.reset_index(drop=True) #Just in case.....

        #Then we merge based on "prositeID, prositeName and alignment_position
        #If no match then shannon&shannonH10 == NaN
        return  pd.merge(DATASET,conservationdf, on=["prositeName","prositeID","alignment_position"],how="outer")













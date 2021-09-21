#!/usr/bin/env python3

import os
import glob
import subprocess
from datetime import datetime
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
import tempfile
import textwrap
import pandas as pd
import numpy as np
from Bio import SeqIO
import json


def load_map(map_path):
    """Load map in json format and return a dic

    Args:
        map_path (str): Path to json file

    Returns:
        map_dict: Map as dictionary.
    """
    with open(map_path) as map_handle:
        map_dict = json.load(map_handle)
    return map_dict


def map_protein_to_taxon(fasta_path):
    """Maps a FASTA file, parses the taxon ID (OX=#) and stores the key, value
    in a dict.

    Args:
        fasta_path ([type]): [description]
    """
    with open(fasta_path) as handle:
        fasta = SeqIO.read(handle, format="fasta")
        uniprot_id = fasta.id.split("|")[1]
        description_fields = fasta.description.split()
        ox = [x for x in description_fields if x.startswith("OX=")]
        if len(ox) == 0:
            print(f"OX not found in header for {uniprot_id}")
            taxon = ""
        else:
            taxon = ox[0].split("=")[1]
        return uniprot_id, taxon


def run_iupred(
    input_folder,
    output_folder,
    iupred_flag,
    regex="*.fasta",
    iupred_path="./software/iupred/",
):
    """Runs IUPred in the specified mode ("short" or "long" flags)
    over multiple, multirecord FASTA files. Sequences are processed
    in temporary files and the output .iup file contains all concatenated
    calculations. Input/output folders and, optionally, a regular expression
    are required.

    Arguments:
        input_folder {str} -- [description]
        output_folder {str} -- [description]
        iupred_flag {str} -- [description]
        iupred_flag {str} -- [description]

    Keyword Arguments:
        regex {str} -- [description] (default: {'*.fasta'})

    Returns:
        {file } --  A single .iup file formatted as following (no header):

    index    residue_number    residue_type    iupred_score    record_id
    """
    # Set IUPRed environment variables for interaction
    my_env = os.environ.copy()
    my_env["IUPred_PATH"] = iupred_path

    input_files = glob.glob(input_folder + regex)
    abs_run_count = 0
    # Start log file
    log_handle = open(output_folder + "iupred-" + iupred_flag + ".log", mode="w+")
    log_handle.write(datetime.now().strftime("%y/%m/%d %H:%M:%S") + "\n")
    log_handle.write(
        "Disorder prediction performed with IUPred" + "\n" + 79 * "-" + "\n"
    )
    log_handle.write("protein_ac\tdataset")
    log_handle.write("\n")

    for file in input_files:
        dataset = file.split("/")[-1].split(".")[0]
        resfile = output_folder + dataset + "_iupred-" + iupred_flag + ".table"
        # Create temporary results file
        tempres = tempfile.NamedTemporaryFile(mode="a+")
        # Create an index to save from which record comes the prediction.
        indexlist = []
        fasta_handle = SeqIO.parse(file, format="fasta")
        # Due to different fasta description formatting schemes,
        # here comes a try block.
        # The most common "|"-separated description is the default parsing
        runcount = 0
        print(79 * "-" + "\n", dataset, "\n" + 79 * "-")
        for record in fasta_handle:
            try:
                protein_ac = record.id.split("|")[1]
            except IndexError:
                protein_ac = record.id
            # Adapt sequence to FASTA format column number specification.
            seq = textwrap.fill(str(record.seq), width=60)
            # Create index entries for all aa
            indexlist.extend([protein_ac] * len(record.seq))
            tempseq = tempfile.NamedTemporaryFile(mode="a+")
            print(
                """Running IUPred ({} mode) \
on record {} in tempfile {}""".format(
                    iupred_flag, protein_ac, tempseq.name
                )
            )
            tempseq.write(">{}\n{}".format(protein_ac, seq))
            # This is very important, makes it go to the beginning of the file!
            tempseq.seek(0)
            # If stdout is set to None, it will stream to the terminal.
            # Useful  for testing purposes
            iupredproc = subprocess.Popen(
                [iupred_path + "iupred", tempseq.name, iupred_flag],
                env=my_env,
                stdout=tempres,
            )
            iupredproc.communicate()
            runcount += 1
            tempseq.close()
            log_handle.write("{}\t{}".format(protein_ac, dataset))
            log_handle.write("\n")
        # Load temporary results file as a pandas DataFrame. Easy parsing.
        abs_run_count += runcount
        df = pd.read_csv(
            filepath_or_buffer=tempres.name,
            header=None,
            sep=" ",
            skipinitialspace=True,
            comment="#",
            names=["position", "residue", "score"],
        )
        indexseries = pd.Series(data=[x for x in indexlist], name="protein_ac")
        df["protein_ac"] = indexseries
        print("Saving IUPred scores as {} ...".format(resfile))
        df.to_csv(path_or_buf=resfile, sep="\t", header=False)
        print(
            """{} IUPred runs performed \
for records in file {}""".format(
                runcount, file
            )
        )
    log_msg = """{} IUPred runs from protein records \
in {} fasta files""".format(
        abs_run_count, len(input_files)
    )
    log_handle.write(79 * "-" + "\n" + log_msg)


def load_iupred_scores(input_folder, regex="*.table"):
    """Parses all .table files containing concatenated IUPred scores in a folder
    produced by the run_iupred function and returns an iterator.
    Additionally, two extra columns are added specifying a) from which
    fasta record collection the scores come from, and b) the IUPred
    prediction mode used during the run.
    Getting a generator avoids loading in memory all scores at the same time.

    Arguments:[type]
        input_folder {string} -- Path to IUPred results folder.

    Keyword Arguments:
        regex {str} -- Optional regular expression to select subsets of results
                       (default: {'*.table'})

    Yields:
        [string, pandas.DataFrame] -- Name of the dataset (e.g. a proteome
        accession number), DataFrame of concatenated disorder predictions.
    """

    filenames = glob.glob(input_folder + regex)
    for f in filenames:
        filename = f.split("/")[-1].split(".")[0].split("_")
        dataset = filename[0]
        disorder_mode = filename[-1]
        # Parse the results file with appropriate header and data types to avoid
        # problems later.
        score = pd.read_csv(
            filepath_or_buffer=f,
            sep="\t",
            names=["index", "residue_number", "amino_acid", "score", "protein_ac"],
            dtype={
                "index": np.dtype(int),
                "residue_number": np.dtype(int),
                "amino_acid": np.dtype(str),
                "score": np.dtype(float),
                "protein_ac": np.dtype(str),
            },
            index_col="index",
        )
        # Add collection and disorder prediction data in new columns
        score["dataset"] = dataset
        score["disorder_predictor"] = disorder_mode
        yield filename, score


def load_disopred_scores(input_folder, regex="*diso"):
    """Parses all .diso files containing protein per-aa scores in a folder and
    returns an iterator. Additionally, an extra column is added specifying from
    which collection the scores come from. Getting a generator avoids loading
    in memory all scores at the same time.

    Arguments:[type] input_folder {string} -- Path to Disopred results folder.

    Keyword Arguments: regex {str} -- Optional regular expression to select
        subsets of results (default: {'*.diso'})

    Yields: [string, pandas.DataFrame] -- Name of the dataset (e.g. a proteome
        accession number), DataFrame of concatenated disorder predictions.
    """

    filenames = glob.glob(input_folder + regex)
    for f in filenames:
        # This works when files are named after the uniprot ID
        protein_ac = f.split("/")[-1].split(".")[0]
        # This works when the string is not properly formatted
        # e.g. sp_UNIPROTID_GENENAME
        if len(protein_ac.split("_")) > 1:
            protein_ac = protein_ac.split("_")[1]
        dataset = f.split("/")[-2]
        disorder_predictor = "disopred31"
        # Parse the results file with appropriate header and data types to avoid
        # problems later.
        score = pd.read_csv(
            filepath_or_buffer=f,
            skiprows=3,
            delim_whitespace=True,
            skipinitialspace=True,
            names=["residue_number", "amino_acid", "symbol", "score"],
            dtype={
                "residue_number": np.dtype(int),
                "amino_acid": np.dtype(str),
                "symbol": np.dtype(str),
                "score": np.dtype(float),
            },
        )
        score.drop(columns=["symbol"], inplace=True)
        # Add collection and disorder prediction data in new columns
        score["protein_ac"] = protein_ac
        score["dataset"] = dataset
        score["disorder_predictor"] = disorder_predictor
        yield [protein_ac, disorder_predictor], score


def calc_disfrac(input_df, threshold=0.5):
    """Take a merged IUPRed input (multiple proteins) and calculate
    the fraction of disordered residues according to a specific threshold for
    each case. A pd.GroupBy object is first generated and then the function
    applied.

    Arguments:
        input_df {string} -- [description]

    Keyword Arguments:
        threshold {float} -- [description] (default: {0.5})

    Returns:
        [type] -- [description]
    """
    grouped = input_df.groupby(["protein_ac", "dataset", "disorder_predictor"])
    # Aggregate IUPred scores using the disorder fraction formula
    disfrac = grouped.agg({"score": lambda x: sum(x >= threshold) / len(x)})
    # Rename score column to disorder fraction
    disfrac = disfrac.rename(columns={"score": "disorder_fraction"})
    disfrac.reset_index(inplace=True)

    return disfrac


def categorize_disorder(scores_df, disfrac_df):
    """Take a merged disorder input and the derived aggregated disorder fractions
    and classify proteins categorically according to Deiana et al.
    (http://dx.doi.org/10.1101/446351)

    The criteria used are:
    A10 = 10% disorder fraction
    A30 = 30% disorder fraction
    B = 22 consecutive disordered amino acids (we need the raw scores!)

    The disordered classes are defined below.

    Arguments:
        scores_df {pd.DataFrame} -- [description]
        disfrac_df {pd.DataFrame} -- [description]

    Returns:
        cat_disfrac {pd.DataFrame} -- [description]
    """
    # First, conditions are checked for each protein

    condA10 = {}
    condA30 = {}
    condB = {}

    grouped_disfrac = disfrac_df.groupby(["protein_ac"])
    grouped_scores = scores_df.groupby(["protein_ac"])

    # Condition A: Disorder fraction
    for gp, df in grouped_disfrac:
        if df["disorder_fraction"].unique() < 0.1:
            condA10[gp] = False
            condA30[gp] = False
        elif df["disorder_fraction"].unique() < 0.3:
            condA10[gp] = True
            condA30[gp] = False
        elif df["disorder_fraction"].unique() >= 0.3:
            condA10[gp] = True
            condA30[gp] = True
        else:
            print("Disorder fraction value not found for {}".format(gp))

    # Condition B:Disordered stretches
    for gp, df in grouped_scores:
        # Define a rolling window of 22 amino acids
        wdw = df["score"].rolling(22)
        # Check how many stretches exist where all scores indicate disorder
        idrs = df[wdw.min() >= 0.5]
        if len(idrs) > 0:
            condB[gp] = True
        else:
            condB[gp] = False

    disorder_type_dict = {}

    # Categorize disorder fractions
    for protein in disfrac_df["protein_ac"]:
        A10 = condA10[protein]
        A30 = condA30[protein]
        B = condB[protein]

        if A10 is False:
            # 'ORD': Ordered protein
            disorder_type_dict[protein] = "ORD"
        elif A10 is True and A30 is False:
            if B is False:
                # 'NDP': Not disordered protein
                disorder_type_dict[protein] = "NDP"
            else:
                # 'PDR': Partially disordered protein
                disorder_type_dict[protein] = "PDR"
        elif A30 is True:
            if B is False:
                # 'FRAG': Fragmentarily disordered protein
                disorder_type_dict[protein] = "FRAG"
            else:
                # 'IDP':Intrinsically disordered protein
                disorder_type_dict[protein] = "IDP"

    # Copy disorder fraction and add classification column
    cat_disfrac = disfrac_df
    cat_disfrac["disorder_category"] = cat_disfrac["protein_ac"].map(disorder_type_dict)

    return cat_disfrac


def concatenate_results(path_list):
    """Concatenates

    Args:
        path_list (list): Paths to result tables

    Returns:
        df: DataFrame containing all concatenated results
    """
    dflist = []
    for path in path_list:
        tmpdf = pd.read_csv(filepath_or_buffer=path, header=0, sep="\t", index_col=0)
        dflist.append(tmpdf)
    df = pd.concat(dflist)

    return df


def load_single_table(path_to_table):
    df = pd.read_csv(path_to_table, header=0, index_col=0, sep="\t")
    df.sort_index(inplace=True)
    return df


def load_multi_tables(path_to_folder, regex="*.table"):
    path_list = glob.glob(path_to_folder + regex)
    dflist = []
    for path in path_list:
        tmpdf = pd.read_csv(path, header=0, index_col=0, sep="\t")
        tmpdf.set_index(tmpdf["protein_ac"], inplace=True)
        dflist.append(tmpdf)
    df = pd.concat(dflist)
    df.sort_index(inplace=True)
    return df


def calc_classfreq(aggregated_df):

    categories = ["IDP", "PDR", "FRAG", "NDP", "ORD"]
    classfreq = aggregated_df.disorder_category.value_counts() / len(aggregated_df)
    classfreq = classfreq.reindex(categories)
    classfreq.fillna(0, inplace=True)
    return classfreq


def calc_mannwithney(merged_df):
    mwu_stat = {}

    for gp, df in merged_df.groupby('effector_type'):
        ref = df[df['collection_type']=='Reference']['disorder_fraction']
        eff = df[df['collection_type']=='Effector']['disorder_fraction']
        reftaxon = df[df['collection_type']=='Reference']['dataset'].unique()[0]
        _, pval = mannwhitneyu(x=ref, y=eff, alternative='less')
        mwu_stat[gp] = (reftaxon, pval)
    mwu_stat_df = pd.DataFrame.from_dict(mwu_stat, orient='index', columns=["Reference taxon", "p-value"])

    return mwu_stat_df

def calc_kolmogorovsmirnov(merged_df):
    ks_stat = {}

    for gp, df in merged_df.groupby('effector_type'):
        ref = df[df['collection_type']=='Reference']['disorder_fraction']
        eff = df[df['collection_type']=='Effector']['disorder_fraction']
        reftaxon = df[df['collection_type']=='Reference']['dataset'].unique()[0]
        _, pval = ks_2samp(ref, eff)
        ks_stat[gp] = (reftaxon, pval)
    ks_stat_df = pd.DataFrame.from_dict(ks_stat, orient='index', columns=["Reference taxon", "p-value"])

    return ks_stat_df



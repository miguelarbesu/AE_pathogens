# Computational analysis of structural disorder in A/E pathogen effectors.

This repository contains the code written for the bioinformatic study in *A pathogen-encoded signaling receptor mediating host-like interactions through intrinsic disorder (temporary title)* (2021).

We use sequence-based predictors to assess and compare the amount of structural disorder
in proteins from intestinal pathogens. Per-amino acid scores were calculated, aggregated
per-protein and used to classify them.

## Description

Effector collections from Enteropathogenic and Enterohemorrhagic *Escherichia coli*
(EPEC and EHEC, respectively) and *Citrobacter rodentium*, the standard small animal A/E
model, were studied and compared to reference proteomes from the corresponding species.
Additionally, the human reference proteome was analyzed to provide a perspective of the
eukaryotic host.

The predictions were done with [DISOPRED
3](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtu744) and [IUPred
1.0](https://doi.org/10.1093/bioinformatics/bti541). Since DISOPRED is more
computationally expensive the calculations were run in a cluster, so **the code here
only runs the IUPred predictions**.

### Organization

Most functions are imported from the `src/aepathbactdis.py` module. Self-descriptive
Python scripts are run *on each fasta folder* for each predictor (DISOPRED and IUpred
1.0 in *short* and  *long* modes).

As a preparatory step, `map_taxa.py` builds json files mapping the corresponding protein
IDs (Uniprot format) to the corresponding taxon numbers and ultimately reference
proteomes.

The processing pipeline then has two stages:

- Prediction of structural disorder of protein sequence collections (`run_iupred.py`).
- Aggregation and classification of the individual proteins (`aggregate_*.py`).

Jupyter notebooks are used for analysis and plotting (see `/notebooks`). The resulting
summary files are saved to the respective aggregated data folders, and plots are saved
to a `/figures` folder (see [Data availability](##Data-availability)). 

## Data availability

The associated data (sequence collections, prediction scores, aggregated values and
classifications) is accessible [here](https://osf.io/3mka9/).

## Installation and dependencies

Assuming you have Python installed, the required libraries can be installed via `pip
install -r requirements.txt`.

IUPred 1.0 is needed to run the predictions. It can be downloaded
[here](https://iupred.elte.hu/Downloads.php) and installed where desired. The path is
optionally provided to `run_iupred.py` as an argument.


## License

This software is licensed under the MIT license. 

## Contributors

- Miguel Arbesú (miguel.arbesu@gmail.com)
- Tiago Cordeiro (tiago.n.cordeiro@gmail.com)
- Andreas Zanzoni (andreas.zanzoni@gmail.com)


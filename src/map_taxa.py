#!/usr/bin/env python3

import sys
import glob
import json
import aepathdisorder

in_folder = sys.argv[1]
out_file = sys.argv[2]

input_files = glob.glob(in_folder + "/*.fasta")
d = {}

for f in input_files:
    uniprot_id, taxon = aepathdisorder.map_protein_to_taxon(f)
    d[uniprot_id] = taxon

with open(out_file, 'w') as out_handle:
    json.dump(d, out_handle, indent=4, sort_keys=True)

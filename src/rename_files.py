#!/usr/bin/env python3

import sys
import os
import glob

# This is a tool to rename files according to a map.
# The file name must contain a key, that will be replaced by a value
# while everything else is kept.
# All happens in place
in_folder = sys.argv[1]
dict_file = sys.argv[2]


# Read the dictionary
d = {}
with open(dict_file) as f:
    for line in f:
        (key, val) = line.split()
        d[key] = val
# List file paths
input_files = glob.glob(in_folder + "/*")
print(f"{len(input_files)} files found")
# Transform file paths 
for f in input_files:
    for key in d.keys():
        if key in f:
            new_f = f.replace(key, d[key])
            os.rename(f, new_f)
            print(f"Moved {f} to {new_f}")
        else:
            pass


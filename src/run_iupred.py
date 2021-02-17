#!/usr/bin/env python3

import sys
import aepathdisorder

in_folder = sys.argv[1]
out_folder = sys.argv[2]

iupred_flags = ['short', 'long']

for flag in iupred_flags:
    aepathdisorder.run_iupred(in_folder,
                              out_folder,
                              flag,
                              iupred_path='./src/3rdparty/iupred/')

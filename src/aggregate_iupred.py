# /usr/bin/env/ python

import sys
import aepathdisorder

# This script takes disorder predictions from protein datasets and:
# - Calculates disorder fraction
# - Categorizes disorder type
# See the module for a more detailed description

in_folder = sys.argv[1]
out_folder = sys.argv[2]

# Load disorder prediction results
iupred_scores = aepathdisorder.load_iupred_scores(in_folder)

# Individually aggregate and save
for filename, score in iupred_scores:
    out_file = '_'.join(filename)+'_disfrac-cat.table'
    out_path = out_folder+out_file
    print('Aggregating disorder scores from {} ...'.format('_'.join(filename)))
    disfrac = aepathdisorder.calc_disfrac(score)
    cat_disfrac = aepathdisorder.categorize_disorder(score, disfrac)

    cat_disfrac.to_csv(path_or_buf=out_path, sep='\t')
    print('Results saved in {}'.format(out_path))
    print(79*'-')

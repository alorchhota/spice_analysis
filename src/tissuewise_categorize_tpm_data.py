'''
This script categorizes rpkm data.
It reads raw rpkm data, groups data according to tissue type and saves in files.
'''

import os
import sys
sys.path.append(os.path.dirname(__file__))
import pandas as pd
import argparse
import gzip
import csv
import re

print(os.path.basename(__file__) +  ': parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-rpkm', '--rpkm',
                    help='path to rpkm data file.',
                    default='/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/data/download/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
                    #default='results/transcript_rpkm.txt')
parser.add_argument('-ann', '--annotation',
                    help='path to sample annotation file.',
                    default='/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/data/download/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
parser.add_argument('-st', '--sample_start_col',
                    help='column number (1-based) from which sample data starts from.',
                    type=int,
                    default=3)
parser.add_argument('-skiprow',
                    help='number of rows to skip from the top in the rpkm file.',
                    type=int,
                    default=2)
parser.add_argument('-sep1',
                    help='separator in the input data file',
                    default='\t')
parser.add_argument('-sep2',
                    help='separator in the annotation file',
                    default='\t')
parser.add_argument('-sep3',
                    help='separator in the output data file',
                    default='\t')
parser.add_argument('-o', '--output_directory',
                    help='output directory to save categorized data.',
                    default='results/isoform')

args = parser.parse_args()

''' setting variables '''
rpkm_fn = args.rpkm
annotation_fn = args.annotation
out_dir = args.output_directory
start_col = args.sample_start_col
n_skip_row = args.skiprow
sep1 = args.sep1
sep2 = args.sep2
sep3 = args.sep3

if start_col < 1:
    raise Exception('sample start col must be a positive integer.\nstart_col:' + str(start_col))

if not os.path.exists(out_dir):
    raise Exception('output directory does not exist.\noutput_directory: ' + out_dir)

def split_with_quote(s, sep='\t', quote='"'):
    return [l for l in csv.reader([s], delimiter=sep, quotechar=quote)][0]


print('reading input files...')
annotation = pd.read_table(annotation_fn, sep=sep2)
# GTEx v8 data does not have SAMPID col, 1st col contains SAMPID
if 'SAMPID' not in annotation.columns.tolist():
    if annotation.columns[0] == 'Unnamed: 0':
        col_names = annotation.columns.tolist()
        col_names[0] = 'SAMPID'
        annotation.columns = col_names
    else:
        raise Exception('could not identify the sample id column.')

annotation = annotation.loc[:, ['SAMPID','SMTSD']]

print('samples without tissue type annotation ...')
for idx, ann in annotation.iterrows():
    if type(ann['SMTSD']) != str:
        print(ann['SAMPID'] + sep3 + str(ann['SMTSD']))


print('creating sample id to tissue type map ...')
sampleid_to_tissue_map = {row['SAMPID']:row['SMTSD'] for idx, row in annotation.iterrows()}

print('creating file for each tissue ...')
# find all tissue types
tissue_types = set(annotation['SMTSD'])
tissue_types = [t for t in tissue_types if type(t) is str]

print('processing header line of rpkm data ...')
# skip row
if rpkm_fn.lower().endswith('gz'):
    fh_rpkm = gzip.open(rpkm_fn, 'rt')
else:
    fh_rpkm = open(rpkm_fn, 'rt')

for i in range(n_skip_row):
    tmp = fh_rpkm.readline()

header_line = fh_rpkm.readline().rstrip()
col_titles = split_with_quote(header_line, sep = sep1, quote= '"')
initial_cols = list(range(start_col-1))

## there are 53 unnamed columns in transcript rpkm file.
unnamed_titles = [title for title in col_titles if len(title.strip())==0]
print("#Unnamed_cols: " + str(len(unnamed_titles)))

print('creating tissue to column indexes map ...')
tissue_to_col_map = {t:[] for t in tissue_types}
for i in range(len(col_titles)):
    sid = col_titles[i]
    if sid in sampleid_to_tissue_map:
        tissue = sampleid_to_tissue_map[sid]
        tissue_to_col_map[tissue].append(i)


print('creating one file per tissue ...')
def get_file_name(tissue, dir=out_dir):
    s = re.sub('[^0-9a-zA-Z\-]+', ' ', tissue)
    s = re.sub('[ ]+', '_', s.lstrip().rstrip())
    return os.path.join(dir, s + '.txt')

fh_tissues = {}
for t in tissue_types:
    fh_t = open(get_file_name(t), 'w')
    fh_tissues[t] = fh_t
    # write header line
    cols = [col_titles[i] for i in initial_cols + tissue_to_col_map[t]]
    header_text = sep3.join(cols)
    fh_t.write(header_text)
    # check if there exists multiple sets of data from same individual
    individuals = [ sid.split('-')[1] for sid in cols[start_col-1:]]
    if len(set(individuals)) != len(cols) - start_col + 1:
        print(t + ': there exists multiple sets of data from same individual!!!!!!!!!!!!!')

print('categorizing tissues ...')
line_count = 0
for line in fh_rpkm:
    line = line.rstrip()
    if len(line) == 0:
        continue
    line_count += 1
    if line_count % 1000 == 0:
        print(line_count)
    sp = split_with_quote(line, sep=sep1, quote='"')
    for t in tissue_types:
        cols = [sp[i] for i in initial_cols + tissue_to_col_map[t]]
        text = '\n' + sep3.join(cols)
        fh_tissues[t].write(text)

print('closing all the files ...')
for t in tissue_types:
    fh_tissues[t].close()

fh_rpkm.close()

print('done! see results in ' + out_dir)

#######################################################################
#### multiple tabs can be used as a delimiter #########################
#### following code checks if all multi-tabs are in same position #####
#### otherwise, the above code will produce error #####################
#######################################################################
def are_multi_tabs_consistent(rpkm_fn):
    fh = open(rpkm_fn, 'r')
    line = fh.readline().rstrip()
    col_titles = line.split('\t')
    unnamed_idx = [i for i in range(len(col_titles)) if col_titles[i]=='']

    def same_unnamed(line):
        line = line.rstrip()
        sp = line.split('\t')
        line_unnamed_idx = [i for i in range(len(sp)) if sp[i]=='']
        return unnamed_idx == line_unnamed_idx

    different_unnamed_all = [not same_unnamed(line) for line in fh]
    n_mismatch = sum(different_unnamed_all)
    if n_mismatch == 0:
        print('all tabs are positioned in consistently.')
    else:
        print('tabs are not consistent. data splitting code needs modification.')
    fh.close()
    return

#/bin/sh
# Do whatever is necessary to load your favorite version of python on the line(s) below

# Module package *always* sends its messages to stderr rather than stdout--even non-errors.
if module list 2>&1 | grep -q anaconda3/2020.11; then :
else
    module load anaconda3/2020.11
    conda activate py38 2>/dev/null
fi

# Now tell me the name of the python interpreter that you want to use, given the above
PYTHON=python3

# PYTHON GENIUSES DO NOT TOUCH ANYTHING BETWEEN HERE, AND WHERE IT SAYS "Put your Python code below ..."
#############################################################################################################

USAGE="$0 bla bla bla
PURPOSE: Run a python script included at the end of the Bourne Shell script."

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

sed '1,/^#Put/d' "$0" > $TMPDIR/$BASENAME.py
chmod +x $TMPDIR/$BASENAME.py

exec $PYTHON $TMPDIR/$BASENAME.py "$@"

conda deactivate py38

exit

#############################################################################################################
#Put your Python code below this line.

import csv
import argparse
from os import isatty
from sys import stdin
from io import StringIO

import numpy as np

from matplotlib import pyplot as plt

from sklearn.metrics import auc
from sklearn.metrics import ndcg_score
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import precision_recall_curve

class CompiledResult:
    
    def __init__(self):
        
        self.Ts = []
        self.Rs = []
        self.Ps = []
        self.node_pairs = []
        self.confidences = []
        self.scores = []
        self.orbit_pairs = []
        
    def add(self, T, R, P, node_pair, confidence, score, orbit_pair):
        
        self.Ts.append(int(T))
        self.Rs.append(int(R))
        self.Ps.append(float(P))
        self.node_pairs.append(node_pair)
        self.confidences.append(float(confidence))
        self.scores.append(float(score))
        self.orbit_pairs.append(orbit_pair)
        
        
    def get_ts(self):
        
        return self.Ts
    
    def get_scores(self):
        
        return self.scores
    
    def get_Rs(self):
        
        return self.Rs
    
    def __len__(self):
        
        return len(self.Ts)

def parse_csv(infile_, has_header):

    result = CompiledResult()
    
    file_rows = csv.reader(infile_, delimiter="\t")
    
    for r_count, file_row in enumerate(file_rows):
            
        if has_header and r_count == 0:
            continue

        T, R, _, P, node_pair = file_row[0].split(" ")
        confidence = float(file_row[1])
        score, _, orbit_pair = file_row[2].split(" ")

        result.add(T, R, P, node_pair, confidence, score, orbit_pair)
        
    return result

def read_file(fpath, has_header=False):
    
    result = None
    
    with open(fpath, 'r') as infile_:
        
        result = parse_csv(infile_, has_header)
        
        infile_.close()
        
    return result
    
def read_stdin(has_header=False):
    
    result = None
    
    infile_ = StringIO(stdin.read())
    
    result = parse_csv(infile_, has_header)
    
    return result
    
def get_auroc(result, out):
    
    assert isinstance(result, CompiledResult)
    
    y_pred = [1]*len(result)
    
    fprs, tprs, threshold = roc_curve(result.Ts, result.confidences)

    if out:

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)

        ax.plot(fprs, tprs)

        plt.savefig(f'{out}_roc.png')
    
    return roc_auc_score(result.Ts, result.confidences)

def get_aupr(result, out):
    
    assert isinstance(result, CompiledResult)
    
    y_pred = [1]*len(result)
    
    precision, recall, threshold = precision_recall_curve(result.Ts, result.confidences)

    if out:

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)

        ax.plot(recall, precision)

        plt.savefig(f'{out}_pr.png')
    
    return auc(recall, precision)

def get_ndcg(result):
    
    assert isinstance(result, CompiledResult)
    
    y_pred = [[1-confidence, confidence] for confidence in result.confidences]
    
    encoder = OneHotEncoder(sparse=False)

    reshaped_ts = np.array(result.Ts).reshape((len(result), 1))
    
    y_true = encoder.fit_transform(reshaped_ts)
    
    return ndcg_score(y_true, y_pred)

def get_parser():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--has-header', action='store_true')
    parser.add_argument('--fpath', type=str, help="Path to the tsv file")
    parser.add_argument('--out', type=str, help="Output path prefix for the images")
    
    return parser

if __name__ == "__main__":
    
    parser = get_parser()
    
    args = parser.parse_args()
    
    if args.fpath:
        
        result = read_file(args.fpath, args.has_header)

    else:

        result = read_stdin(args.has_header)


    print(f"AUROC: {get_auroc(result, args.out)}")

    print(f"AUPR: {get_aupr(result, args.out)}")

    print(f"NDCG: {get_ndcg(result)}")


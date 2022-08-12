import argparse
import numpy as np
from copy import deepcopy
import pandas as pd
import re
from Bio import SeqIO
import os, glob
from fuzzysearch import find_near_matches


class primer_list:
    def __init__(self, primer_seq, mean_score):
        """

        :rtype: object
        """
        self.primer_seq = primer_seq
        self.mean_score = mean_score

class exact_match_primer:
    def __init__(self, primer_seq, window_l, n_match):
        self.primer_seq = primer_seq
        self.window_l = window_l
        self.n_match = n_match


def import_files(target_path, off_path):
    targ_fasta_paths = glob.glob(os.path.join(target_path, '*.fna'))
    off_fasta_paths = glob.glob(os.path.join(off_path, '*.fna'))

    # create lists of sequences (target and off target)
    # and labels (0 = target, !0 = off target)
    aln = []
    labs = []
    for fasta_path in targ_fasta_paths:
        print(fasta_path)
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            labs.append(0)
            aln.append(str(seq_record.seq))
    for fasta_path in off_fasta_paths:
        print(fasta_path)
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            labs.append(1)
            aln.append(str(seq_record.seq))

    target = [i for i, val in enumerate(labs) if val == 0]  # target seq has label 0
    off = [i for i, val in enumerate(labs) if val != 0]
    return (aln, labs, target, off)


# select reference genome
# remove reference from target list
def select_ref(aln, target, index=0):
    ref_ind = target[index]
    ref = aln[ref_ind]
    aln_l = deepcopy(target)
    aln_l.pop(0)
    return(aln_l, ref)

def iter_window(ref, window, slide_ind, aln, aln_l):
    j = 0
    seq_l = []
    while j < len(ref) - window:
        primer_seq = ref[j:j + window]
        new_primer = detect_exact_seqs(primer_seq, aln_l, aln)
        seq_l.append(new_primer)
        j = j + slide_ind
    return(seq_l)

def detect_fuzzy_matches(seq_l, aln_l, off, aln, max_dif):
    # for low matches and no matches, get hamming dist
    primers_l = []
    matches=[] #number of exact matches
    target_snps=[] # avg hamming distance primer-target
    off_snps=[] # avg hamming distance primer-off
    if len(seq_l) > 0:
        for i in range(len(seq_l)):
            primer_seq = seq_l[i].primer_seq
            new_target = calc_nuc_dif(primer_seq, aln_l, aln, max_dif)
            new_off_target = calc_nuc_dif(primer_seq, off, aln, max_dif)
            primers_l.append(primer_seq)
            matches.append(seq_l[i].n_match / len(aln_l))
            target_snps.append(new_target.mean_score)
            off_snps.append(new_off_target.mean_score)
    return(primers_l, target_snps, matches, off_snps)

def write_out(primers_l, target_snps, matches, off_snps, filename):
    df = pd.DataFrame({'primers': primers_l,
                       'n_target_seq_match': matches,
                       'n_mismatches_primer_target': target_snps,
                       'n_mismatches_primer_off_target': off_snps})
    df.to_csv(filename, index=False)


# finds exact matches of window sized sequence fragment created from reference genome
# compared first to all targets. If an exact match is found in > threshold fraction of
# targets, and the fragment is not found in any of the off target sequences,
# if not, None is returned

def detect_exact_seqs(primer_seq, aln_l, aln):
    window_l = []
    for i in range(len(aln_l)):
        seq2n = aln_l[i]
        s2 = aln[seq2n]
        if re.search(primer_seq, s2):
            window_l.append(seq2n)
    le = len(window_l)
    new_match = exact_match_primer(primer_seq, window_l, le)
    return(new_match)


# calculates pairwise nuc differences (score) for a given primer
# and all other target or off target sequences
def calc_nuc_dif(primer_seq, aln_l, aln, max_dif):
    mismatch_l = []
    for i in range(len(aln_l)):
        seq2n = aln_l[i]
        s2 = aln[seq2n]
        matches = find_near_matches(primer_seq, s2, max_l_dist=max_dif)
        dist = []
        if matches != []:
            for match in matches:
                dist.append(match.dist)
            mismatch_l.append(min(dist))
        else:
            mismatch_l.append(max_dif)
    mean_score = np.mean(mismatch_l)
    primer_scores = primer_list(primer_seq, mean_score)
    return primer_scores


def main():
    parser = argparse.ArgumentParser(description='Choose primer locations')
    parser.add_argument('-target',
                        action="store",
                        help="directory containing target fastas",
                        dest="targ", type=str,
                        default='/Users/kristinaceres/CXL/ASFV/panaroo_blast_target_genes/D250R')
    parser.add_argument('-off',
                        action="store",
                        help="directory containing off target fastas",
                        dest="off", type=str,
                        default='/Users/kristinaceres/CXL/ASFV/off_target_genomes')
    parser.add_argument('-ref_index',
                        action="store",
                        help="index of reference genome in target directory",
                        dest="index", type=int,
                        default=0)
    parser.add_argument('-window',
                        action="store",
                        help="target region size",
                        dest="window", type=int,
                        default=300)
    parser.add_argument('-slide',
                        action="store",
                        help="step to slide window to look for next potential target region",
                        dest="slide", type=int,
                        default=10)
    parser.add_argument('-max_dif',
                        action="store",
                        help="maximum allowed differences between primers and targets",
                        dest="max_dif", type=int,
                        default=5)
    parser.add_argument('-out',
                        action="store",
                        help="filename for output primer sequences",
                        dest="out", type=str,
                        default="potential_primers_out.csv")
    args = parser.parse_args()
    aln, labs, target, off = import_files(target_path=args.targ, off_path=args.off)
    aln_l, ref = select_ref(aln, target, index=args.index)
    seq_l = iter_window(ref, window=args.window, slide_ind=args.slide, aln=aln, aln_l=aln_l)
    primers_l, target_snps, matches, off_snps = detect_fuzzy_matches(seq_l, aln_l, off, aln, max_dif=args.max_dif)
    write_out(primers_l, target_snps, matches, off_snps, args.out)


if __name__ == '__main__':
    main()

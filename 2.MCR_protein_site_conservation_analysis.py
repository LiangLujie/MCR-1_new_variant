from Bio import SeqIO
import sys
import numpy as np
import math

def loadAASeq(infile):
    seq = {}
    for i in SeqIO.parse(infile, 'fasta'):
        seq[i.id] = list(i.seq)
    return seq

def AAHash():
    AARow = {"A": 0,
             "R": 1,
             "N": 2,
             "D": 3,
             "C": 4,
             "Q": 5,
             "E": 6,
             "G": 7,
             "H": 8,
             "I": 9,
             "L": 10,
             "K": 11,
             "M": 12,
             "F": 13,
             "P": 14,
             "S": 15,
             "T": 16,
             "W": 17,
             "Y": 18,
             "V": 19,
             "-": 20}
    return AARow


def countFreq(seq, gaps, posi):
    if Seq[seq][posi] == '-':
        gaps += 1
    else:
        for key in Seq.keys():
            aa = Seq[key][posi]
            row = int(AA[aa])
            col = posi - gaps
            Count[row, col] += 1
    return gaps


def Quantile(Vec):
    Vlen = len(Vec)
    SortIndex = Vec.argsort()
    VRank = np.zeros(Vlen)
    forwardIndex = list()
    forwardValue = 0
    for i in range(0, Vlen):
        if Vec[SortIndex[i]] == forwardValue:
            forwardIndex.append(SortIndex[i])
            if i == Vlen - 1:
                for j in forwardIndex:
                    VRank[j] = i + 1
        else:
            if len(forwardIndex) >= 2:
                for j in forwardIndex:
                    VRank[j] = i
            else:
                if i > 0:
                    VRank[SortIndex[i - 1]] = i
            if i == Vlen - 1:
                VRank[SortIndex[i]] = i + 1
            forwardIndex = [SortIndex[i]]
            forwardValue = Vec[SortIndex[i]]
    PercValue = VRank * 100 / Vlen
    return PercValue


if __name__ == '__main__':
    in_msa = sys.argv[1]
    nseq = sys.argv[2]
    Seq = loadAASeq(in_msa)
    AA = AAHash()
    seqlenAll = len(Seq[nseq]) 
    seqlenInt = seqlenAll - Seq[nseq].count('-')  
    Count = np.zeros((21, seqlenInt), dtype=float) + 0.000001 
    gaps = 0
    for i in range(0, seqlenAll):
        gaps = countFreq(nseq, gaps, i)
    colSum = sum(Count[..., 0])
    CountRate = Count / colSum
    Shannon = -(CountRate * np.log2(CountRate)).sum(axis=0)
    Rindex = math.log2(21) - Shannon
    Rrelat = (Rindex - Rindex.min()) / (Rindex.max() - Rindex.min())
    Rquant = Quantile(Rindex)

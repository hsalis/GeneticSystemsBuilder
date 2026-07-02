#Commonly used Model Functions for use with the Non-Repetitive Parts Calculator (Maker Mode)
import re

def modelFunction_GC_Content(seq, GC_percentage_minimum, GC_percentage_maximum):
    GC_count = seq.count('G') + seq.count('C')
    return GC_percentage_minimum <= (GC_count * 100.0) / len(seq) <= GC_percentage_maximum

def modelFunction_Exclude_Sequences(seq, compiledRegExp):
    match = compiledRegExp.search(seq)
    if match is None:
        return (True, None)
    else:
        return (False, match.start() )

#### Model functions for designing non-repetitive coding sequences

AA2Codon = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '*' : ['TAA','TGA'],
    'O' : ['TAG']
}

Codon2AA = {}
for (AA, codonList) in list(AA2Codon.items()):
    for codon in codonList:
        Codon2AA[codon] = AA

def modelFunction_AASeq(seq, aa_seq):
    seqlen = len(seq)

    if seqlen % 3 == 0:
        AA_pos = seqlen // 3 - 1
        lastCodon = seq[-3:]
        if Codon2AA[lastCodon] == aa_seq[AA_pos]:
#            print "AA %s matches expected AA %s" % (Codon2AA[lastCodon], aa_seq[AA_pos])
#            print "AA matches at position %s" % seqlen
            return (True, None)
        else:
#            print "AA %s does not match expected AA %s" % (Codon2AA[lastCodon], aa_seq[AA_pos])
#            print "Going back to position %s" % str(seqlen - 1)
            return (False, seqlen - 1 )
    else:
        return (True, None)







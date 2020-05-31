#!/usr/bin/python3

'''
This script takes a file with the genes expressed in a tissue and then returns only those from a list of mean and std. dev. values.

e.g.
getGenes.py /cygdrive/c/Users/evansp1/Dropbox/dropbox_evolution_paper/ExpressionData/Adipose_Subcutaneous_euroids_expr_expressed0.1in0.1pop.txt /cygdrive/c/Users/evansp1/Dropbox/dropbox_evolution_paper/ExpressionData/Adipose_Subcutaneous_rpkm_meanRPKM.txt > outputfile
'''

import sys

#main
if len(sys.argv) != 3:
    sys.stderr.write(__doc__+"\n")
    sys.exit(0)

#import and process gene list
genedict = {}
fin = open(sys.argv[1], 'r')
head = fin.readline()
for line in fin:
    line = line.strip().split()
    genedict[line[0]] = 1
fin.close()

#import and process expression mean file
fin = open(sys.argv[2], 'r')
head = fin.readline()
head = head.strip()
print(head)
for line in fin:
    line = line.strip().split()
    if line[0] in genedict:
        print("\t".join(line))
fin.close()
sys.exit(0)

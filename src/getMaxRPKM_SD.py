#!/usr/bin/python3

'''
This script takes a list of files with mean and standard deviation RPKM values for each gene and returns the tissue that has the maximum mean RPKM and its corresponding standard deviation

e.g. ./getMaxRPKM_SD.py filelist > outputfile
'''

import sys

#main
if len(sys.argv) != 2:
    sys.stderr.write(__doc__+"\n")
    sys.exit(0)

#import and process file list
fin = open(sys.argv[1], 'r')
count = 0
filelist = []
for line in fin:
    print(count)
    line = line.strip()
    filelist[count] = line
    count += 1
fin.close()

#open and process files
genedict = {}
for file in filelist:
    fin = open(file, 'r')
    head = fin.readline()
    for line in fin:
        line = line.strip().split()
        if line[0] in genedict:
            currentvalue = genedict[line[0]]
            if currentvalue[0] < line[1]:
                genedict[line[0]] = (line[1], line[2])
        else:
            genedict[line[0]] = (line[1], line[2])
    fin.close()

#print values
for gene in genedict:
    print(gene+"\t"+"\t".join(genedict[gene]))

sys.exit(0)

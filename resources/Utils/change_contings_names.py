#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import csv

inmap=sys.argv[1]
inputfile=open(inmap)
outmap=sys.argv[2]
contigs_dict={}
contig_num=0
with open(outmap,"w", newline="") as f:
    writer=csv.writer(f,delimiter='\t')
    for line in csv.reader(inputfile, delimiter='\t'):
        contig=line[0]
        if (contig in contigs_dict):
            line[0]=contigs_dict.get(contig)
            writer.writerow(line)
        else:
            contig_num+=1
            contigs_dict[contig]=contig_num
            line[0]=contigs_dict.get(contig)
            writer.writerow(line)

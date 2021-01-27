#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

# This scripts creates the parameter file for the convertf utility

if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")
    # initializing variables name from snakemake objects
    out_dir=snakemake.wildcards.backbone
    target=snakemake.wildcards.target
    outfile=snakemake.output[0]
    # creating the parameter file
    with open(outfile, "w") as f:
        out_line="genotypename:\tresults/tmp/"+out_dir+"/"+target+".ped\n"
        f.write(out_line)
        out_line="snpname:\tresults/tmp/"+out_dir+"/"+target+".map\n"
        f.write(out_line)
        out_line="indivname:\tresults/tmp/"+out_dir+"/"+target+".ped\n"
        f.write(out_line)
        out_line="outputformat:\tEIGENSTRAT\n"
        f.write(out_line)
        out_line="genotypeoutname:\tresults/tmp/"+out_dir+"/"+target+".eigenstratgeno\n"
        f.write(out_line)
        out_line="snpoutname:\tresults/tmp/"+out_dir+"/"+target+".snp\n"
        f.write(out_line)
        out_line="indivoutname:\tresults/tmp/"+out_dir+"/"+target+".ind\n"
        f.write(out_line)
        out_line="familynames:\tNO"
        f.write(out_line)

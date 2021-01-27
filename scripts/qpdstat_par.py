#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script creates the parameter file for the convertf utility

if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")
    # Defining variables for convenience via accessing snakemake objects
    out_dir=snakemake.wildcards.backbone
    target=snakemake.wildcards.target
    outfile=snakemake.output[0]
    size=snakemake.params.blocksize
    # Writing the parameter file
    with open(outfile, "w") as f:
        out_line="genotypename:\tresults/tmp/"+out_dir+"/"+target+".eigenstratgeno\n"
        f.write(out_line)
        out_line="snpname:\tresults/tmp/"+out_dir+"/"+target+".snp\n"
        f.write(out_line)
        out_line="indivname:\tresults/tmp/"+out_dir+"/"+target+"_correct.ind\n"
        f.write(out_line)
        out_line="poplistname:\tresults/tmp/"+out_dir+"/"+target+"_keep_list.txt\n"
        f.write(out_line)
        out_line="blgsize:\t"+size+"\n"
        f.write(out_line)

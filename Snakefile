#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
import pandas as pd
startTime = datetime.now()

# Loading the config file
configfile: "arctic_dogs_fox_ref_simple_bb_config.yaml"

# Defining container
container: "docker://continuumio/miniconda3"


# Parsing the file containing the targets and creating a list of targets name
TARGET_LIST=list(pd.read_csv(config["target_list"]["tl"],header=None,skip_blank_lines=True)[0])

# Defining parameters/constants:
# 1) blg is the length (in Morgan) of the blocks for the jackknifing procedure
# employed by qpdstats (Admixtools) to calculate the standard error
# around the mean estimate of th D-statistics
# 2) n_sd_for_outliers is the number of standard deviation around the mean that
# is used by Admixture Graph (R package) to define the span of the C.I. around
# each of the observed d-statistic values and consequently to determine whether
# or not a fitted value represents an outlier under the model in question.
blg="0.005"
n_sd_for_outliers=6

######### Modules #########
include: "rules/subsetting.smk"
include: "rules/dstats.smk"
include: "rules/modelling.smk"

######### Target rules #########
rule all:
    input:
        "results/tmp/subset_all.done",
        "results/tmp/dstats_all.done",
        "results/tmp/graphs_all.done",
        "results/tmp/heat_map_all.done",
        expand(
        "results/Admixture_Graphs/full_heat_map_{analysis_name}.pdf",
        analysis_name=config["analysis_name"])

print(datetime.now() - startTime)

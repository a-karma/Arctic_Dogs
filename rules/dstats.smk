############################
### Step 2: D-statistics ###
###         (5 rules)    ###
############################

# This rule creates a configuration file needed by Admixtools convertf utility
rule convertf_par:
    input:
        bb_path=lambda wildcards: config["backbones"][wildcards.backbone]
    log:
        "results/log/{backbone}/convertf_par_{target}.log"
    output:
        temp("results/tmp/{backbone}/convertf_{target}.par")
    script:
        "../scripts/convertf_par.py"

# This rule uses the parameter file created by the previous rule
# and converts each .ped and .map file (produced by the
# subsetting.smk module) into eigenstraat format
rule convertf:
    input:
        par_file="results/tmp/{backbone}/convertf_{target}.par",
        in_file1="results/tmp/{backbone}/{target}.ped",
        in_file2="results/tmp/{backbone}/{target}.map"
    log:
        "results/log/{backbone}/convertf_run_{target}.log"
    output:
        temp(multiext("results/tmp/{backbone}/{target}",".snp",".eigenstratgeno",".ind")),
        indfile=temp("results/tmp/{backbone}/{target}_correct.ind")
    params:
        ind="results/tmp/{backbone}/{target}.ind"
    conda:
        "../envs/admixtool.yaml"
    shell:
        "convertf -p {input.par_file} 2>{log};"
        """awk -v OFS='\\t' '{{print "\\t",$1,$2,$1}}' {params.ind} > {output.indfile}"""

# This rule creates a parameter file required by the program qpdstat (Admixtools)
rule qpdstat_par:
    input:
        bb_path=lambda wildcards: config["backbones"][wildcards.backbone]
    log:
        "results/log/{backbone}/qpdstat_par_{target}.log"
    output:
        temp("results/tmp/{backbone}/qpdstat_{target}.par")
    params:
        blocksize=blg
    script:
        "../scripts/qpdstat_par.py"

# This rule uses the newly converted files in eigenstraat format and the Admixtools suit
# to calculate d-statistics for any possible quadruplet for the set {backbone + target}
rule qpdstat:
    input:
        pop_list="results/tmp/{backbone}/{target}_keep_list.txt",
        indfile="results/tmp/{backbone}/{target}_correct.ind",
        snpfile="results/tmp/{backbone}/{target}.snp",
        genofile="results/tmp/{backbone}/{target}.eigenstratgeno",
        par="results/tmp/{backbone}/qpdstat_{target}.par"
    log:
        "results/log/{backbone}/qpdstat_run_{target}.log"
    output:
        "results/D-stats/{backbone}/dstats_table_{target}.tsv"
    conda:
        "../envs/admixtool.yaml"
    shell:
        "qpDstat -p {input.par} > {log}; "
        """grep "result:" {log} | """
        """awk '{{OFS="\\t" ; print $2,$3,$4,$5,$6,$7}}' """
        " | cat resources/header.txt - > {output}"

# This final rule requests that the d-statistic calculations are
# done for any combination of wildcards: target/backbone and
# that were stored in the appropriate file
rule par_agg:
    input:
        expand("results/D-stats/{backbone}/dstats_table_{target}.tsv", target=TARGET_LIST, backbone=config["backbones"])
    log:
        "results/log/dstats_all.log"
    output:
        temp("results/tmp/dstats_all.done")
    shell:
        "touch {output} 2> {log}"

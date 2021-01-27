##########################
### Step 1: Subsetting ###
###         (3 rules)  ###
##########################

# This rule creates a file containing the list of samples we want to analyse:
# one backbone + a single target
rule k_list_one:
    input:
        bb_path=lambda wildcards: config["backbones"][wildcards.backbone]
    log:
        "results/log/{backbone}/keep_list_{target}.log"
    output:
        temp("results/tmp/{backbone}/{target}_keep_list.txt")
    shell:
        """awk '{{OFS="\\t"; print $1,$1 }}' {input.bb_path} > {output}; """
        """echo {wildcards.target} | awk '{{OFS="\\t"; print $1,$1 }}' >> {output} 2>{log}"""

# This rule subsets the master ped file extracting the individuals required
# for the analysis according to the list created by the previous rule.
rule plink_subset_one:
    input:
        "results/tmp/{backbone}/{target}_keep_list.txt"
    log:
        "results/log/{backbone}/plink_{target}.log"
    output:
        temp(multiext("results/tmp/{backbone}/{target}",".ped",".log",".map",".nosex"))
    params:
        m_set_prefix=config["master_file_set"]["master_no_ext"],
        target_prefix=lambda w, output: os.path.splitext(output[0])[0],
        chrset=config["chr_set"]
    conda:
        "../envs/plink.yaml"
    shell:
        "plink"
        " --file {params.m_set_prefix}"
        " --keep {input}"
        " --{params.chrset}"
        " --allow-no-sex"
        " --recode"
        " --out {params.target_prefix} ;"
        "cp results/tmp/{wildcards.backbone}/{wildcards.target}.log {log}"

# This final rule request that the subsetting is done for any combination of
# wildcards: target/backbone
rule subset_agg:
    input:
        expand("results/tmp/{backbone}/{target}.ped", target=TARGET_LIST, backbone=config["backbones"])
    log:
        "results/log/subset_all.log"
    output:
        temp("results/tmp/subset_all.done")
    shell:
        "touch {output} 2> {log}"

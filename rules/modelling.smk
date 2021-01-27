###################################
### Step 3: Admixture modelling ###
###            (7 rules)        ###
###################################

# This rule takes the d-stat calculations (performed by
# the dstats.smk module) as input and produces three output files
# for any combination of one backbone + a single targetself.
# Outputs consist of:
#   1) a text file containing the summary of each fitted model
#   given by the package admixture graphs
#   2) a pdf file containing a schematic representation of each
#   fitted graphs and how it fits the data
#   3) a temporary file with the number of outliers under each graph
rule modelling_single_sample:
    input:
        d_table="results/D-stats/{backbone}/dstats_table_{target}.tsv",
        pop_names=lambda wildcards: config["backbones"][wildcards.backbone]
    log:
        model_fitting_log="results/log/{backbone}/admixture_modelling_{target}.log"
    output:
        mod_summary="results/Admixture_Graphs/{backbone}/summary_{target}.txt",
        graphs="results/Admixture_Graphs/{backbone}/graphs_{target}.pdf",
        outliers=temp("results/tmp/{backbone}/outliers_{target}.tsv")
    params:
        sample="{target}",
        bb_prefix="{backbone}",
        models_to_fit=config["model_scripts"]["adm_models_to_fit"],
        sigma=n_sd_for_outliers
    conda:
        "../envs/adm.yaml"
    script:
        "../scripts/adm_modelling.R"

# This rule requests that the model fitting is done for
# any combination of wildcards: target/backbone
rule models_agg:
    input:
        expand("results/Admixture_Graphs/{backbone}/graphs_{target}.pdf",
        target=TARGET_LIST, backbone=config["backbones"])
    log:
        "results/log/graphs_all.log"
    output:
        temp("results/tmp/graphs_all.done")
    shell:
        "touch {output} 2> {log}"

# This rule gather all the temporary files containing the number of outliers under each model
# for a given combination of target/backbone and creates the full outliers table
rule full_outliers_tab:
    input:
        full_t_h="resources/full_table_header.txt",
        outl_files=
        expand("results/tmp/{backbone}/outliers_{target}.tsv",
        target=TARGET_LIST,backbone=config["backbones"])
    log:
        expand("results/log/outliers_full_table_{analysis_name}.log",
        analysis_name=config["analysis_name"])
    output:
        expand("results/Admixture_Graphs/outliers_full_table_{analysis_name}.tsv",
        analysis_name=config["analysis_name"])
    shell:
        "cat {input.full_t_h} {input.outl_files} > {output} 2>{log}"

# This rule extracts the outliers for every target and a given backbone
# producing an outliers table that serves as input for the following rule
rule single_backbone_outl_tab:
    input:
        expand("results/Admixture_Graphs/outliers_full_table_{analysis_name}.tsv",
        analysis_name=config["analysis_name"])
    log:
        "results/log/{backbone}/outliers_table.tsv"
    output:
        "results/Admixture_Graphs/{backbone}/outliers_table.tsv"
    shell:
        """
        grep "{wildcards.backbone}" {input} | cat resources/full_table_header.txt - |
        awk '{{ for (i=2; i<NF; i++) printf $i "\t"; print $NF}}' > {output} 2>{log}
        """

# This rule generates an outlier heat map for a given backbones
rule single_heat_map:
    input:
        bb_outl_tab="results/Admixture_Graphs/{backbone}/outliers_table.tsv",
        bb_list=lambda wildcards: config["backbones"][wildcards.backbone]
    log:
        hm_log="results/log/{backbone}/heat_map.log"
    output:
        hm_s_out=report("results/Admixture_Graphs/{backbone}/heat_map.pdf",
        caption="../report/simple_heat_map.rst",
        category="Backbone specific analysis")
    conda:
        "../envs/adm.yaml"
    script:
        "../scripts/heat_map_v2.R"

# This rule makes sure that an heat map will be generated for every backbone
rule heat_map_all:
    input:
        expand("results/Admixture_Graphs/{backbone}/heat_map.pdf",
        backbone=config["backbones"])
    log:
        hm_all_log="results/log/heat_map_all.log"
    output:
        temp("results/tmp/heat_map_all.done")
    shell:
        "touch {output} 2> {log}"

# This final rule generates an aggregated heat map based on the number of outliers
# under each model averaged across backbones
rule full_heat_map:
    input:
        outl_full_tab=expand(
        "results/Admixture_Graphs/outliers_full_table_{analysis_name}.tsv",
        analysis_name=config["analysis_name"])
    log:
        full_hm_log=expand("results/log/full_heat_map_{analysis_name}.log",
        analysis_name=config["analysis_name"])
    output:
        hm_f_out=report(
        expand("results/Admixture_Graphs/full_heat_map_{analysis_name}.pdf",
        analysis_name=config["analysis_name"]),
        caption="../report/full_heat_map.rst",
        category="Main Results"),
        graphs=report(
        expand("results/Admixture_Graphs/adm_graphs_schematics_{analysis_name}.pdf",
        analysis_name=config["analysis_name"]),
        caption="../report/admixture_models.rst",
        category="Main Results")
    conda:
        "../envs/adm.yaml"
    params:
        models_schematics=config["model_scripts"]["adm_graphs"]
    script:
        "../scripts/full_heat_map_v2.R"

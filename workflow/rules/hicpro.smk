import os

rule hicpro_align:
    output:
        bwt2_R1_bwt2merged_bam=temp("Results/{sample}/bowtie_results/bwt2/{sample}/{sample}_R1_mm10.bwt2merged.bam"),
        bwt2_R2_bwt2merged_bam=temp("Results/{sample}/bowtie_results/bwt2/{sample}/{sample}_R2_mm10.bwt2merged.bam")
    conda:
        config["conda_env"]
    params:
        rawdata_dir = os.path.abspath(config["input_dir"]),
        #rawdata_dir = lambda wildcards: os.path.abspath(config["input_dir"]) / wildcards.sample,  # 将相对路径转换为绝对路径
        config_hicpro = "./config/config-hicpro.txt",
        results_dir = os.path.abspath("Results")
    shell:
        """
        # mapping
        HiC-Pro -c {params.config_hicpro} -i {params.rawdata_dir}/{wildcards.sample} -o Results/{wildcards.sample} -s mapping
        
        # deleting intermediate
        sleep 10
        echo "deleting bwt2_global and bwt2_local of {wildcards.sample}"
        rm -r {params.results_dir}/{wildcards.sample}/bowtie_results/bwt2_global/{wildcards.sample} {params.results_dir}/{wildcards.sample}/bowtie_results/bwt2_local/{wildcards.sample}
        """

rule hicpro_proc_hic:
    input:
        bwt2_R1_bwt2merged_bam="Results/{sample}/bowtie_results/bwt2/{sample}/{sample}_R1_mm10.bwt2merged.bam",
        bwt2_R2_bwt2merged_bam="Results/{sample}/bowtie_results/bwt2/{sample}/{sample}_R2_mm10.bwt2merged.bam"
    output:
        bwt2_pairs_bam = "Results/{sample}/bowtie_results/bwt2/{sample}/{sample}_mm10.bwt2pairs.bam",
        DumpPairs=temp("Results/{sample}/hic_results/data/{sample}/{sample}_mm10.bwt2pairs.DumpPairs"),
        validPairs="Results/{sample}/hic_results/data/{sample}/{sample}_mm10.bwt2pairs.validPairs"
    conda:
        config["conda_env"]
    params:config_hicpro = "./config/config-hicpro.txt"
    shell:
        """
        HiC-Pro -c {params.config_hicpro} -i Results/{wildcards.sample}/bowtie_results/bwt2 -o Results/{wildcards.sample} -s proc_hic
        """

rule hicpro_merge:
    input:
        bwt2_pairs_bam = "Results/{sample}/hic_results/data/{sample}/{sample}_mm10.bwt2pairs.validPairs"
    output:
        validPairs="Results/{sample}/hic_results/data/{sample}/{sample}.allValidPairs"
    conda:
        config["conda_env"]
    params:
        config_hicpro = "./config/config-hicpro.txt",
        collect_dir = "Results/00collect",
        validpairs = lambda wildcards: os.path.abspath(f"Results/{wildcards.sample}/hic_results/data/{wildcards.sample}/{wildcards.sample}.allValidPairs")
    shell:
        """
        HiC-Pro -c {params.config_hicpro} -i Results/{wildcards.sample}/hic_results/data -o Results/{wildcards.sample} -s merge_persample
        mkdir -p {params.collect_dir}
        ln -sf {params.validpairs} {params.collect_dir}
        """

rule hicpro_quality_checks:
    input:
        validPairs="Results/{sample}/hic_results/data/{sample}/{sample}.allValidPairs"
    output:
        fig="Results/{sample}/hic_results/pic/{sample}/plotHiCContactRanges_{sample}.pdf"
    conda:
        config["conda_env"]
    params:config_hicpro = "./config/config-hicpro.txt"
    shell:
        """
        HiC-Pro -c {params.config_hicpro} -i Results/{wildcards.sample}/hic_results/data -o Results/{wildcards.sample} -s quality_checks
        """

rule summary:
    input:
        hicpro_finish=expand("Results/{sample}/hic_results/data/{sample}/{sample}.allValidPairs", sample = samples)
    output:
        summary="Results/summary.txt"
    conda:
        config["conda_env"]
    params:
        sample_list=samples,
        results_dir = os.path.abspath("Results")
    script:
        "../scripts/summary.py"



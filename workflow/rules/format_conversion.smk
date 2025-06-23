rule allVP2hic:
    input:
        allValidPairs = "Results/{sample}/hic_results/data/{sample}/{sample}.allValidPairs"
    output:
        hic = "Results/00collect/{sample}.allValidPairs.hic"
    conda:
        config["conda_env"]
    group: "format_conversion"
    params:
        hic2juice = config["hic2juice"],
        genomesize = config["genomesize"],
        juicertool = config["juicertool"],
        enzyme_bed = config["enzyme_bed"],
        output_dir = "Results/00collect"
        #output_dir= lambda wildcards: f"Results/{wildcards.sample}/hic_results/data/{wildcards.sample}"
    shell:
        """
        {params.hic2juice} -i {input.allValidPairs} -g {params.genomesize} -j {params.juicertool} -r {params.enzyme_bed} -o {params.output_dir}
        """

rule hic2cool:
    input:
        hic = "Results/00collect/{sample}.allValidPairs.hic"
    output:
        cool = "Results/00collect/{sample}_{resolution}.cool"
    conda:
        config["conda_env"]
    group: "format_conversion"
    params:
        cool_name = "Results/00collect/{sample}.cool" # hic的输出会自动在给定的name后加上resolution
        #cool_name = "Results/{sample}/hic_results/data/{sample}/{sample}.cool"
    shell:
        """
        # convert
        hicConvertFormat -m {input.hic} --inputFormat hic --outputFormat cool -r {wildcards.resolution} -o {params.cool_name}
        # balance
        cooler balance {output.cool}
        """
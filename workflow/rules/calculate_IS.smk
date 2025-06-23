rule calculate_IS:
    input:
        cool = "Results/00collect/{sample}_{resolution}.cool"
    output:
        IS = "Results/00collect/{sample}_{resolution}.insulation_scores.tsv"
    conda:
        config["conda_env"]
    params:
        windows=config["windows"]
    shell:
        """
        cooltools insulation -p 20 -o {output.IS} --chunksize 20000000 {input.cool} {params.windows}
        """

# rule IS_visulization:
#     input:
#         IS = "Results/{sample}/hic_results/data/{sample}/{sample}_{resolution}.insulation_scores.tsv"
#     output:
#         filtered_bed = "Results/{sample}/hic_results/data/{sample}/{sample}_{resolution}_filter.bed",
#         bw = "Results/{sample}/hic_results/data/{sample}/{sample}_{resolution}.IS.bw"
#     group:"processing_group"
#     params:
#         genomesize=config["genomesize"]
#     shell:
#         """
#         awk '$1 ~ /^[1-9]$|^1[0-9]$|^[XY]$/ && $6 != "nan" {print "chr"$1, $2, $3, $6}' {input.IS} > {output.filtered_bed}
#         bedGraphToBigWig {output.filtered_bed} {params.genomesize} {output.bw}
#         """
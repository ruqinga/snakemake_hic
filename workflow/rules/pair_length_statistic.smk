rule pair_length_statistic:
    input:
        hicpro_finish=expand("Results/{sample}/hic_results/data/{sample}/{sample}.allValidPairs", sample = samples)
    output:
        cumulative_plot="Results/statistic/cumulative_plot.png",
        distance_histogram="Results/statistic/distance_histogram.png",
        pair_statistic="Results/statistic/pair_statistic.tsv"
    params:
        ext=".allValidPairs",
        input_dir="Results/00collect"
    conda:
        config["conda_env"]
    group: "global_process"
    shell:
        """
        python workflow/scripts/pair_statistic.py  --indir {params.input_dir} --outplot {output.cumulative_plot} --outbar {output.distance_histogram} --outtsv {output.pair_statistic}
        """
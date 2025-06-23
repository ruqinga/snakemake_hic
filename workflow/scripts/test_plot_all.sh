#!/bin/bash

# 定义样本列表
samples=("E11" "E13_M" "E13_F")
plot_IS="/home_data/home/slst/leixy2023/data/project/ref_data/hic_PGC_PRJNA484983/workflow/scripts/plot_IS.py"
plot_CTCF="/home_data/home/slst/leixy2023/data/project/ref_data/hic_PGC_PRJNA484983/workflow/scripts/plot_CTCF_stackup.py"

# 遍历每个样本
for sample in "${samples[@]}"; do
    # 打印正在处理的样本
    echo "Processing sample: $sample"

    # 定义基础路径和文件名
    cool_file="${sample}_40000.cool"
    insulation_file="${sample}.insulation_scores.tsv"
    output_plot="./fig/plot_IS_${sample}_l90.png"
    ctcf_stack_output="./fig/plot_CTCF_stack_${sample}.png"
    bw_file="/home_data/home/slst/leixy2023/data/project/CTCF/20250415_all_from_huanhuan/Results/06_visualization/${sample}_rep1_unique_100.bw"

    # 执行 IS 图像绘制
    python "$plot_IS" -c "$cool_file" -is "$insulation_file" -o "$output_plot" -w 1000000 -s 10_500_000 -l 90
    if [ $? -eq 0 ]; then
        echo "plot IS finish for $sample"
    else
        echo "Error in plotting IS for $sample"
        continue
    fi

    # 执行 CTCF stack 绘制
    python "$plot_CTCF" -is "$insulation_file" -o "$ctcf_stack_output" -w 1000000 -r 4000 -bw "$bw_file"
    if [ $? -eq 0 ]; then
        echo "plot CTCF stack finish for $sample"
    else
        echo "Error in plotting CTCF stack for $sample"
        continue
    fi
done

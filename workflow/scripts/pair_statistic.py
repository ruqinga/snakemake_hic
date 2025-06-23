#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="分析 allValidPairs 文件：统计 intra/inter-chr 配对，并绘制距离分析图")
    parser.add_argument('--indir', required=True, help='输入目录')
    parser.add_argument('--ext', default='.allValidPairs', help='文件后缀，默认 .allValidPairs')
    parser.add_argument('--outplot', default='cumulative_plot.png', help='输出累计图像')
    parser.add_argument('--outbar', default='distance_histogram.png', help='输出区间柱状图')
    parser.add_argument('--outtsv', default='pair_statistic.tsv', help='输出 TSV 统计文件')
    return parser.parse_args()

def valid_chromosomes():
    return [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]

# 距离区间
BINS = [0, 5_000, 10_000, 20_000, 50_000, 100_000, 250_000, 500_000, 1_000_000, 2_500_000, float('inf')]
LABELS = ["<5k", "5-10k", "10-20k", "20-50k", "50-100k", "100-250k", "250-500k", "500k-1M", "1M-2.5M", "≥2.5M"]

def process_file(filepath, chromosomes):
    df = pd.read_csv(filepath, sep='\t', header=None, usecols=[1, 2, 4, 5])
    df.columns = ['chr1', 'pos1', 'chr2', 'pos2']

    # 过滤目标染色体
    df = df[df['chr1'].isin(chromosomes) & df['chr2'].isin(chromosomes)]

    same_chr_mask = df['chr1'] == df['chr2']
    same_chr_count = same_chr_mask.sum()
    diff_chr_count = len(df) - same_chr_count

    diffs = np.abs(df.loc[same_chr_mask, 'pos1'] - df.loc[same_chr_mask, 'pos2'])

    # 分箱统计
    bin_counts = pd.cut(diffs, bins=BINS, labels=LABELS, right=False).value_counts(sort=False)

    return same_chr_count, diff_chr_count, diffs, bin_counts

def plot_cumulative_diffs(diff_dict, outplot):
    plt.figure(figsize=(8, 6))
    for fname, diffs in diff_dict.items():
        log_diffs = np.log10(diffs + 1)
        sorted_vals = np.sort(log_diffs)
        cdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        plt.plot(sorted_vals, cdf, label=os.path.basename(fname), linewidth=1)

    plt.xlabel('log10(Absolute Position Difference + 1)')
    plt.ylabel('Cumulative Fraction')
    plt.title('Intra-chromosomal Distance (Log-scaled CDF)')
    plt.legend(fontsize='small')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outplot)
    print(f"[✓] 累积分布图已保存：{outplot}")

def plot_histogram(summary_df, outbar):
    plt.figure(figsize=(10, 6))
    bin_labels = LABELS
    bar_data = summary_df[bin_labels].astype(int)

    bar_data.index = summary_df['File']

    ax = bar_data.plot(kind='bar', figsize=(10, 6), colormap='tab20')

    # 在每个条形上显示数值
    for p in ax.patches:
        ax.annotate(f'{p.get_height()}',
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    xytext=(0, 5),
                    textcoords='offset points',
                    ha='center', va='bottom', fontsize=10)

    plt.ylabel('Read Pair Count')
    plt.title('Intra-chromosomal Distance Distribution')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outbar)
    print(f"[✓] 柱状图已保存：{outbar}")

# def plot_histogram(summary_df, outbar):
#     plt.figure(figsize=(10, 6))
#     bin_labels = LABELS
#     bar_data = summary_df[bin_labels].astype(int)
#
#     bar_data.index = summary_df['File']
#     bar_data.T.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab20')
#     plt.ylabel('Read Pair Count')
#     plt.title('Intra-chromosomal Distance Distribution')
#     plt.xticks(rotation=45)
#     plt.tight_layout()
#     plt.savefig(outbar)
#     print(f"[✓] 柱状图已保存：{outbar}")

def main():
    args = parse_args()
    chromosomes = valid_chromosomes()

    summary_rows = []
    diff_dict = {}

    for fname in os.listdir(args.indir):
        if not fname.endswith(args.ext):
            continue
        fpath = os.path.join(args.indir, fname)
        try:
            same_count, diff_count, diffs, bin_counts = process_file(fpath, chromosomes)
            summary_row = {
                'File': fname,
                'Intra': same_count,
                'Inter': diff_count
            }
            summary_row.update(bin_counts.to_dict())
            summary_rows.append(summary_row)
            if same_count > 0:
                diff_dict[fname] = diffs
            print(f"{fname}: intra = {same_count}, inter = {diff_count}")
        except Exception as e:
            print(f"[警告] 文件处理失败 {fname}: {e}")

    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.fillna(0)
    summary_df.to_csv(args.outtsv, sep='\t', index=False)
    print(f"[✓] 统计表已保存：{args.outtsv}")

    if diff_dict:
        plot_cumulative_diffs(diff_dict, args.outplot)
        plot_histogram(summary_df, args.outbar)

if __name__ == '__main__':
    main()

import bbi
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from skimage.filters import threshold_li, threshold_otsu
import argparse


def plot_boundary_signal(insulation_table_path, resolution, windows, bw_path, output_file_name,method):
    insulation_table = pd.read_csv(insulation_table_path, sep='\t')

    thresholds= {}

    for w in windows:
        boundary_strength = insulation_table[f'boundary_strength_{w}'].dropna().values
        if method == 'otsu':
            thresholds[w] = threshold_otsu(boundary_strength)
        elif method == 'li':
            thresholds[w] = threshold_li(boundary_strength)

    w = windows[0]
    top_boundaries = insulation_table[insulation_table[f'boundary_strength_{w}'] >= thresholds[w]]

    if not top_boundaries['chrom'].iloc[0].startswith('chr'):
        top_boundaries.loc[:, 'chrom'] = 'chr' + top_boundaries['chrom'].astype(str)

    flank = 50000
    nbins = 100
    stackup = bbi.stackup(
        bw_path,
        top_boundaries.chrom,
        top_boundaries.start + resolution // 2 - flank,
        top_boundaries.start + resolution // 2 + flank,
        bins=nbins,
        exact=True
    )

    plt.switch_backend('Agg')  # 设置 Matplotlib 使用非交互式后端
    f, ax = plt.subplots(figsize=[7, 5])
    ax.plot(np.nanmean(stackup, axis=0))
    ax.set(
        xticks=np.arange(0, nbins + 1, 10),
        xticklabels=(np.arange(0, nbins + 1, 10) - nbins // 2) * flank * 2 / nbins / 1000,
        xlabel='Distance from boundary, kbp',
        ylabel='CTCF enrichment'
    )
    output_fig = output_file_name + "_" + method + ".png"
    plt.savefig(output_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {output_fig}")

    # 导出表格
    output_table = output_file_name + "_" + method + ".bed"
    top_boundaries_first_three_columns = top_boundaries.iloc[:, :3]
    top_boundaries_first_three_columns.to_csv(output_table, sep='\t', index=False, header=False)
    print(f"Saved boundaries table to: {output_table}")


# 命令行执行支持
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot mean signal around insulation boundaries")

    parser.add_argument('-is', '--insulation_file', type=str, required=True,
                        help="Path to the insulation scores file (TSV format)")
    parser.add_argument('-o', '--output_file_name', type=str, required=True,
                        help="output_file_name")
    parser.add_argument('-w', '--windows', type=int, nargs='+', default=[1000000],
                        help="List of window sizes (e.g., --windows 10000 20000)")
    parser.add_argument('-r', '--resolution', type=int, default=5000,
                        help="cool 文件的resolution")
    parser.add_argument('-m', '--method', type=str, default="otsu",
                        help="用于计算边界阈值的方法，otsu更严格，li更宽松")
    parser.add_argument('-bw', '--bw_dir', type=str, required=True,
                        help="Path to the bigWig file")

    args = parser.parse_args()

    # 使用 args 传参调用主函数
    plot_boundary_signal(
        insulation_table_path=args.insulation_file,
        resolution=args.resolution,
        windows=args.windows,
        bw_path=args.bw_dir,
        output_file_name=args.output_file_name,
        method=args.method
    )


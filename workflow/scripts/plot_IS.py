#conda activate cooltools

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd

import cooler
import cooltools.lib.plotting
from cooltools import insulation
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.4'):
    raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
                         "please check your cooltools version and update to the latest")


def parse_args():
    parser = argparse.ArgumentParser(description="plot IS")

    # 添加命令行参数
    parser.add_argument('-c', '--cooler_file', type=str, required=True, help="Path to the cooler file")
    parser.add_argument('-is', '--insulation_file', type=str, required=True,
                        help="Path to the insulation scores file (tsv)")
    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help="Path to the output PNG file (e.g., output/plot.png)")
    parser.add_argument('-w','--windows', type=int, nargs='+', default=[1000000],
                        help="Comma-separated list of window range values (e.g., --windows 10000 20000 30000)")
    parser.add_argument('-s','--start', type=int, default=10_500_000,
                        help="start region")
    parser.add_argument('-l', '--length', type=int, default=90,
                        help="填windows的倍数")



    return parser.parse_args()

# 定义绘图辅助函数
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

def format_ticks(ax, x=True, y=True, rotate=True):
    bp_formatter = EngFormatter('b')
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

def main():
    args = parse_args()
    # load cool data
    clr = cooler.Cooler(args.cooler_file)  # 分辨率是文件里指定的
    resolution = clr.binsize  # 直接从文件提取

    # IS calculate
    # windows = [3 * resolution, 5 * resolution, 10 * resolution, 25 * resolution]
    # insulation_table = insulation(clr, windows, verbose=True)

    # 也可以直接load
    insulation_table = pd.read_csv(args.insulation_file, sep='\t')
    windows = [1000000]  # insulation 的窗口大小，匹配你的文件

    # 确定绘图区间
    start = args.start
    end = start + args.length * windows[0]
    region = ('2', start, end)


    #### 绘图 ####
    plt.switch_backend('Agg') # 设置 Matplotlib 使用非交互式后端
    # 绘制 hic 热图
    # 获取绘图区间数据
    data = clr.matrix(balance=True).fetch(region)
    norm = LogNorm(vmax=0.1, vmin=0.001) # 标准化参数
    # 绘制热图
    f, ax = plt.subplots(figsize=(18, 6))
    im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
    ax.set_aspect(0.5)
    ax.set_ylim(0, 10 * windows[0])
    format_ticks(ax, rotate=False)
    ax.xaxis.set_visible(False)
    # 添加图例-颜色条
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
    plt.colorbar(im, cax=cax)

    # 绘制下半部分的IS折线图
    # 选取 region 内的 insulation 分数
    ins_region = bioframe.select(insulation_table, region)
    ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
    ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0, 1, 5))))
    ins_ax.plot(
        ins_region[['start', 'end']].mean(axis=1).values,
        ins_region['log2_insulation_score_' + str(windows[0])].values,
        label=f'Window {windows[0]} bp'
    )
    # 添加图例
    ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)
    format_ticks(ins_ax, y=False, rotate=False)
    ax.set_xlim(region[1], region[2])

    # 保存
    plt.tight_layout()
    plt.savefig(args.output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {args.output_file}")

if __name__ == "__main__":
    main()


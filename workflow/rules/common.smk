import os

# 解析input，生成output
class SampleProcessor:
    def __init__(self, config):
        """
        初始化 SampleProcessor。

        参数:
            config (dict): 包含 reads 和 bamCoverage 的配置信息。
        """
        self.inputdir = config.get("input_dir", [])
        self.resolution = config.get("resolution",[])
        self.sample_names = []
        self.get_sample_names()

    def get_sample_names(self):
        """
        从给定的输入目录（inputdir）生成 sample_names，提取其下所有二级目录的名称
        """
        # 如果没有提供 inputdir 或者输入的 inputdir 无效，sample_names 设为 NA
        if self.inputdir is None or not os.path.isdir(self.inputdir):
            self.sample_names = ["NA"]
        else:
            # 获取 inputdir 下的所有子目录（一级目录）
            self.sample_names = [d for d in os.listdir(self.inputdir) if os.path.isdir(os.path.join(self.inputdir, d))]


    def generate_targets(self, sample):
        """
        生成输出结果
        """
        return [
            #f"Results/{sample}/bowtie_results/bwt2/{sample}/{sample}_mm10.bwt2pairs.bam",
            # f"Results/{sample}/hic_results/data/{sample}/{sample}_mm10.bwt2pairs.validPairs",
            f"Results/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
            f"Results/summary.txt",
            f"Results/{sample}/hic_results/pic/{sample}/plotHiCContactRanges_{sample}.pdf",
            f"Results/statistic/cumulative_plot.png",
            #f"Results/00collect/{sample}_{self.resolution}.insulation_scores.tsv"
        ]

    def get_all_targets(self):
        """
        获取所有样本的目标文件路径。

        返回:
            list[str]: 所有样本的目标路径列表
        """
        all_paths = []
        for sample in self.sample_names:
            paths = self.generate_targets(sample)
            all_paths.extend(paths)
        return all_paths

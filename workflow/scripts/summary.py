import re
from pathlib import Path
import pandas as pd
import numpy as np

# dir_path = Path("/home_data/home/slst/leixy2023/data/project/ref_data/hic_PGC_PRJNA484983/Results")  # Passed via params
# samples = ["E11_rep1","E11_rep2", "E13_F_rep1", "E13_F_rep2", "E13_M_rep1", "E13_M_rep2"]
# out_file = dir_path / 'summary.txt'

dir_path = Path(snakemake.params.results_dir)
samples = snakemake.params.sample_list
out_file = Path(snakemake.output.summary)

def extract_map(file_path, read_type):
    if read_type == 'R1':
        patterns = ['total_R1', 'mapped_R1']
    elif read_type == 'R2':
        patterns = ['total_R2', 'mapped_R2']
    else:
        raise ValueError("read_type must be 'R1' or 'R2'")

    data = {}
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            for key in patterns:
                match = re.search(fr'{key}\s+(\d+)', content)
                data[key] = int(match.group(1)) if match else 0
    except Exception as e:
        print(f"Error reading or parsing {file_path}: {e}")
    return data

def extract_total_pairs(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            return int(re.search(r'Total_pairs_processed\s+(\d+)', content).group(1))
    except Exception as e:
        print(f"Error reading or parsing {file_path}: {e}")
        return 0


def extract_valid_pairs(file_path):
    result = {}
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            result = {
                'valid_pairs': int(re.search(r'valid_interaction\s+(\d+)', content).group(1)),
                'valid_pairs_rmdup': int(re.search(r'valid_interaction_rmdup\s+(\d+)', content).group(1)),
                'cis_interaction': int(re.search(r'cis_interaction\s+(\d+)', content).group(1)),
                'cis_longRange': int(re.search(r'cis_longRange\s+(\d+)', content).group(1)),
            }
    except Exception as e:
        print(f"Error reading or parsing {file_path}: {e}")
        result = {k: 0 for k in ['valid_pairs', 'valid_pairs_rmdup', 'cis_interaction', 'cis_longRange']}
    return result

def extrac_inter_intra(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, usecols=[1, 2, 4, 5])
    df.columns = ['chr1', 'pos1', 'chr2', 'pos2']

    # 统计inter和intra互作
    same_chr_mask = df['chr1'] == df['chr2']
    diff_chr_count = (~same_chr_mask).sum()  # Counts rows where chr1 != chr2
    #same_chr_count = same_chr_mask.sum()  # Counts rows where chr1 == chr2

    # 对于染色体内互作，即cis，计算>20kb的数量
    intra_long_count = np.sum(np.abs(df.loc[same_chr_mask, 'pos1'] - df.loc[same_chr_mask, 'pos2']) > 20_000)

    return diff_chr_count, intra_long_count


def main():

    all_data = []

    for sample in samples:
        stats_path = dir_path / sample / "hic_results/stats" / sample
        map_R1 = stats_path / f"{sample}_R1.mmapstat"
        map_R2 = stats_path / f"{sample}_R2.mmapstat"
        total_pairs_path = stats_path / f"{sample}.mpairstat"
        valid_pairs_path = stats_path / f"{sample}_allValidPairs.mergestat"

        map_data_R1 = extract_map(map_R1, "R1")
        map_data_R2 = extract_map(map_R2, "R2")
        total_pairs = extract_total_pairs(total_pairs_path)
        valid_data = extract_valid_pairs(valid_pairs_path)

        total_reads = map_data_R1.get('total_R1', 0) + map_data_R2.get('total_R2', 0)
        mapped_reads = map_data_R1.get('mapped_R1', 0) + map_data_R2.get('mapped_R2', 0)
        mapping_rate = mapped_reads / total_reads if total_reads else 0

        valid_pairs = valid_data['valid_pairs']
        valid_pairs_rmdup = valid_data['valid_pairs_rmdup']
        duplicates = (valid_pairs - valid_pairs_rmdup) / valid_pairs if valid_pairs else 0
        Cis_long_Rate = valid_data['cis_longRange'] / valid_pairs_rmdup if valid_pairs_rmdup else 0

        allValidPairs_path = dir_path / "00collect" / f"{sample}.allValidPairs"
        inter,intra_long = extrac_inter_intra(allValidPairs_path)

        all_data.append({
            "sample": sample,
            "raw_reads": total_reads,
            "mapped_reads": mapped_reads,
            "mapping_rate": mapping_rate,
            "Total_Pairs": total_pairs,
            "Valid_Pairs": valid_pairs,
            "Duplicates": duplicates,
            "Non-redundant_valid_pairs": valid_pairs_rmdup,
            "Inter-chromosomal": inter,
            "Intra-chromosomal > 20kb": intra_long
        })

    df = pd.DataFrame(all_data)

    # 保留三位小数
    float_cols = ["Duplicates", "Cis_long_Rate"]
    for col in float_cols:
        if col in df.columns:
            df[col] = df[col].astype(float).round(4)

    df.to_csv(out_file, sep='\t', index=False)
    print(f"Saved summary to {out_file}")


if __name__ == "__main__":
    main()

import re
from pathlib import Path
import pandas as pd

# dir_path = Path("/home_data/home/slst/leixy2023/data/project/ref_data/hic_PGC_PRJNA484983/Results")  # Passed via params
# samples = ["E11_rep1","E11_rep2", "E13_F_rep1", "E13_F_rep2", "E13_M_rep1", "E13_M_rep2"]
# out_file = dir_path / 'summary.txt'

dir_path = Path(snakemake.params.results_dir)
samples = snakemake.params.sample_list
out_file = Path(snakemake.output.summary)

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


def main():

    all_data = []

    for sample in samples:
        sample_path = dir_path / sample / "hic_results/stats" / sample
        total_pairs_path = sample_path / f"{sample}.mpairstat"
        valid_pairs_path = sample_path / f"{sample}_allValidPairs.mergestat"

        total_pairs = extract_total_pairs(total_pairs_path)
        valid_data = extract_valid_pairs(valid_pairs_path)

        valid_pairs = valid_data['valid_pairs']
        valid_pairs_rmdup = valid_data['valid_pairs_rmdup']
        duplicates = (valid_pairs - valid_pairs_rmdup) / valid_pairs if valid_pairs else 0
        Cis_long_Rate = valid_data['cis_longRange'] / valid_pairs_rmdup if valid_pairs_rmdup else 0

        all_data.append({
            "sample": sample,
            "Total_Pairs": total_pairs,
            "Valid_Pairs": valid_pairs,
            "Duplicates": duplicates,
            "Non-redundant_valid_pairs": valid_pairs_rmdup,
            "Cis_Interaction": valid_data['cis_interaction'],
            "Cis_LongRange": valid_data['cis_longRange'],
            "Cis_long_Rate": Cis_long_Rate
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

#PBS -N multi_hic_IS.pbs
#PBS -l nodes=1:ppn=20
#PBS -S /bin/bash
#PBS -j oe
#PBS -q pub_fat
#PBS -l walltime=60:00:00

source ~/.bashrc
conda activate hicpro

# ==== Configuration ====
hic2juice=~/software/hicpro/bin/utils/hicpro2juicebox.sh
genomesize=/home_data/home/slst/leixy2023/data/database/mm10/archive/mm10.chrom.sizes
juicertool=~/software/juicer/scripts/common/juicer_tools.jar
enzyme_bed=/public/slst/home/zhenghui/Denghuanhuan/HiC_analysis/mm10_MboI.bed

inputdir=/home_data/home/slst/leixy2023/data/project/ref_data/hic_PGC_PRJNA484983/Results/valid_pairs_collect
outputdir=./

cd "$outputdir" || exit 1

samples=(E11 E13_M E13_F)

# ==== Main Loop ====
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"

    valid_pairs_file="${inputdir}/${sample}_mm10.bwt2pairs.validPairs"

    # Step 1: Generate .hic
    $hic2juice -i "$valid_pairs_file" -g "$genomesize" -j "$juicertool" -r "$enzyme_bed" -o "$outputdir"

    # Step 2: Convert to .cool
    hic_file="${sample}_mm10.bwt2pairs.validPairs.hic"
    cool_file="${sample}_40000.cool"
    hicConvertFormat -m "$hic_file" --inputFormat hic --outputFormat cool -r 40000 -o "$cool_file"

    # Step 3: Balance matrix
    cooler balance "$cool_file"

    # Step 4: Calculate insulation score
    cooltools insulation -p 4 -o "${sample}.100k.insulation_scores.tsv" --chunksize 20000000 "$cool_file" 1000000
done

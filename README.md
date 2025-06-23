# snakemake_chip
对log文件夹加上保护，避免在删除整个文件夹时被删除

配置环境除了env里通过conda安装的，还需要另外安装hicpro和juicer

创建符合输入格式的软连接
```shell
bash workflow/scripts/create_input_link.sh .test/data rawdata/input_ln
```

生成config-hicpro.txt
```shell
snakemake generate_hicpro_config
```

运行hicpro
```shell
# 本地运行
snakemake -np --use-conda --config input_dir="rawdata/input_ln" --jobs 1
# 提交pbs
nohup snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "python workflow/scripts/submit_job.py --config config/cluster_config.yaml --seqtype "hicpro" --sample {wildcards} --rule {rule}"  \
    --latency-wait 60 \
    --jobs 1 \
    --use-conda \
    --config input_dir="rawdata/input_ln" \
    > snakemake_hicpro_rep2.log 2>&1 &
    
nohup snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "python workflow/scripts/submit_job.py --config config/cluster_config.yaml --seqtype "hicpro" --sample {wildcards} --rule {rule}"  \
    --latency-wait 60 \
    --jobs 2 \
    --use-conda \
    --config input_dir="rawdata/inputdir" \
    > snakemake_hicpro_mESC.log 2>&1 &
    
    
    
nohup snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "python workflow/scripts/submit_job.py --config config/cluster_config.yaml --seqtype "hicpro" --sample {wildcards} --rule {rule}"  \
    --latency-wait 60 \
    --jobs 1 \
    --use-conda \
    --config input_dir="rawdata/input_ln_rep1" \
    > snakemake_hicpro_E11_250526.log 2>&1 &
    

```


```shell
HiC-Pro -c config/config-hicpro.txt -i Results/hic_results/data -o Results -s merge_persample
```
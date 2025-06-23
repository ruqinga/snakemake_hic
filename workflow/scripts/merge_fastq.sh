#!/bin/bash

echo "---Starting the merging process---"

# 检查命令行参数是否正确
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_folder> <output_folder>"
    exit 1
fi

input_folder="$1"
output_folder="$2"

# 检查输入文件夹是否存在
if [ ! -d "$input_folder" ]; then
    echo "Input folder does not exist: $input_folder"
    exit 1
fi

# 创建输出文件夹（如果不存在）
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

# 创建一个关联数组来存储文件前缀和对应的文件列表
declare -A file_map_1
declare -A file_map_2

# 使用 mapfile 将 find 命令的输出读入数组
mapfile -t files < <(find "$input_folder" -name "*.fastq.gz")

# 遍历所有文件
for file in "${files[@]}"; do
    # 提取文件名的前缀（去除 _batch[num].fastq.gz 后缀）
    filename=$(basename "$file")
    
    # 检查文件名是否包含 _batch 后跟数字
    if [[ "$filename" =~ _batch[0-9]+ ]]; then
        # 如果包含，则提取 prefix
        prefix="${filename%%_batch[0-9]*}"
    else
        # 如果不包含，则输出不处理信息
        echo "不处理$file"
        continue
    fi

    # 根据文件名后缀区分两端 (_1 和 _2)
    if [[ "$filename" =~ _1.fastq.gz ]]; then
        # 如果是 _1.fastq.gz 文件，添加到 file_map_1
        file_map_1["$prefix"]+="$file "
    elif [[ "$filename" =~ _2.fastq.gz ]]; then
        # 如果是 _2.fastq.gz 文件，添加到 file_map_2
        file_map_2["$prefix"]+="$file "
    fi
done

# 遍历每个前缀，分别合并 _1 和 _2 文件
for prefix in "${!file_map_1[@]}"; do
    # 对 _1 文件进行合并
    merge_files_1=$(echo "${file_map_1[$prefix]}" | xargs)
    echo -e "\nMerging the following _1 files into $output_folder/${prefix}_1.fastq.gz:"
    for merge_file in $merge_files_1; do
        echo "  $merge_file"
    done
    # 使用 zcat 合并 _1 文件并重定向到最终文件
    zcat $merge_files_1 | pigz -p 20 > "$output_folder/${prefix}_1.fastq.gz"

    # 对 _2 文件进行合并
    merge_files_2=$(echo "${file_map_2[$prefix]}" | xargs)
    echo -e "\nMerging the following _2 files into $output_folder/${prefix}_2.fastq.gz:"
    for merge_file in $merge_files_2; do
        echo "  $merge_file"
    done
    # 使用 zcat 合并 _2 文件并重定向到最终文件
    zcat $merge_files_2 | pigz -p 20 > "$output_folder/${prefix}_2.fastq.gz"
done

echo -e "\nAll files have been merged and compressed."



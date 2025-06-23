#!/bin/bash


## 在有子文件夹的情况下将input dir下的文件都转换为软链接
## 检查参数
#if [ -z "$1" ]; then
#    echo "用法: $0 <inputdir> [outputdir]"
#    exit 1
#fi
#
#inputdir="$1"
#outputdir="${2:-{inputdir}/rawdata_ln}"
#
#mkdir -p "$outputdir"
#
#for sample in "$inputdir"/*; do
#    [ -d "$sample" ] || continue  # 如果不是目录就跳过
#    sample_name=$(basename "$sample")
#    inner_dir="$outputdir/$sample_name/$sample_name"
#    mkdir -p "$inner_dir"
#
#    for fq in "$sample"/*.fastq.gz; do
#        [ -f "$fq" ] || continue  # 如果没有匹配文件就跳过
#        ln -s "$(realpath "$fq")" "$inner_dir/"
#    done
#done


# 在没有子文件夹的情况下转换

#!/bin/bash

# 检查是否提供了输入目录
if [ -z "$1" ]; then
    echo "用法: $0 <inputdir> [outputdir]"
    exit 1
fi

# 设置输入输出目录
inputdir="$1"
outputdir="${2:-${inputdir}/rawdata_ln}"

# 创建输出目录（如果不存在的话）
if [ ! -d "$outputdir" ]; then
    mkdir -p "$outputdir"
    echo "创建文件夹：$outputdir"
fi

# 遍历输入目录下所有 .fastq.gz 文件
for file in "$inputdir"/*.fastq.gz; do
    # 提取文件名前缀（去除 _R1.fastq.gz 或 _R2.fastq.gz 部分）
    filename=$(basename "$file")
    # 检查文件名是否以 _R1.fastq.gz 或 _R2.fastq.gz 结尾
    if [[ "$filename" =~ _R1.fastq.gz$ ]]; then
        prefix="${filename%_R1.fastq.gz}"
    elif [[ "$filename" =~ _R2.fastq.gz$ ]]; then
        prefix="${filename%_R2.fastq.gz}"
    else
        echo "警告：文件名 $filename 不符合预期格式，跳过处理。"
        continue
    fi

    # 创建以前缀命名的文件夹
    prefix_folder="$outputdir/$prefix/$prefix"
    if [ ! -d "$prefix_folder" ]; then
        mkdir -p "$prefix_folder"
        #echo "创建文件夹：$prefix_folder"
    fi

    # 创建软连接到文件夹中
    ln -s "$(realpath "$file")" "$prefix_folder/$filename"
    #echo "创建软连接：$prefix_folder/$filename -> $file"
done

echo "软连接创建完成"
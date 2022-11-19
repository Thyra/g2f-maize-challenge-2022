#!/bin/bash

file_loc=$1
file_name=$2
miss_thr=$3 # float : Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
mac_thr=$4 # integer : Include only sites with Minor Allele Count greater than or equal to the "--mac" value 

# put everthing processed in /proj/g2f-maize-challenge-2022/processed_data
pro_dir_path="/qg-10/data/AGR-QG/Gogna/g2f-maize-challenge-2022/processed_data"

if [[ ! -d "$pro_dir_path" ]]
then
  mkdir "$pro_dir_path"
fi

# create a easy to work with format of vcf

if [[ ! -f  "$pro_dir_path/$file_name.vcf.gz"  ]]
then
  unzip "$file_loc/$file_name.vcf.zip" -d "$pro_dir_path"
  bgzip "$pro_dir_path/$file_name.vcf"
  tabix -p vcf "$pro_dir_path/$file_name.vcf.gz" 
else
  echo "vcf.gz exists"
fi

# create site specific stats
if [[ ! -f  "$pro_dir_path/$file_name.lmiss"  ]]
then
  vcftools --gzvcf "$pro_dir_path/$file_name.vcf.gz" \
    --missing-site \
    --stdout 2> "$pro_dir_path/$file_name.lmiss.log" \
    > "$pro_dir_path/$file_name.lmiss"
else
  echo "lmiss file exists"
fi

# create sample specific stats
if [[ ! -f  "$pro_dir_path/$file_name.imiss"  ]]
then
  vcftools --gzvcf "$pro_dir_path/$file_name.vcf.gz" \
    --missing-indv \
    --stdout 2> "$pro_dir_path/$file_name.imiss.log" \
    > "$pro_dir_path/$file_name.imiss"
else
  echo "imiss file exists"
fi

if [[ ! -f  "$pro_dir_path/$file_name.het"  ]]
then
  vcftools --gzvcf "$pro_dir_path/$file_name.vcf.gz" \
    --het \
    --stdout 2> "$pro_dir_path/$file_name.het.log" \
    > "$pro_dir_path/$file_name.het"
else
  echo "het file exists"
fi


# filter the vcf for missing values, monomorphic markers and minor allele count
if [[ ! -f  "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.vcf.gz"  ]]
then
  vcftools --gzvcf "$pro_dir_path/$file_name.vcf.gz" \
    --max-missing $miss_thr \
    --mac $mac_thr \
    --recode --recode-INFO-all \
    --stdout 2> "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.vcf.gz.err" \
    | bgzip -c > "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.vcf.gz"
fi

# recode to numeric
if [[ ! -f  "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.012"  ]]
then
  vcftools --gzvcf "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.vcf.gz" \
    --012 \
    --recode --recode-INFO-all \
    --stdout 2> "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.012.err" \
    | bgzip -c > "$pro_dir_path/$file_name.miss.$miss_thr.mac.$mac_thr.012"
fi
# g2f-maize-challenge-2022
https://eval.ai/web/challenges/challenge-page/1878/overview

# Processing steps

## Genodata
It was first processed in bash with scripts/process_geno_data.sh and then minimally formatted to a usable csv with numeric genotype calls. This is stored at processed_data/geno_processed.miss.1.mac.1.biallelic.txt.xz 

## Others
training data was processed with scripts/preprocessing.ipynb and put as combined imputed data at processed_data/combined_mat.csv.xz. The same script also has code used to create processed_data/train_test_split_v2.json.xz

# Usage 
Since the object sizes are becoming large. I compressed all processed files with xz (https://linux.die.net/man/1/xz). You would need to download and ucompress it for further use. 
the processed_data/combined_mat.csv.xz has no geno data but you can connect the same (processed_data/geno_processed.miss.1.mac.1.biallelic.txt.xz) with it using Hybrid <-> ge_dta_Hybrid key <-> value pair. 

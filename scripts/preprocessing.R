library(tidyverse)
library(qs)
library(hablar)
library(asreml)
options("scipen"=10, "digits"=2)

# paths -------------------------------------------------------------------

tmp_at <- '/proj/g2f-maize-challenge-2022/tmp_data'
dump_at <- '/proj/g2f-maize-challenge-2022/dump'
source_data <- '/proj/g2f-maize-challenge-2022/source_data/Training_Data'
processed_data <- '/proj/g2f-maize-challenge-2022/processed_data'

# Phenodata -----------------------------------------------------------
## overview

trait_data <- read.csv(sprintf("%s/1_Training_Trait_Data_2014_2021.csv", source_data), sep = ",")
#Env - Environment (combination of location of evaluation and year) that will match fields in environmental data sets
#Year	- Year of evaluation.
#Field-Location	- G2F Field-location name.
#Experiment	- Experiment name.
#Replicate	- Large-scale field block.
#Block	- Smaller-scale field block nested within Replicate.     
#Plot	- Designation of individual experimental unit.
#Range	- Designation of field range of the plot (ranges are organized perpendicular to corn rows).
#Pass	- Designation of field pass of the plot (passes are organized parallel to corn rows. Combination of range and pass form coordinate grid system describing location of each plot within the field).
#Hybrid -	Hybrid name that will match genotype file if it has been genotyped.
#Hybrid_orig_name -	Hybrid name as in the DOI.
#Hybrid_Parent1 -	Hybrid’s parental inbred line 1.
#Hybrid_Parent2	- Hybrid ‘s parental inbred line 2.
#Plot_Area_ha	- Calculated plot area in hectares.
#Date_Planted -	Date the plot was planted. (Month/Day/Year)
#Date_Harvested	- Date the plot was harvested. (Month/Day/Year)
#Stand_Count_plants -	Number of plants, per plot, at harvest.
#Pollen_DAP_days -	Number of days after planting that 50% of plants in the plot began shedding pollen.
#Silk_DAP_days	- Number of days after planting that 50% of plants in the plot had visible silks.
#Plant_Height_cm	- Measured as the distance between the base of a plant and the ligule of the flag leaf (centimeter).
#Ear_Height_cm	- Measured as the distance from the ground to the primary ear bearing node (centimeter).
#Root_Lodging_plants	- Number of plants per plot that show root lodging.
#Stalk_Lodging_plants	- Number of plants per plot broken between ground level and top ear node at harvest.
#Yield_Mg_ha	- Grain yield in Mg per ha at 15.5% grain moisture, using plot area without alley (Mg/ha).
#Grain_Moisture	- Water content in grain at harvest (percentage).
#Twt_kg_m3	- Shelled grain test weight (kg/m3), a measure of grain density.

missing_data_overview <- trait_data %>% group_by(Env) %>% 
  summarize(vals = n(),
            missing = sum(is.na(Yield_Mg_ha)),
            missing_prop = missing/vals) %>%
  arrange(desc(missing_prop)) # shows proportion of missing data points for grain yield
-------------------------------------------------
# Model selection and fit
# random factors BLUEs
BLUE_rf1 <- sprintf("~ %s", "Replicate")
BLUE_rf2 <- sprintf("~ %s", "Block")
BLUE_rf3 <- sprintf("~ %s + %s:%s", "Replicate", "Block", "Replicate")
BLUE_rf4 <- sprintf("~ %s + %s:%s + %s", "Replicate", "Block", "Replicate", "units")

# random factors BLUPs
BLUP_rf1 <- sprintf("~ %s + %s", "Hybrid", "Replicate")
BLUP_rf2 <- sprintf("~ %s + %s", "Hybrid", "Block")
BLUP_rf3 <- sprintf("~ %s + %s + %s:%s",  "Hybrid", "Replicate", "Block", "Replicate")
BLUP_rf4 <- sprintf("~ %s + %s + %s:%s + %s",  "Hybrid", "Replicate", "Block", "Replicate", "units")

# residual structure
r1 <- sprintf("~ ar1v(%s):ar1(%s)", "Range", "Pass")
r2 <- sprintf("~ ar1(%s):idv(%s)", "Range", "Pass")

#asreml options
asreml.options(maxit = 50,             
               workspace = "512mb",       
               pworkspace = "512mb",
               trace = F,
               extra = 10)

n_env <- trait_data %>% distinct(Env) %>% pull(Env) # 217 environments

non_mono_miss <- function(x){
  out <- ifelse(sum(is.na(x)) == length(x)| length(unique(x)) == 1, NA, length(unique(x)))
  return(out)
}

any_miss <- function(x){
  out <- ifelse(sum(is.na(x)) >1, TRUE, FALSE)
  return(out)
}

check_uq_count <- function(x, y){
  out <- ifelse(length(count(as.character(x), as.character(y)) %>% filter(n > 1)), FALSE, TRUE)
  return(out)
}

overview_trial <- trait_data %>% group_by(Env) %>% 
  mutate(range_pass = paste0(Range, "_", Pass)) %>%
  summarize(reps = non_mono_miss(Replicate),
            miss_reps = ifelse(!is.na(reps), any_miss(Replicate), NA),
            blocks = non_mono_miss(Block),
            miss_blocks = ifelse(!is.na(blocks), any_miss(Block), NA),
            ranges = non_mono_miss(Range),
            miss_ranges = ifelse(!is.na(ranges), any_miss(Range), NA),
            passes = non_mono_miss(Pass),
            miss_passes = ifelse(!is.na(passes), any_miss(Pass), NA),
            rep_confounds_block = ifelse(sum(Replicate == Block) == n(), TRUE, FALSE),
            range_pass_uq = ifelse(length(which(table(range_pass) > 1)) > 0, FALSE, TRUE),
            range_pass_uq_resi = ifelse(length(unique(Range))*length(unique(Pass)) != n(), FALSE, TRUE),
            .groups = "drop") %>% arrange(rep_confounds_block, ranges, passes) %>%
  mutate(BLUP_fixed = "Yield_Mg_ha ~ 1",
         BLUE_fixed = "Yield_Mg_ha ~ Hybrid",
         id = row_number())

env_info_replicate <- overview_trial %>% filter(rep_confounds_block == TRUE, is.na(ranges)) %>%
  bind_rows(overview_trial %>% filter(Env == "TXH1-Early_2018")) %>%  # the exception
  mutate(BLUP_random = ifelse(!is.na(reps), BLUP_rf1, BLUP_rf2),
         BLUP_residual = "~ units",
         BLUE_random = ifelse(!is.na(reps), BLUE_rf1, BLUE_rf2),
         BLUE_residual = "~ units")

env_info_replicate_range_pass <- overview_trial %>% filter(rep_confounds_block == TRUE, !is.na(ranges)) %>%
  mutate(BLUP_random = BLUP_rf1,
         BLUP_residual = ifelse(range_pass_uq & range_pass_uq_resi, r1, "~ units"),
         BLUE_random = BLUE_rf1,
         BLUE_residual =  ifelse(range_pass_uq & range_pass_uq_resi, r1, "~ units"))

env_info_replicate_block <- overview_trial %>% filter(rep_confounds_block == FALSE, !is.na(reps)) %>%
  bind_rows(overview_trial %>% filter(rep_confounds_block == FALSE, is.na(reps), !is.na(ranges)))

env_info_replicate_block_set_1 <- env_info_replicate_block %>% 
  filter(!is.na(reps), !is.na(blocks), !is.na(ranges), !is.na(passes)) %>%
  filter(Env != "OHH1_2019") %>% 
  mutate(BLUP_random = BLUP_rf3,
         BLUP_residual = ifelse(range_pass_uq & range_pass_uq_resi & miss_ranges == FALSE, r1, "~ units"),
         BLUE_random = BLUE_rf3,
         BLUE_residual =  ifelse(range_pass_uq & range_pass_uq_resi & miss_ranges == FALSE, r1, "~ units"))

env_info_replicate_block_set_2 <- env_info_replicate_block %>% 
  anti_join(env_info_replicate_block_set_1, by = "id") %>%
  filter(is.na(reps) | is.na(blocks), !is.na(ranges), !is.na(passes), miss_ranges == FALSE, miss_passes == FALSE) %>% 
  mutate(BLUP_random = ifelse(!is.na(reps), BLUP_rf1, BLUP_rf2),
         BLUP_residual = ifelse(range_pass_uq & range_pass_uq_resi, r1, "~ units"),
         BLUE_random = ifelse(!is.na(reps), BLUP_rf1, BLUP_rf2),
         BLUE_residual =  ifelse(range_pass_uq & range_pass_uq_resi, r1, "~ units")) # here residuals may be mapped but block or rep is missing

env_info_replicate_block_set_3 <- env_info_replicate_block %>% 
  anti_join(env_info_replicate_block_set_1, by = "id") %>%
  anti_join(env_info_replicate_block_set_2, by = "id") %>% 
  mutate(BLUP_random = NA,
         BLUP_residual = NA,
         BLUE_random = NA,
         BLUE_residual =  NA) 

model_overview <- env_info_replicate %>%
  bind_rows(env_info_replicate_range_pass) %>%
  bind_rows(env_info_replicate_block_set_1) %>%
  bind_rows(env_info_replicate_block_set_2) %>%
  bind_rows(env_info_replicate_block_set_3) %>%
  arrange(id)

# prepare data
trait_data_tibble <- trait_data %>%
  as_tibble() %>%
  convert(fct(Hybrid,  Replicate, Block, Range, Pass)) %>%
  arrange(Replicate, Block, Range, Pass)

an.error.occured <- FALSE
tryCatch( { result <- log("not a number"); print(res) }
          , error = function(e) {an.error.occured <<- TRUE})
print(an.error.occured)


# get results
output <- list()
caught <- NULL
for(i in model_overview$Env[1:10]){
  # reset variables
  trait_data_subset <- NULL
  overview <- NULL
  multiple <- NULL
  data_out <- list()
  model_BLUP <- NULL
  model_BLUE <- NULL
  sigma.g <- NULL
  sigma.e <- NULL
  Geno <- NULL
  error <- FALSE
  
  # get specific data
  trait_data_subset <- trait_data_tibble %>% filter(Env == i) %>% arrange(Replicate, Block, Range, Pass)
  overview <- model_overview %>% filter(Env == i) %>% as.data.frame()

  tryCatch( {
    model_BLUP <- asreml(fixed = as.formula(overview$BLUP_fixed),
                         random = as.formula(overview$BLUP_random),
                         residual = as.formula(overview$BLUP_residual),
                         data = trait_data_subset)
    model_BLUE <- asreml(fixed = as.formula(overview$BLUE_fixed),
                         random = as.formula(overview$BLUE_random),
                         residual = as.formula(overview$BLUE_residual),
                         data = trait_data_subset)
    catch <- i
    multiple <- length(unique(trait_data_subset$Replicate))
    
    if(is.null(catch)){caught <- catch}else{caught <- c(caught, catch)}
    sum_mat <- summary(model_BLUP)$varcomp
    
    # heritability
    sigma.g <- sum_mat["Hybrid","component"]
    if(grepl("Range|Pass", overview$BLUP_residual)){
      sigma.e <- mean(sum_mat[grep("Range|Pass", rownames(sum_mat), value = T), "component"])
    } else if (grepl("units", overview$BLUP_residual)) {
      sigma.e <- sum_mat["units!R","component"]
    } else {
      stop()
    }
    
    # BLUEs
    Geno <- predict.asreml(model_BLUE, classify = "Hybrid",)$pvals[,1:2]
    colnames(Geno)[2] <- i
    
    # export results
    data_out[["BLUEs"]] <- as.data.frame(Geno, row.names = NULL)
    data_out[["rep"]] <- sigma.g/(sigma.g+sigma.e/multiple)
  }, error = function(e) {
    error <<- TRUE
  })
  
  if(error) {
    data_out[["BLUEs"]] <- NA
    data_out[["rep"]] <- i
    print(sprintf("%s of %s done. there was an error.", overview$id, nrow(model_overview)))
  }
  output[[i]] <- data_out
  
  # print info
  if((overview$id %% 25) ==0){
    print(sprintf("%s of %s done-------------------------------------", overview$id, nrow(model_overview)))
  }
}

rep <- do.call(rbind, lapply(output, function(x) x[["rep"]]))

model_overview_with_rep <- data.frame("Env" = rownames(rep),
                                      "repetability" = rep) %>%
  left_join(model_overview, by = "Env") %>%
  select(-all_of(colnames(model_overview)[2:12])) %>%
  select(-id)

missing_rep <- model_overview_with_rep %>% filter(is.na(rep)) %>%
  pull(Env) %>% as.vector() # 29
missing_rep[which(missing_rep %in% caught)]

# prepare output
BLUEs <- as.data.frame(do.call(rbind, lapply(output, function(x) {
  if(is.data.frame(x[["BLUEs"]])){
    data <- x[["BLUEs"]]
    Env <- colnames(data)[2]
    data$Env <- Env
    colnames(data) <- c("Hybrid", "Yield_Mg_ha", "Env")
  } else {
    data <- data.frame("Hybrid" = NA,
                       "Yield_Mg_ha" = NA,
                       "Env" = x[["rep"]])
  }
  return(data)
  })))
row.names(BLUEs) <- NULL

# Genodata ----------------------------------------------------------------
system(sprintf('zcat %s/5_Genotype_Data_All_Years.vcf.gz | grep -v -P "#" | cut -f1-2 > %s/marker_names', processed_data, processed_data))
system(sprintf("zcat %s/5_Genotype_Data_All_Years.vcf.gz | head -n 30 | grep '#' | tail -n 1 > %s/sample_names", processed_data, processed_data))

marker_names <- read.table(sprintf("%s/marker_names", processed_data), header = F) # unfiltered - 437214
sample_names <- strsplit(readLines(sprintf("%s/sample_names", processed_data)), "\t")[[1]]
sample_names_df <- data.frame("Hybrid" = sample_names[10:length(sample_names)]) # 4928

# check overlap between genotypic and phenotypic data

overlap_check_geno_pheno <- trait_data %>% distinct(Hybrid) %>% 
  mutate(present = ifelse(Hybrid %in% sample_names_df$Hybrid, TRUE, FALSE))
overlap_check_geno_pheno %>% count(present)
# FALSE  260 # these hybrids has no genotypic data available. Most likely these genotypes need to be filtered from the genodata. Else we need to do a manual search
# TRUE 4423

miss_stats <- read.table(sprintf("%s/5_Genotype_Data_All_Years.lmiss", processed_data), header = T) %>%
  mutate(miss_prop = N_MISS/N_DATA) %>% filter(miss_prop == 0)

#Var1   Freq
# 0   49642  # the one to start working
#(0,0.1] 359131 # here onwards some sort of imputation will be needed. I will initially upload these markers. 
#(0.1,0.2]  18679 
#(0.2,0.3]   5166
#(0.3,0.4]   2385
#(0.4,0.5]   1142
#(0.5,0.6]    612
#(0.6,0.7]    261
#(0.7,0.8]    107
#(0.8,0.9]     75
#(0.9,1]     14

# Env data ----------------------------------------------------------------

env_data <- read.csv(sprintf("%s/4_Training_Weather_Data_2014_2021.csv", source_data), header = T) # 16 variables
overlap_check_env_pheno <- trait_data %>% distinct(Env) %>% 
  mutate(present = ifelse(Env %in% env_data$Env, TRUE, FALSE))
overlap_check_env_pheno %>% count(present)
#FALSE   5 # these env will be removed from phenodata
#TRUE 212

# Soil data ----------------------------------------------------------------

soil_data <- read.csv(sprintf("%s/3_Training_Soil_Data_2015_2021.csv", source_data), header = T) 
overlap_check_soil_pheno <- trait_data %>% distinct(Env) %>% 
  mutate(present = ifelse(Env %in% soil_data$Env, TRUE, FALSE))
overlap_check_soil_pheno %>% count(present)
# FALSE  76 # these have also to be removed 
# TRUE 141 # these will also have to be removed

# EC data  ----------------------------------------------------------------
EC_data <- read.csv(sprintf("%s/6_Training_EC_Data_2014_2021.csv", source_data), header = T) 
overlap_check_EC_pheno <- trait_data %>% distinct(Env) %>% 
  mutate(present = ifelse(Env %in% EC_data$Env, TRUE, FALSE))
overlap_check_EC_pheno %>% count(present)
#FALSE   52 # these env will be removed from phenodata
#TRUE 165

# Meta data ---------------------------------------------------------------
meta_data <- read.csv(sprintf("%s/2_Training_Meta_Data_2014_2021_utf_encoded.csv", source_data), header = T) 
overlap_check_meta_pheno <- trait_data %>% distinct(Env) %>% 
  mutate(present = ifelse(Env %in% meta_data$Env, TRUE, FALSE))
overlap_check_meta_pheno %>% count(present) # all present


# Get overview of availability --------------------------------------------

model_overview_availability <- model_overview %>% distinct(Env) %>%
  mutate(present_env = ifelse(Env %in% env_data$Env, TRUE, NA),
         present_soil = ifelse(Env %in% soil_data$Env, TRUE, NA),
         present_EC = ifelse(Env %in% EC_data$Env, TRUE, NA),
         present_meta = ifelse(Env %in% meta_data$Env, TRUE, NA)) %>%
  filter(!is.na(present_env), !is.na(present_soil), !is.na(present_EC), !is.na(present_meta)) # 109 env

# todo: make data available for 109 envs to git 

# Libraries and packages --------------------------------------------------
packs <- c("tidyverse", "qs", "hablar", "data.table", "AGHmatrix"
           , "jsonlite", "RSpectra", "foreach", "doParallel"
           , "Metrics", "BGLR")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))

if (!any(success)) {
  missing_packs <- names(success)[!success]
  install.packages(missing_packs)
  sapply(missing_packs, require, character.only = TRUE)
  cat(sprintf("%s\npackages were absent. These were installed and loaded", paste0(missing_packs, collapse = "\t")))
} else {
  print("All packages were present and are loaded!")
}

setDTthreads(10)
options("scipen"=10, "digits"=2)

yield_rmse <- function(run, cv_data, eig_g_at,
                       eig_e_at, eig_ge_at, pheno_data, 
                       save_loc, nIter=100, burnIn=10){
  run_name <- names(cv_data)[run]
  log_at <-  sprintf("%s/logs", save_loc)
  data_at <-  sprintf("%s/data", save_loc)
  result_at <- sprintf("%s/results", save_loc)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  if(!dir.exists(data_at)){dir.create(data_at, recursive = T)}
  if(!dir.exists(result_at)){dir.create(result_at, recursive = T)}
  
  # get decomposed matrices
  Kg_eigen <- qread(eig_g_at)
  Ke_eigen  <- qread(eig_e_at)
  Kge_eigen <- qread(eig_ge_at)
  
  # get pheno_data
  phenoGE <- pheno_data
  phenoGE$obs <- phenoGE$Yield_Mg_ha
  phenoGE$obs[unlist(cv_data[[run]][["test"]])] <- NA
  
  cat(sprintf("phenoGE defined for run %s", run), 
      file = sprintf("%s/run_%s.log", log_at, run),
      sep = "\n")
  
  ## set ETA and run model
  ETA <- list(G = list(X=Kg_eigen, model='BRR'),
              E = list(X=Ke_eigen, model='BRR'),
              GE = list(X=Kge_eigen, model='BRR'))
  
  cat(sprintf("ETA defined for run %s", run), 
      file = sprintf("%s/run_%s.log", log_at, run),
      sep = "\n",
      append = T)
  
  ## set BLAS threads
  #RhpcBLASctl::blas_set_num_threads(1)
  
  t0 <- Sys.time()
  model_fit <- BGLR(y=phenoGE$obs,
                    ETA=ETA,
                    nIter=nIter,
                    burnIn=burnIn,
                    verbose = FALSE,
                    saveAt=sprintf("%s/run_%s", data_at, run))
  t1 <- Sys.time()
  cat(sprintf("model fitting took %s minutes", difftime(t1, t0, units='mins')), 
      file = sprintf("%s/run_%s.log", log_at, run),
      append = T,
      sep = "\n")
  
  phenoGE$pred <- model_fit$yHat
  phenoGE$run_name <- run_name
  
  qsave(phenoGE, sprintf("%s/run_%s.qs", result_at, run))
  
  mean_rmse <- phenoGE %>% filter(is.na(obs)) %>% group_by(Env) %>% 
    summarize(RMSE = rmse(Yield_Mg_ha, pred), .groups = "drop") %>%
    pull(RMSE) %>% mean()
  mean_r2 <- phenoGE %>% filter(is.na(obs)) %>% group_by(Env) %>% 
    summarize(r2 = cor(Yield_Mg_ha, pred), .groups = "drop") %>%
    pull(r2) %>% mean()
  
  cat(sprintf("run_name = %s ; Mean RMSE (over environments) = %s", run_name, mean_rmse),
      file = sprintf("%s/run_%s.log", log_at, run),
      append = T,
      sep = "\n")
  
  # produce output
  out <- cbind("run" = run_name,
               "G" = model_fit$ETA$G$varB,
               "E"= model_fit$ETA$E$varB,
               "GxE"= model_fit$ETA$GE$varB,
               "err"= model_fit$varE,
               "Mean_RMSE" = mean_rmse,
               "Mean_cor" = mean_r2)
  
  return(out)
}

eigen_decomp <- function(inc_mat = NULL, mat, save_at){
  t1 <- Sys.time()
  mat_eigen <- eigen(mat)
  mat_eigen$vectors <- mat_eigen$vectors[, mat_eigen$values>1e-8]
  mat_eigen$values <- mat_eigen$values[mat_eigen$values>1e-8]
  mat_PC <- sweep(mat_eigen$vectors,MARGIN=2,STATS=sqrt(mat_eigen$values),FUN='*')
  t2 <- Sys.time()
  print(sprintf("decomp took %s minutes", difftime(t2, t1, units='mins')))
  if(!is.null(inc_mat)){
    final_mat <- inc_mat %*% mat_PC
  } else {
    final_mat <- mat_PC
  }
  qs::qsave(final_mat, save_at)
  print(sprintf("saved at %s", save_at))
  return(final_mat)
}

generate_eigen_decomps <- function(G_mat, E_mat, pheno_data, eig_g_at, eig_e_at, eig_ge_at){
  mat_names <- paste0(pheno_data$Env, ":", pheno_data$Hybrid)
  if (!file.exists(eig_g_at)){
    Zg <- model.matrix(~ -1 + Hybrid, pheno_data)
    colnames(Zg) <- gsub('Hybrid', '', colnames(Zg), perl = T)
    G_mat_ordered <- G_mat[colnames(Zg), colnames(Zg)]
    
    # decompose and save
    decomp_g <- eigen_decomp(inc_mat = Zg,
                             mat = G_mat_ordered,
                             save_at = eig_g_at)
    rm(Zg, G_mat_ordered, decomp_g)
  }
  
  if(!file.exists(eig_e_at)){
    Ze <- model.matrix(~ -1 + Env, pheno_data)
    colnames(Ze) <- gsub('Env', '', colnames(Ze), perl = T)
    E_mat_ordered <- E_mat[colnames(Ze), colnames(Ze)]
    
    # decompose and save
    decomp_e <- eigen_decomp(inc_mat = Ze,
                             mat = E_mat_ordered,
                             save_at = eig_e_at)
    rm(Ze, E_mat_ordered, decomp_e)
  }
  
  if(!file.exists(eig_ge_at)){
    Zg <- model.matrix(~ -1 + Hybrid, pheno_data)
    colnames(Zg) <- gsub('Hybrid', '', colnames(Zg), perl = T)
    G_mat_ordered <- G_mat[colnames(Zg), colnames(Zg)]
    Kg <- Zg %*% tcrossprod(G_mat_ordered, Zg)
    colnames(Kg) <- rownames(Kg) <- mat_names
    
    Ze <- model.matrix(~ -1 + Env, pheno_data)
    colnames(Ze) <- gsub('Env', '', colnames(Ze), perl = T)
    E_mat_ordered <- E_mat[colnames(Ze), colnames(Ze)]
    Ke <- Ze %*% tcrossprod(E_mat_ordered, Ze)
    colnames(Ke) <- rownames(Ke) <- mat_names
    
    Kge <- Kg*Ke
    
    # decompose and save
    Kge_eigen <- eigen_decomp(mat = Kge, 
                              save_at = eig_ge_at)
    rm(Zg, G_mat_ordered, Kg, Ze, E_mat_ordered, Ke, Kge, Kge_eigen)
  }
  
  output <- list()
  output[["eig_g_at"]] <- eig_g_at
  output[["eig_e_at"]] <- eig_e_at
  output[["eig_ge_at"]] <- eig_ge_at
  return(output)
}

#source("https://raw.githubusercontent.com/cran/CovCombR/master/R/CovComb.R")
#source("https://raw.githubusercontent.com/denizakdemir/CovCombR/master/R/Hmatfunc.R")
# paths -------------------------------------------------------------------

paths <- list(
  "tmp_at" = '/proj/g2f-maize-challenge-2022/tmp_data',
  "dump_at" = '/proj/g2f-maize-challenge-2022/dump',
  "source_data" = '/proj/g2f-maize-challenge-2022/source_data/Training_Data',
  "submission_data" = '/proj/g2f-maize-challenge-2022/source_data/Testing_Data',
  "processed_data" = '/proj/g2f-maize-challenge-2022/processed_data',
  "results" = '/proj/g2f-maize-challenge-2022/processed_data'
)

# Load data ---------------------------------------------------------------

## phenodata
pheno_data_at <- sprintf("%s/pheno_data.qs", paths[["processed_data"]])
if (!file.exists(pheno_data_at)){
  
  data <- setDF(fread(sprintf("%s/combined_mat_v3_BLUES_env.csv", paths[["processed_data"]])))
  all_cols <- colnames(data)
  
  sl <- all_cols[grep("sl_", all_cols)]
  wt <- all_cols[grep("wt_", all_cols)]
  ec <- all_cols[grep("ec_", all_cols)]
  
  pheno_data <- data[, all_cols[!(all_cols %in% c(sl, ec, wt))]]
  qsave(pheno_data, pheno_data_at)
} else {
  pheno_data <- qread(pheno_data_at)
}

## env_data
E_mat_at <- sprintf("%s/E_mat.qs", paths[["processed_data"]])
if (!file.exists(E_mat_at)){
  #data <- setDF(fread(sprintf("%s/combined_mat_v3.csv", paths[["processed_data"]])))
  #all_cols <- colnames(data)
  #ec <- all_cols[grep("ec_", all_cols)]
  
  # weather data derived kinship 
  #meta_data_sub <-  read.csv(sprintf("%s/2_Testing_Meta_Data_2022.csv", paths[["submission_data"]])) %>% filter(Date_Planted != "")
  weather_data_sub <- read.csv(sprintf("%s/4_Testing_Weather_Data_2022.csv", paths[["submission_data"]])) %>% 
    bind_rows(read.csv(sprintf("%s/4_Training_Weather_Data_2014_2021.csv", paths[["source_data"]]))) %>%
    pivot_longer(cols = c(3:18), names_to = "variable", values_to = "val") %>%
    filter(!is.na(val)) %>% mutate(month = substr(Date, 5, 6), ec = paste0(variable, "_" ,month)) 
    
  weather_data_bad <- weather_data_sub %>% count(Env, month, variable) %>% 
    arrange(Env, variable) %>% filter(n < 20)
  weather_data_sub_fil <- weather_data_sub %>% anti_join(weather_data_bad, by = c("Env", "variable", "month")) %>%
    group_by(Env, ec) %>% summarize(mean_val = mean(val), .groups = "drop") %>% 
    pivot_wider(id_cols = "Env", names_from = "ec", values_from = mean_val)
  
  miss_cols <- apply(weather_data_sub_fil[, -1], 2, function(x) sum(is.na(x)))
  weather_data_sub_fil_no_miss <- weather_data_sub_fil[names(which(miss_cols == 0))]
  
  weather_data_scaled_mat <- scale(weather_data_sub_fil_no_miss[, -1], center = T, scale = T) 
  rownames(weather_data_scaled_mat) <- weather_data_sub_fil$Env
  weather_data_kinship <- weather_data_scaled_mat %*% t(weather_data_scaled_mat)
  E_mat <- weather_data_kinship/(sum(diag(weather_data_kinship))/nrow(weather_data_kinship))
  
  # EC derived kinship
  #ec_train <- read.csv("~/g2f-maize-challenge-2022/source_data/Training_Data/6_Training_EC_Data_2014_2021.csv", header = T)
  #ec_sub <- read.csv("~/g2f-maize-challenge-2022/source_data/Testing_Data/6_Testing_EC_Data_2022.csv", header = T)
  #ec_data <- ec_train %>% bind_rows(ec_sub) %>% distinct(Env, .keep_all = T)
  #
  ##sl_data <- data[, sl] %>% distinct()
  ##wt_data <- data[, wt] %>% pivot_longer(!wt_dta_Env, names_to = "names", values_to = "value") %>% 
  ##  mutate(variable = gsub("(.*)\\_\\d{4}", "\\1", names, perl = T), 
  ##         day = gsub(".*\\_(\\d{4})", "\\1", names, perl = T)) %>%
  ##  select(-names) %>% pivot_wider(id_cols = c("wt_dta_Env", "day"), names_from = "variable",
  ##                                 values_from = "value")
  #
  #ec_scaled <- scale(ec_data[, -1], scale = T, center = T)
  #rownames(ec_scaled) <- ec_data$Env
  #
  #miss_info <- apply(ec_scaled, 2, function(x) sum(is.na(x)))
  #non_miss_cols <- miss_info[which(miss_info == 0)]
  #
  #ec_no_miss <- ec_scaled[, names(non_miss_cols)]
  #
  #ec_mat <- ec_no_miss %*% t(ec_no_miss)
  #
  #E_mat <- ec_mat/(sum(diag(ec_mat))/nrow(ec_mat))
  #
  #E_mat_imputed <- CovComb(Klist = list(E_mat, weather_data_kinship_scaled))
  
  qsave(E_mat, E_mat_at)
} else {
  E_mat <- qread(E_mat_at)
}

## genetic data
G_mat_at <- sprintf("%s/G_mat.qs", paths[["processed_data"]])
if (!file.exists(G_mat_at)){
  geno_data <- fread(sprintf("%s/geno_processed.miss.1.mac.1.biallelic.txt", paths[["processed_data"]]), 
                     header = T, sep = ",")
  
  geno_data_df <- setDF(geno_data)
  rm(geno_data)
  
  geno_data_mat <- as.matrix(geno_data_df[, 2:ncol(geno_data_df)])
  rownames(geno_data_mat) <- geno_data_df$Hybrid
  
  G_mat <- Gmatrix(geno_data_mat, method = "VanRaden")
  qsave(G_mat, G_mat_at)
} else {
  G_mat <- qread(G_mat_at)
}

#D.matrix <-function(M){
#  # from SnpReady
#  N <- nrow(M)
#  m <- ncol(M)
#  p <- colMeans(M)/2
#  D <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
#  return(D) 
#}
#
#geno_data_dom <- D.matrix(geno_data_mat)
#
#D_mat <- Gmatrix(geno_data_dom, method = "VanRaden")
#  

## format phenodata for preds
pheno_data <- pheno_data %>% filter(Env %in% rownames(E_mat))

## cv data
cv_data <- read_json(sprintf("%s/train_test_split_v4.json", paths[["processed_data"]]))
cv_data[["submission"]] <- list("train" = which(pheno_data$type == "train"), "val" = NA, "test" = which(pheno_data$type == "submission"))

## submission temp
orig_sub <- read.csv(sprintf("%s/1_Submission_Template_2022.csv", paths[["submission_data"]]))

# Generate matrices -------------------------------------------------------
# decomp for train data
eig_train <- generate_eigen_decomps(G_mat = G_mat,
                                    E_mat = E_mat,
                                    pheno_data = pheno_data[which(pheno_data$type == "train"), ],
                                    eig_g_at = sprintf("%s/Kg_eigen_train.qs", paths[["processed_data"]]),
                                    eig_e_at = sprintf("%s/Ke_eigen_train.qs", paths[["processed_data"]]),
                                    eig_ge_at = sprintf("%s/Kge_eigen_train.qs", paths[["processed_data"]])) # 5 hours

# decop for submission_data
eig_sub <- generate_eigen_decomps(G_mat = G_mat,
                                  E_mat = E_mat,
                                  pheno_data = pheno_data,
                                  eig_g_at = sprintf("%s/Kg_eigen_sub_v2.qs", paths[["processed_data"]]),
                                  eig_e_at = sprintf("%s/Ke_eigen_sub_v2.qs", paths[["processed_data"]]),
                                  eig_ge_at = sprintf("%s/Kge_eigen_sub_v2.qs", paths[["processed_data"]])) # 7 hours

#source("https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/getGEenriched.R")


# Predictions  -------------------------------------------------------

## cross validations

instance <- format(Sys.time(), format = "%H_%M_%d_%m_%y")
save_loc <- sprintf("%s/mixed_model/at_%s", paths[["dump_at"]], instance)
if(!dir.exists(save_loc)){dir.create(save_loc, recursive = T)}

## In parallel

cl <- makeCluster(length(cv_data)/2, outfile=sprintf("%s/%s", save_loc, "cluster.log"))
registerDoParallel(cl)
system.time(out_raw <- foreach(i=1:length(cv_data),
                               .packages = c("dplyr", "BGLR", "qs", "Metrics"))
            %dopar% yield_rmse(run = i,
                               cv_data = cv_data,
                               pheno_data = pheno_data[which(pheno_data$type == "train"), ],
                               eig_g_at = eig_train[["eig_g_at"]],
                               eig_e_at = eig_train[["eig_e_at"]],
                               eig_ge_at = eig_train[["eig_ge_at"]],
                               save_loc = save_loc,
                               nIter=15000, 
                               burnIn=2000))
stopCluster(cl)

results <- list()

results[["out_raw_df"]] <- as.data.frame(do.call(rbind, out_raw))

## for the submission set alone

pred_sub <- yield_rmse(run = 81,
                       cv_data = cv_data,
                       pheno_data = pheno_data,
                       eig_g_at = eig_sub[["eig_g_at"]],
                       eig_e_at = eig_sub[["eig_e_at"]],
                       eig_ge_at = eig_sub[["eig_ge_at"]],
                       save_loc = save_loc,
                       nIter=15000, 
                       burnIn=2000)
save_loc <- "~/g2f-maize-challenge-2022/dump/mixed_model/at_21_31_13_01_23"
read_data <- paste0(save_loc, "/results/", list.files(sprintf("%s/results", save_loc)))

# Read in results

all_res_files <- paste0(save_loc, "/results/", list.files(sprintf("%s/results", save_loc)))
results[["raw_data"]] <- do.call(rbind, lapply(all_res_files, qread))

## prepare submission
dennis_data <- read.csv("~/g2f-maize-challenge-2022/source_data/2168fc11-8a3d-472f-bcbd-7c380b489da9 (1).csv", 
                        header = T, col.names = c("Env", "Hybrid", "pred_dennis"))
to_submit <- results$raw_data %>% filter(type == "submission") %>% 
  select(Env, Hybrid, pred) %>% 
  right_join(orig_sub %>% select(-Yield_Mg_ha), by = c("Env", "Hybrid")) %>%
  rename(Yield_Mg_ha = pred) %>%
  left_join(dennis_data, by = c("Env", "Hybrid")) %>%
  mutate(Yield_Mg_ha = ifelse(is.na(Yield_Mg_ha), pred_dennis, Yield_Mg_ha)) %>%
  select(-pred_dennis)
  
write.csv(to_submit, "~/g2f-maize-challenge-2022/results/submit_v1_with_d_data.csv", row.names = F)

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
  ETA <- list(G = list(V=Kg_eigen$vectors,d=Kg_eigen$values, model='RKHS'),
              E = list(V=Ke_eigen$vectors,d=Ke_eigen$values, model='RKHS'),
              GE = list(V=Kge_eigen$vectors,d=Kge_eigen$values, model='RKHS'))
  
  cat(sprintf("ETA defined for run %s", run), 
      file = sprintf("%s/run_%s.log", log_at, run),
      sep = "\n",
      append = T)
  
  ## set BLAS threads
  RhpcBLASctl::blas_set_num_threads(1)
  
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
               "G" = model_fit$ETA$G$varU,
               "E"= model_fit$ETA$E$varU,
               "GxE"= model_fit$ETA$GE$varU,
               "err"= model_fit$varE,
               "Mean_RMSE" = mean_rmse,
               "Mean_cor" = mean_r2)
  
  return(out)
}

eigen_decomp <- function(mat, save_at){
  t1 <- Sys.time()
  eigen <- eigen(mat)
  t2 <- Sys.time()
  print(sprintf("decomp took %s minutes", difftime(t2, t1, units='mins')))
  qs::qsave(eigen, save_at)
  print(sprintf("saved at %s", save_at))
  return(eigen)
}
# paths -------------------------------------------------------------------

paths <- list(
  "tmp_at" = '/proj/g2f-maize-challenge-2022/tmp_data',
  "dump_at" = '/proj/g2f-maize-challenge-2022/dump',
  "source_data" = '/proj/g2f-maize-challenge-2022/source_data/Training_Data',
  "processed_data" = '/proj/g2f-maize-challenge-2022/processed_data',
  "results" = '/proj/g2f-maize-challenge-2022/processed_data'
)

# Load data ---------------------------------------------------------------

## phenodata
pheno_data_at <- sprintf("%s/pheno_data.qs", paths[["processed_data"]])
if (!file.exists(pheno_data_at)){
  
  data <- setDF(fread(sprintf("%s/combined_mat_v2.csv", paths[["processed_data"]])))
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
  data <- setDF(fread(sprintf("%s/combined_mat_v2.csv", paths[["processed_data"]])))
  all_cols <- colnames(data)
  ec <- all_cols[grep("ec_", all_cols)]
  ec_data <- data[, ec] %>% distinct()
  
  #sl_data <- data[, sl] %>% distinct()
  #wt_data <- data[, wt] %>% pivot_longer(!wt_dta_Env, names_to = "names", values_to = "value") %>% 
  #  mutate(variable = gsub("(.*)\\_\\d{4}", "\\1", names, perl = T), 
  #         day = gsub(".*\\_(\\d{4})", "\\1", names, perl = T)) %>%
  #  select(-names) %>% pivot_wider(id_cols = c("wt_dta_Env", "day"), names_from = "variable",
  #                                 values_from = "value")
  
  ec_scaled <- scale(ec_data[, -1], scale = T, center = T)
  rownames(ec_scaled) <- ec_data$ec_dta_Env
  
  miss_info <- apply(ec_scaled, 2, function(x) sum(is.na(x)))
  non_miss_cols <- miss_info[which(miss_info == 0)]
  
  ec_no_miss <- ec_scaled[, names(non_miss_cols)]
  
  ec_mat <- ec_no_miss %*% t(ec_no_miss)
  
  E_mat <- ec_mat/(sum(diag(ec_mat))/nrow(ec_mat))
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

## cv data
cv_data <- read_json(sprintf("%s/train_test_split_v3.json", paths[["processed_data"]]))

# Generate matrices -------------------------------------------------------
mat_names <- paste0(pheno_data$Env, ":", pheno_data$Hybrid)

# Generate design matrices
eig_g_at <- sprintf("%s/Kg_eigen.qs", paths[["processed_data"]])
if (!file.exists(eig_g_at)){
  Zg <- model.matrix(~ -1 + Hybrid, pheno_data)
  colnames(Zg) <- gsub('Hybrid', '', colnames(Zg), perl = T)
  G_mat_ordered <- G_mat[colnames(Zg), colnames(Zg)]
  Kg <- Zg %*% tcrossprod(G_mat_ordered, Zg)
  colnames(Kg) <- rownames(Kg) <- mat_names
  Kg_eigen <- eigen_decomp(Kg, eig_g_at)
}

eig_e_at <- sprintf("%s/Ke_eigen.qs", paths[["processed_data"]])
if(!file.exists(eig_e_at)){
  Ze <- model.matrix(~ -1 + Env, pheno_data)
  colnames(Ze) <- gsub('Env', '', colnames(Ze), perl = T)
  E_mat_ordered <- E_mat[colnames(Ze), colnames(Ze)]
  Ke <- Ze %*% tcrossprod(E_mat_ordered, Ze)
  colnames(Ke) <- rownames(Ke) <- mat_names
  Ke_eigen <- eigen_decomp(Ke, eig_e_at)
}

eig_ge_at <- sprintf("%s/Kge_eigen.qs", paths[["processed_data"]])
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
  Kge_eigen <- eigen_decomp(Kge, eig_ge_at)
}

#source("https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/getGEenriched.R")

# perform predictions

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
                               pheno_data = pheno_data,
                               eig_g_at = eig_g_at,
                               eig_e_at = eig_e_at,
                               eig_ge_at = eig_ge_at,
                               save_loc = save_loc))
stopCluster(cl)

results <- list()

results[["out_raw_df"]] <- as.data.frame(do.call(rbind, out_raw))

# Read in results

all_res_files <- paste0(save_loc, "/results/", list.files(sprintf("%s/results", save_loc)))
results[["raw_data"]] <- do.call(rbind, lapply(all_res_files, qread))

qsave(results, sprintf("%s/mixed_model_res.qs", paths[["results"]]))
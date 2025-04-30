library(argparse)
library(bigsnpr)
library(tidyverse)
library(tictoc)
library(ggplot2)
library(ggpubr)
library(yardstick)
library(pROC)
library(assertthat)
library(logr)
library(tidylog, warn.conflicts = FALSE)
library(glue)

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# https://stackoverflow.com/questions/68620638/how-to-execute-r-inside-snakemake
parser <- ArgumentParser(description = 'Launch the summary results')

parser$add_argument('--bedfile', help= 'name of the bed file without the extension. The script will check if <bedfile>.rds exists. If not, it will create it using snp_readBed' ) # nolint
parser$add_argument('--input-file', help= 'name of the input file containing the train and test samples together with the covariates' ) 
parser$add_argument("--snps", help = "file containing the list of the SNPs to perform bigsnpr on. E.g. clumped or pruned")
parser$add_argument("--diagnosis", help = "diagnosis name", type = "character") # e.g. skizo, skizospek, adhd, bipol, majordd,
parser$add_argument("--cohort", help = "cohort name", type = "character")
parser$add_argument("--control", help = "control name", type = "character") # nopsych_kontrol2015I
parser$add_argument("--k", help = "k fold cross validation", type = "double")
parser$add_argument("--covariates", help = "covariates to use", type = "character", nargs = "+", default = NULL)
parser$add_argument("--alphas", help = "alpha value for the ridge regression", type = "double", nargs = "+", default = NULL)

parser$add_argument("--variable", help = "variable to predict, can be either categorical or continuous", type = "character", choices = c("cat", "cont"))

parser$add_argument("--suffix", help = "suffix to add to the output file", type = "character", default = "")

parser$add_argument("--output", help = "output file", type = "character")
parser$add_argument("--input-type", help = "input type", type = "character", choices = c("onlygen", "onlycov", "gencov"))

parser$add_argument("--model-file", help = "model file", type = "character")
parser$add_argument("--predictions-file-test", help = "predictions file for test ", type = "character")
parser$add_argument("--predictions-file-train", help = "predictions file for train ", type = "character")
parser$add_argument("--metrics-file", help = "metrics file", type = "character")
parser$add_argument("--roc-curve-file", help = "roc curve file", type = "character")
parser$add_argument("--distr-cc-file", help = "distr cc file", type = "character")
# parser$add_argument("--confusion-matrix-file", help = "confusion matrix file", type = "character")
parser$add_argument("--valid-perf-file", help = "validation performance file", type = "character")






xargs <- parser$parse_args()



# Assign each of the xargs to a variable
bedfile <- xargs$bedfile
input  <- xargs$input_file
snps_file <- xargs$snps
diagnosis <- xargs$diagnosis
control <- xargs$control
column_diag <- paste0(diagnosis, "_", control)
covariates <- xargs$covariates
k <- xargs$k
alphas <- xargs$alphas
variable <- xargs$variable
suffix <- xargs$suffix
cohort <- xargs$cohort
input_type <- xargs$input_type
output_folder <- paste0(diagnosis, "_", control, "_", cohort, "_", input_type, "/")
model_name <- paste0(diagnosis, "_", control, "_", cohort, "_", input_type)
model_file <- xargs$model_file
predictions_file_train <- xargs$predictions_file_train
predictions_file_test <- xargs$predictions_file_test
metrics_file <- xargs$metrics_file
roc_curve_file <- xargs$roc_curve_file
distr_cc_file <- xargs$distr_cc_file
valid_perf_file <- xargs$valid_perf_file

log_print("Asserting stuff")

# assert that the previous variables are all character
assert_that(is.character(model_file))
assert_that(is.character(predictions_file_train))
assert_that(is.character(predictions_file_test))
assert_that(is.character(metrics_file))
assert_that(is.character(roc_curve_file))
assert_that(is.character(distr_cc_file))
assert_that(is.character(valid_perf_file))

date <- format(Sys.Date(), "%d_%m_%Y")
log  <- file.path(paste0("bigsnpr_", column_diag, "_", cohort, "_", input_type, "_", date, "_logr.log"))
lf  <- log_open(log, autolog = TRUE, show_notes = FALSE)

log_print("bedfile")
log_print(bedfile)
log_print(paste0(bedfile, ".rds"))

sep("bigsnpr logging")

log_print("Checking if .rds file exists")
log_print("Checking if .rds file exists")
# if bedfile + ".rds" file exists, read it with snp_attach and if not, read it with snp_readBed
if (file.exists(paste0(bedfile, ".rds"))) {
  log_print("File exists. Attaching it with snp_attach")
  bedfile <- snp_attach(paste0(bedfile, ".rds"))
} else {
  log_print("Reading bed file using snp_readBed")
  bedfile <- snp_readBed(bedfile, backingfile = "/home/leocob/igpv/backup/" + cohort + "/" + cohort + "_genotypes_plink/")
}

# if snps_file is not NULL, read it and assign it to snps
if (!is.null(snps_file)) {
  snps <- scan(snps_file, what = "character")
  log_print("Reading snps file")
} else {
  snps <- cols_along(X)
  log_print("No snps file provided. Using all the snps")
}

# extract only the clumped snps
bedfile_rsids  <- bedfile$map$marker.ID


# get indices of the clumped snps
bedfile_clumped_indices  <- which(bedfile_rsids %in% snps)



train_test <- read_csv(input, guess_max = Inf)  %>%  rename(family.ID = ID)



# Create mapping from row index to family.ID
rows_indices <- as_tibble(rows_along(bedfile$fam)) %>% rename(row_index = value)
sample_ids <- as_tibble(bedfile$fam %>% select(family.ID))
rows_to_id_mapping <- bind_cols(rows_indices, sample_ids)
# 80,771

train_test <- left_join(train_test, rows_to_id_mapping) %>% select(row_index, split, everything())

# sort train_test row_index
train_test <- train_test %>% arrange(row_index)

ind_train <- train_test %>% filter(split == "train") %>% select(row_index) %>% pull()
ind_test <- train_test %>% filter(split == "test") %>% select(row_index) %>% pull()

# ind.set creation for the kfold cross validation (CMSA)
valid_sets  <- train_test %>% 
filter(split == "train")  %>% 
mutate(valid_set = as.integer(str_extract(kfold, "(?<=_)[0-9]+"))) %>% 
pull(valid_set)

# Remove sex column from fam
bedfile$fam <- bedfile$fam[-5]

# Labels of the training and test set
y_train <- train_test %>% filter(split == "train") %>% pull(column_diag)
y_test <- train_test %>% filter(split == "test") %>% pull(column_diag)






# if y_train is categorical, run logistic regression
fun <- if (variable == "cat") big_spLogReg else big_spLinReg

covariates_train <- train_test %>% filter(split == "train") %>% select(all_of(covariates))
covariates_test <- train_test %>% filter(split == "test") %>% select(all_of(covariates))



covariates_fun_train <- if (input_type == "onlygen") NULL else covar_from_df(covariates_train)

covariates_fun_test <- if (input_type == "onlygen") NULL else covar_from_df(covariates_test)


if(input_type == "onlygen") {
  log_print(glue("No covariates used, as input_type is {input_type}"))
  covariates_fun_train <- NULL
  covariates_fun_test <- NULL
  pf_covar  <- NULL
  covariates  <- ""

} else if(input_type == "gencov") {
  log_print(glue("Input_type is {input_type} so covariates are used and are: {paste(colnames(covariates_train), collapse = ', ')}"))
  covariates_fun_train  <- covar_from_df(covariates_train)
  covariates_fun_test  <- covar_from_df(covariates_test)
  # Unpenalized covariates, following https://github.com/privefl/paper2-PRS/blob/master/response-snpnet/code/run-bigstatsr.R
  pf_covar  <- rep(0, ncol(covariates_fun_train))
  log_print(head(covariates_fun_train))
} else if(input_type == "onlycov") {
  # exit the script and throw an error saying it's not implemented yet
  stop(glue("Input_type is {input_type}, but this is not implemented yet"))
}


log_print("Starting penalized regression")

time <- system.time(
  model <- fun(X = bedfile$genotypes,
              y01.train = y_train,
              ind.train = ind_train,
              ind.sets = valid_sets,
              ind.col = bedfile_clumped_indices,
              covar.train = covariates_fun_train,
              pf.covar = pf_covar,
              alphas = alphas,
              K = k,
              warn = FALSE))

log_print("Penalized regression finished")
#saveRDS model with model_file

saveRDS(model, paste0(output_folder, "model_", model_name, ".rds"))
saveRDS(model, model_name)

# Save validation performances plot
valid_perf_plot <- plot(model)



#####################
# PREDICTION ON TEST SET
#####################
preds <- predict(model,
                 X = bedfile$genotypes,
                 ind.row = ind_test,
                 covar.row = covariates_fun_test)


# Save best model together with its best parameters
best_model_params <- summary(model, best.only = TRUE)

preds_train  <- predict(model,
                        X = bedfile$genotypes,
                        ind.row = ind_train,
                        covar.row = covariates_fun_train)

#####################
# Writing predictions to dataframe
#####################
predictions <- tibble(truth = y_test, pred = preds) %>% 
  mutate(truth = as.factor(truth),
         pred = 1 - pred,
         diag = model_name)  %>% 
         mutate(predicted = as.factor(ifelse(pred > 0.5, 1, 0))) %>% 
         bind_cols(row_index = ind_test) %>%
         left_join(rows_to_id_mapping) %>% 
         rename(ID = family.ID)

# double check that they are the same 
ind_test <- train_test %>% filter(split == "test") %>% select(row_index) %>% pull()
IDs_test  <- train_test %>% filter(split == "test") %>% pull(family.ID)
assert_that(sum(predictions$row_index %in% ind_test) == length(ind_test))
assert_that(sum(predictions$ID %in% IDs_test) == length(IDs_test))


predictions_train  <- tibble(truth = y_train, pred = preds_train) %>% 
  mutate(truth = as.factor(truth),
         pred = 1 - pred,
         diag = model_name)  %>% 
         mutate(predicted = as.factor(ifelse(pred > 0.5, 1, 0))) %>% 
         bind_cols(row_index = ind_train) %>%
         left_join(rows_to_id_mapping) %>% 
        rename(ID = family.ID)

#####################
# Calculate AUC with 10,000 bootstraps
#####################

# TODO: should have added y_test to the predictions dataframe
# Convert list to transposed dataframe and then to tibble
metrics_auc <- as_tibble(t(data.frame(AUCBoot(preds, y_test, nboot = 10)))) %>% 
mutate(model = model_name,
       k = k,
       covariates = paste(covariates, collapse = ","))


#####################
# Curves cases vs controls
#####################
distr_cases_controls <- qplot(preds, fill = as.logical(y_test),
      geom = "density", alpha = I(0.4)) +
  labs(fill = "Case?") +
  theme_bigstatsr() +
  theme(legend.position = c(0.52, 0.8)) +
  ggtitle(paste("bigsnpr", model_name)) +
  theme(plot.title = element_text(hjust = 0.5))




#####################
# ROC curve
#####################
roc_curve <- predictions %>%
  roc_curve(truth, pred) %>%
  ggplot(., aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  labs(
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(paste("bigsnpr", model_name)) +
      theme(plot.title = element_text(hjust = 0.5))




#####################
# SAVE PLOTS AND METRICS
#####################

# Save model
saveRDS(model, model_file)
log_print(paste(model_file, "written"))

# log_print(valid_perf_file, "going to be written")
ggsave(valid_perf_file, valid_perf_plot, useDingbats = FALSE, width = 9, height = 9)
# change extension of valid_perf_file from .pdf to .png
ggsave(gsub(".pdf", ".png", valid_perf_file), valid_perf_plot, width = 9, height = 9)
log_print(paste(valid_perf_file, "written"))


write_csv(metrics_auc, metrics_file)
log_print(paste(metrics_file, "written"))

ggsave(distr_cc_file, distr_cases_controls, useDingbats = FALSE, width = 9, height = 9)
ggsave(gsub(".pdf", ".png", distr_cc_file), distr_cases_controls, width = 9, height = 9)
log_print(paste(distr_cc_file, "written"))

# write predictions to csv file
write_csv(predictions, predictions_file_test)
log_print(paste(predictions_file_test, "written"))

# write predictions to csv file
write_csv(predictions_train, predictions_file_train)
log_print(paste(predictions_file_train, "written"))

ggsave(roc_curve_file, roc_curve, useDingbats = FALSE, width = 9, height = 9)
ggsave(gsub(".pdf", ".png", roc_curve_file), roc_curve, width = 9, height = 9)
log_print(paste(roc_curve_file, "written"))


log_print("script finished")

log_close()
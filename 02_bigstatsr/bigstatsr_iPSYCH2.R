#-
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
library(glue)
library(tidylog, warn.conflicts = FALSE)

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# https://stackoverflow.com/questions/68620638/how-to-execute-r-inside-snakemake
parser <- ArgumentParser(description = 'Launch the summary results')

parser$add_argument('--bedfile', help= 'name of the bed file without the extension. The script will check if <bedfile>.rds exists. If not, it will create it using snp_readBed' ) # nolint
parser$add_argument("--input-file", help= 'name of the input file containing the train and test samples together with the covariates')
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
parser$add_argument("--confusion-matrix-file", help = "confusion matrix file", type = "character")
parser$add_argument("--valid-perf-file", help = "validation performance file", type = "character")





# paste0("model_", model_name, ".rds"


xargs <- parser$parse_args()



# Assign each of the xargs to a variable
bedfile <- xargs$bedfile
# train <- xargs$train
# test <- xargs$test
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
# write command line execution of the r script in bash

# add an argparse for model, predictions, metrics, roc_curves, distr_cc, confusion_matrix
model_file <- xargs$model_file
predictions_file_train <- xargs$predictions_file_train
predictions_file_test <- xargs$predictions_file_test
metrics_file <- xargs$metrics_file
roc_curve_file <- xargs$roc_curve_file
distr_cc_file <- xargs$distr_cc_file
confusion_matrix_file <- xargs$confusion_matrix_file

log  <- file.path(paste0("bigsnpr_", column_diag, "_", cohort, "_", input_type, ".logr"))
lf  <- log_open(log, autolog = TRUE, show_notes = FALSE)
sep("bigsnpr logging")
model  <- readRDS(model_file)

log_print("bedfile")
log_print(bedfile)
log_print(paste0(bedfile, ".rds"))



log_print(glue("Checking if {bedfile}.rds file exists"))
# if bedfile + ".rds" file exists, read it with snp_attach and if not, read it with snp_readBed
if (file.exists(paste0(bedfile, ".rds"))) {
  log_print("File exists. Attaching it with snp_attach")
  bedfile <- snp_attach(paste0(bedfile, ".rds"))
} else {
  log_print("Reading bed file using snp_readBed")
  if (cohort == "iPSYCH2015i") {
    bedfile <- snp_readBed(glue("{bedfile}.bed"), backingfile = glue("/home/leocob/igpv/backup/{cohort}/{cohort}_genotypes_plink/merged_multiallelic/updated_rsid/"))
    # bedfile  <- snp_readBed(bedfile, backingfile = sub_bed(bedfile))
  } else if (cohort == "iPSYCH2012") {
    bedfile <- snp_readBed(glue("{bedfile}.bed"), backingfile = "/home/leocob/igpv/backup/" + cohort + "/" + cohort + "_genotypes_plink/")
  }
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



train_test <- read_csv(input, guess_max = Inf) %>% rename(family.ID = ID)

# Remove the columns from "F0000" to "R5010"



# Create mapping from row index to family.ID
rows_indices <- as_tibble(rows_along(bedfile$fam)) %>% rename(row_index = value)
sample_ids <- as_tibble(bedfile$fam %>% select(family.ID))
rows_to_id_mapping <- bind_cols(rows_indices, sample_ids)
# 80,771

train_test <- left_join(train_test, rows_to_id_mapping) %>% select(row_index, split, everything())

# sort train_test row_index
train_test <- train_test %>% arrange(row_index)


if(cohort == "iPSYCH2012"){
  ind_train <- train_test %>% filter(split == "train") %>% select(row_index) %>% pull()
  ind_test <- train_test %>% filter(split == "test") %>% select(row_index) %>% pull()

  # Labels of the training and test set
  y_train <- train_test %>% filter(split == "train") %>% pull(column_diag)
  y_test <- train_test %>% filter(split == "test") %>% pull(column_diag)
} else if (cohort == "iPSYCH2015i") {
  ind_train = NULL
  ind_test <- train_test %>% select(row_index) %>% pull()
  y_test <- train_test %>% pull(column_diag)
}

# ind.set creation for the kfold cross validation (CMSA)
valid_sets  <- train_test %>% 
filter(split == "train")  %>% 
mutate(valid_set = as.integer(str_extract(kfold, "(?<=_)[0-9]+"))) %>% 
pull(valid_set)

valid_sets_df  <- train_test %>% 
  filter(split == "train")  %>% 
  mutate(valid_set = as.integer(str_extract(kfold, "(?<=_)[0-9]+"))) %>% 
  rename(ID = family.ID) %>% 
  select(ID, row_index, valid_set)

# Remove sex column from fam
bedfile$fam <- bedfile$fam[-5]



truth_df  <- train_test %>% 
select(row_index, family.ID, column_diag) %>% 
rename(truth = column_diag) %>% 
mutate(truth = as.factor(truth))





# if y_train is categorical, run logistic regression
fun <- if (variable == "cat") big_spLogReg else big_spLinReg

# Here the training happened in an other script, on iPSYCH2012. So the covariates here are all of the test
covariates_train <- train_test %>% filter(split == "train") %>% select(all_of(covariates))
covariates_test <- train_test %>% filter(split == "test") %>% select(all_of(covariates))
covariates_train_test  <- train_test %>% select(all_of(covariates))

###############
# Comment out for real data
###############
# ind_train <- ind_train[1:500]
# ind_test <- ind_test[1:200]
# y_train <- y_train[1:500]
# y_test <- y_test[1:200]
# covariates_train <- covariates_train[1:500, ]
# covariates_test <- covariates_test[1:200, ]
# bedfile$genotypes <- as_FBM(bedfile$genotypes[,1:1000])
# input_type <- "onlygen"
########################################################


covariates_fun_train <- if (input_type == "onlygen") NULL else covar_from_df(covariates_train)

covariates_fun_test <- if (input_type == "onlygen") NULL else covar_from_df(covariates_test)

covariates_fun_train_test <- if (input_type == "onlygen") NULL else covar_from_df(covariates_train_test)

covariates_message  <- paste("The covariates are:", paste(colnames(covariates_fun_train), collapse = ", "))
log_print(covariates_message)
log_print("Starting penalized regression")


# Save validation performances plot
valid_perf_plot <- plot(model)


#####################
# PREDICTION ON TEST SET or iPSYCH2015i
#####################
preds <- predict(model,
                 X = bedfile$genotypes,
                 ind.row = ind_test,
                # ind.col = attr(model, "ind.col"),
                ind.col = bedfile_clumped_indices,
                 covar.row = covariates_fun_train_test)




row_ids_test  <- as.integer((attr(preds, "names")))


predictions  <- tibble(row_index = row_ids_test, pred = preds) %>% 
  mutate(one_minus_pred = 1 - pred,
          diag = model_name) %>% 
  mutate(predicted = as.factor(ifelse(pred > 0.5, 1, 0))) %>% 
  left_join(rows_to_id_mapping) %>% 
  left_join(truth_df) %>% 
  rename(ID = family.ID) %>% 
  select(ID, pred, one_minus_pred, truth, predicted, row_index, diag)

# write predictions to csv file
write_csv(predictions, predictions_file_test)
log_print(glue("predictions_file_test {predictions_file_test} written"))



#####################
# Calculate AUC with 10,000 bootstraps
#####################

# # TODO: should have added y_test to the predictions dataframe
# # Convert list to transposed dataframe and then to tibble
metrics_auc <- as_tibble(t(data.frame(AUCBoot(preds, y_test, nboot = 10)))) %>% 
# metrics_auc <- as_tibble(t(data.frame(AUCBoot(predictions$pred, as.integer(as.character(predictions$truth)))))) %>% 
mutate(model = model_name,
       k = k,
       covariates = paste(covariates, collapse = ","))


# #####################
# # Curves cases vs controls
# #####################
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
mutate(truth = as.factor(truth)) %>% 
  roc_curve(truth, one_minus_pred) %>%
  # roc_curve(truth, pred) %>%
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
    ggtitle(paste("Test bigsnpr", model_name)) +
      theme(plot.title = element_text(hjust = 0.5))



#####################
# SAVE PLOTS AND METRICS
#####################




write_csv(metrics_auc, metrics_file)
log_print(glue("metrics_file {metrics_file} written"))

ggsave(distr_cc_file, distr_cases_controls, useDingbats = FALSE, width = 9, height = 9)
ggsave(gsub(".pdf", ".png", distr_cc_file), distr_cases_controls, width = 9, height = 9)
log_print(glue("distr_cc_file {distr_cc_file} written"))

# write predictions to csv file
write_csv(predictions, predictions_file_test)
log_print(glue("predictions_file_test {predictions_file_test} written"))


ggsave(roc_curve_file, roc_curve, useDingbats = FALSE, width = 9, height = 9)
ggsave(gsub(".pdf", ".png", roc_curve_file), roc_curve, width = 9, height = 9)
log_print(glue("roc_curve_file {roc_curve_file} written"))


log_print("script finished")

log_close()
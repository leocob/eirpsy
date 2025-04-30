import os, os.path
# import ruamel.yaml
import yaml
import json


# Ideally I want a rule that writes the csv files for each of the diagnosis based on the control I want
# I want to create a list of directories such as diagnosis_control_cohort etc like in the eir-snakemake 
metadata_path = config.get("metadata_path")
# genotype_path = config.get("genotype_path")
# bedfilename = config.get("bedfilename")

control = config.get("control")
cohort = config.get("cohort")
diagnoses_to_run = config.get("diagnoses_to_run")
covariates = config.get("covariates")
test_size = config.get("test_size")
input_types = list(config.get("input_types"))
snps = config.get("snps")
genotype_path=os.path.join("/faststorage/jail/project/igpv/backup", cohort, cohort + "_genotypes_plink/")

if cohort == "iPSYCH2012":
    bedfile_name=cohort + ".PhaseBEAGLE5.1PhaseStates560ImputeBEAGLE5.1.SNP_SAMPLE_QC.UpdatedRSID.1.merged.biallelic"
    bedfile = genotype_path + bedfile_name
elif cohort == "iPSYCH2015i":
    bedfile_name=cohort + ".PhaseBEAGLE5.1PhaseStates560ImputeBEAGLE5.1.SNP_SAMPLE_QC.UpdatedRSID.1.merged.rsid_updated"
    bedfile = os.path.join(genotype_path, "merged_multiallelic/updated_rsid/", bedfile_name)
else:
    raise ValueError("Cohort not found")




# Create list of folders concatenating diagnosis, control,cat logs cohort, input_types
folders = []
for diagnosis in diagnoses_to_run:
    for input_type in input_types:
        # e.g. skizo_nopsychkontrol2015I_iPSYCH2012_onlygen
        folders.append(f"{diagnosis}_{control}_{cohort}_{input_type}")

# Create list of files contained in folders, concatenating diagnosis, control, cohort, input_types
output_files = []
for diagnosis in diagnoses_to_run:
    for input_type in input_types:
    # concatenate folder name with file name
        # output_files.append(f"{diagnosis}_{control}_{cohort}_{input_type}/model_{diagnosis}_{control}_{cohort}_{input_type}.rds")
        output_files.append(f"{diagnosis}_{control}_{cohort}_{input_type}/test_predictions_{diagnosis}_{control}_{cohort}_{input_type}.csv") 


wildcard_constraints:   
    diagnosis="skizo|autism|skizospek|majordd|bipol|adhd|singledd|recurdd",
    control="nodiagkontrol2015I",
    cohort="iPSYCH2012|iPSYCH2015i",
    input_type="onlygen|onlycov|gencov"

rule all:
    input:
        output_files

"""Run bigsnpr::big_spLogReg or big_spLinReg on the genotype data and the phenotype data"""
rule bigstatsr:
    input:
        input_file = metadata_path + "{cohort}_EurUnrel/{diagnosis}/all_splits_genotypes_{diagnosis}_{control}_{cohort}_EurUnrel.csv",
        model_file="{diagnosis}_{control}_iPSYCH2012_{input_type}/model_{diagnosis}_{control}_iPSYCH2012_{input_type}.rds"

    output:
        predictions_file_test="{diagnosis}_{control}_{cohort}_{input_type}/test_predictions_{diagnosis}_{control}_{cohort}_{input_type}.csv",
        predictions_file_train="{diagnosis}_{control}_{cohort}_{input_type}/train_predictions_{diagnosis}_{control}_{cohort}_{input_type}.csv",
        metrics_file = "{diagnosis}_{control}_{cohort}_{input_type}/metrics_{diagnosis}_{control}_{cohort}_{input_type}.csv",
        roc_curve_file = "{diagnosis}_{control}_{cohort}_{input_type}/roc_curve_{diagnosis}_{control}_{cohort}_{input_type}.pdf",
        distr_cc_file = "{diagnosis}_{control}_{cohort}_{input_type}/distr_cc_{diagnosis}_{control}_{cohort}_{input_type}.pdf"
        # confusion_matrix_file = "{diagnosis}_{control}_{cohort}_{input_type}/confusion_matrix_{diagnosis}_{control}_{cohort}_{input_type}.pdf",
        # valid_perf_file = "{diagnosis}_{control}_{cohort}_{input_type}/valid_perf_{diagnosis}_{control}_{cohort}_{input_type}.pdf"


    log:
        "logs/bigstatsr_preds_train_{diagnosis}_{control}_{cohort}_{input_type}.log"
    
    resources: 
        project = config["project"],
        partition = "normal",
        time = "11:55:00",
        mem = "50G",
        threads = "1"


    params:
        snps = config["snps"],
        bedfile = bedfile,
        covariates=config["covariates"],
        # diagnosis=wildcards.diagnosis,
        control=control,
        cohort=cohort,
        k=config["k"],
        alphas = config["alphas"],
        variable = config["variable"]
        # input_type = wildcards.input_type

    conda:
        "bigsnpr_env"


    shell:
        """
        mkdir -p {diagnosis}_{control}_{cohort}_{input_type}
        Rscript bigstatsr_pred_valids_iPSYCH2015i.R --bedfile {params.bedfile} \
        --input-file {input.input_file} \
        --snps {params.snps} \
        --input-type {wildcards.input_type} \
        --covariates {params.covariates} \
        --diagnosis {wildcards.diagnosis} \
        --control {params.control} \
        --cohort {params.cohort} \
        --k {params.k} \
        --alphas {params.alphas} \
        --variable {params.variable} \
        --model-file {input.model_file} \
        --predictions-file-test {output.predictions_file_test} \
        --predictions-file-train {output.predictions_file_train} \
        --metrics-file {output.metrics_file} \
        --roc-curve-file {output.roc_curve_file} \
        --distr-cc-file {output.distr_cc_file} \
        2> {log}
        """

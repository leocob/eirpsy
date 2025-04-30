"""
Started writing on the 02/01/2023
Snakemake pipeline to perform K-fold cross validation for hyperparameter tuning with EIR
"""


import os, os.path
import yaml
import glob
import json
from datetime import datetime
import sys
from helpers import determine_configs_to_write
metadata_path = config.get("metadata_path")
# Not used
# genotype_path = config.get("genotype_path")
# bedfilename = config.get("bedfilename")
# bedfile = genotype_path + bedfilename
control = config.get("control")
cohort = config.get("cohort")
diagnoses_to_run = config.get("diagnoses_to_run")
diagnoses_concat = "_".join(diagnoses_to_run) # Concatenate the diagnoses to run to use during Multi Task (MT) learning
covariates = config.get("covariates")
# input_types = list(config.get("input_type"))
k = config.get("kfold_splits")
# input_type = list(config.get("input_type").keys())[0] # this returns onlygen, onlycov or gencov
input_types = list(config.get("input_type").keys())
lof = config.get("lof")
num_combs = config.get("num_combs")

# If num_comb is 1, NO hyperparameter tuning will be performed. Write default batch_size to file. Needed for the pipeline
# if num_combs == 1:
#     with open("hyper_comb_1.txt", "w") as f:
#         f.write("--globals.batch_size=64")

range_hyper_combs = range(1, num_combs+1)
# range_hyper_combs = [26]

# join path of metadata_path and bedfilename
# bedfile = os.path.join(metadata_path, bedfilename)
############## for rule write_configs ##############
# configs_to_write, input_configs_to_write = determine_configs_to_write(input_type)
all_configs_to_write = set()
all_input_configs_to_write = set()
for input_type in input_types:
    configs_to_write, input_configs_to_write = determine_configs_to_write(input_type)
    # append each element of the list configs_to_write to the all_configs_to_write list


    all_configs_to_write.update(configs_to_write)
    all_input_configs_to_write.update(input_configs_to_write)


# print("input_type: ", input_type)
# print("configs_to_write: ", configs_to_write)
# print("input_configs_to_write: ", input_configs_to_write)
suffix_output_folder = config.get("suffix_output_folder")

task_type = config.get("task_type")
assert task_type in ["ST","MT"], "task_type should be either ST or MT"


best_hypers_dict = {
    "adhd": "hypercomb_27.txt",
    "autism": "hypercomb_15.txt",
    "bipol": "hypercomb_17.txt",
    "mdd": "hypercomb_15.txt",
    "skizo": "hypercomb_30.txt"
}

####################################################

wildcard_constraints:   
    diagnosis="skizo|autism|skizospek|majordd|bipol|adhd|recurdd|singledd|diab1|diab2|height|BMI|blength|bweight",
    # diagnosis="adhd",
    control="nodiagkontrol2015I|nodiag|nocontrol",
    cohort="iPSYCH2012_EurUnrel",
    input_type="onlygen|onlycov|gencov",
    threshold="0.05" 
# TODO: keep for Informed SNP set at 0.2, uncomment when running on subsets of SNPs


localrules: all, random_hyper_combs, average_validations, write_configs, process_configs_for_eirtrain_pretrained, create_combinations_input_configs, process_configs_for_eirpredict
rule all:
    input:

        expand("test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_global_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        tabdata = config["tabdata"],
        diagnosis = config["diagnoses_to_run"],
        control = config["control"],
        cohort = config["cohort"],
        task_type = config["task_type"],
        input_type = list(config["input_type"].keys())[0],
        suffix_output_folder = config["suffix_output_folder"],
        hypercomb=range_hyper_combs,
        fold=range(1, k+1),
        threshold=config["thresholds"]),
        expand("test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/eirpredict.done",
        tabdata = config["tabdata"],
        diagnosis = config["diagnoses_to_run"],
        control = config["control"],
        cohort = config["cohort"],
        task_type = config["task_type"],
        input_type = list(config["input_type"].keys())[0],
        suffix_output_folder = config["suffix_output_folder"],
        hypercomb=range_hyper_combs,
        fold=range(1, k+1),
        threshold=config["thresholds"])


#####################

# Set up random search of hyperparameters and write each combination to file
rule random_hyper_combs:
    # input:
    #     rules.k_fold_cv_splits.output
    output:
        expand("hyper_comb_{i}.txt", i=range_hyper_combs)

    resources:
        partition = "normal",
        time = "30:00",
        mem = "1G",
        project = config["project"],
        threads = 1

    params:
        num_combs = config.get("num_combs") # 5

    log:
        "logs/random_hyper_combs.log"
    # TODO: add conda environment
    # TODO: add log
    run:
        import random
        import numpy as np
        import pandas as pd

        combinations_df = pd.read_csv("/home/leocob/igpv/backup/metadata/hyperparameter_combinations.csv")

        # create dict from df
        combinations_dict = combinations_df.to_dict(orient='index')
        combinations_dict = {k+1:v for k,v in combinations_dict.items()}


                # Write each combination to a file
        for key, value in combinations_dict.items():
            with open(f'hyper_comb_{key}.txt', 'a') as file:
                # It's a nested dictionary so I need to iterate over the inner dictionary
                for hyper, option in value.items():
                    # if option == 0.0, convert it to 0
                    if option == 0.0:
                        option = 0
                    file.write(f'{hyper}={option}' + '\n')
#####################

# TODO: re-implement write_configs rule
# Write globals, input, fusion, output.yaml configs using default values from EIR's paper based on diagnosis, control, configname, cohort, task_type (ST or MT), input_type, suffix_output_folder
rule write_configs:
    input:
        configuration = "main_config.yaml"
        # task_option = "ST_option.txt"

    output:
        # ST learning config files
        # TODO: put it in a function
        expand("{diagnosis}_{control}_{cohort}_{task_type}_{configname}_{input_type}_{suffix_output_folder}.yaml", 
        diagnosis = config["diagnoses_to_run"], 
        control = config["control"], 
        # configname = configs_to_write,  
        configname = all_configs_to_write if input_type == "gencov" else all_configs_to_write - {"input_tabular_config"} if input_type == "onlygen" else all_configs_to_write - {"input_gln_config"},
        cohort = config["cohort"], 
        task_type = config["task_type"], 
        # input_type = list(config["input_type"].keys())[0],
        input_type = list(config["input_type"].keys()), 
        suffix_output_folder = config["suffix_output_folder"]) if config["task_type"] == "ST" else 
        # MT learning config files, where predict all diagnoses at once, so we need to write a config file having all the diagnoses
        expand("{diagnoses}_{control}_{cohort}_{task_type}_{configname}_{input_type}_{suffix_output_folder}.yaml", 
        diagnoses = diagnoses_concat,
        control = config["control"],
        configname = all_configs_to_write,
        cohort = config["cohort"],
        task_type = config["task_type"],
        input_type = list(config["input_type"].keys())[0], 
        suffix_output_folder = config["suffix_output_folder"])
    params: 
        globalconfig = config["global_config"],
        input_type = config["input_type"],
        fusionconfig = config["fusion_config"],
        outputconfig = config["output_config"],
        input_gln_config = config["input_gln_config"],
        input_tabular_config = config["input_tabular_config"],
        suffix_output_folder = config["suffix_output_folder"]

    resources:
        partition = "normal",
        time = "30:00",
        mem = "1G",
        project = config["project"],
        threads = "1",
        # threads = 1,
        # gres = False
    log:
        "logs/write_configs/write_configs_" + datetime.now().strftime("%Y-%m-%d.%H-%M-%S") + ".log"

    run:

        import sys

        def write_yaml_conf(output_config_name, config_file):
            with open(output_config_name + ".yaml", 'w') as outfile:
                yaml.dump(dict(config_file), outfile)


        if task_type == "ST":
            # List of output folders for the input of rule all
            # output_folders = [f"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}" for diagnosis in diagnoses_to_run]

            for diagnosis in diagnoses_to_run:
                for input_type in input_types:
                    configs_to_write, input_configs_to_write = determine_configs_to_write(input_type)
                    
                    for config_name in configs_to_write:

                        # TODO: should I add this? suffix_output_folder
                        # TODO: I'm actually not changing these in the config files
                        output_config_name = f"{diagnosis}_{control}_{cohort}_{task_type}_{config_name}_{input_type}_{suffix_output_folder}"

                        # Trial merging the branch and committing
                        # Trial to add in the new merged branch
                        if config_name == "globalconfig":
                            config[config_name]["output_folder"] = f"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}"

                        # Trial to add in the new merged branch
                        # if config_name == "outputconfig":
                            # params[config_name]["output_info"]["output_name"] = f"{diagnosis}_{control}"
                            # params[config_name]["output_info"]["output_source"] = "train_" + diagnosis + "_" + control + cohort + "EurUnrel.csv"

                            # params[config_name]["output_type_info"]["target_cat_columns"] = [diagnosis + "_" + control] # as a list
                            
                        if config_name == "input_tabular_config":
                            config[config_name]["input_info"]["input_source"] = os.path.join(metadata_path, cohort, diagnosis, "genotypes", "csv", "train_valid_" + diagnosis + "_" + control + "_" + cohort + ".csv")
                            # TODO: add LoF condition
                            # append all of the ENSGs to the main_config.yaml
                            # If lof ?? true
                            # params[config_name]["input_type_info"]["input_con_columns"] += main_config["lof_genes"]
                            # if lof is true
                            # if params[config_name]["lof"]:
                                # extract ENSG genes
                                # append them to the dictionary

                        # The json.loads(json.dumps()) is needed to convert a nested OrderedDict to a dict
                        # So that the yaml file can be written properly
                        # source: https://www.geeksforgeeks.org/how-to-convert-a-nested-ordereddict-to-dict/
                        write_yaml_conf(output_config_name, json.loads(json.dumps(config[config_name])))

        elif task_type == "MT":
            to_append = [control, cohort, task_type, suffix_output_folder]
            new_list = diagnoses_to_run + to_append
            output_folder = "_".join(new_list)

            diagnoses_concat = "_".join(diagnoses_to_run)
            for config_name in params.configs_to_write:
                                
                output_config_name = f"{diagnoses_concat}_{control}_{cohort}_{task_type}_{config_name}_{input_type}_{suffix_output_folder}"
                # output_config_name = diagnoses_concat + "_" + control + "_" + config_name + "_" + cohort + "_" + task_type + "_" + input_type

                if config_name == "globalconfig":
                    config[config_name]["output_folder"] = diagnoses_concat + "_" + control + "_" + cohort + "_" + task_type + "_" + input_type + "_" + suffix_output_folder

                if config_name == "outputconfig":
                    # params[config_name]["output_info"]["output_name"] = diagnoses_concat + "_" + control
                    # params[config_name]["output_info"]["output_source"] = "train_alldiag" + "_" + control + cohort + "EurUnrel.csv"

                    config[config_name]["output_type_info"]["target_cat_columns"] = [diagnosis + "_" + control for diagnosis in diagnoses_to_run]
                    
                if config_name == "input_tabular_config":
                    config[config_name]["input_info"]["input_source"] = "train_alldiag" + "_" + control + cohort + ".csv"

                
                write_yaml_conf(output_config_name, json.loads(json.dumps(params[config_name])))
#####################
# 


rule eirtrain:
    input:
        # "write_configs_done.txt",
        subset_snps_file = "/home/leocob/igpv/backup/iPSYCH2012_EurUnrel/genome_wide_PRS/psych_P_thresh_{threshold}_iPSYCH2012_EurUnrel_clumped_r2_0.5.rsids", # Informed SNP set 0.2 threshold
        # subset_snps_file = "/home/leocob/igpv/backup/isec_iPSYCH/isec_iPSYCH2015_clumped_200k.rsids" # Uninformed SNP set
        hyperparameters = "hyper_comb_{diagnosis}.txt",
        train = os.path.join(metadata_path, cohort, "{diagnosis}", "genotypes", "csv", "train_valid_{diagnosis}_{control}_{cohort}.csv"),
        valid_ids = os.path.join(metadata_path, cohort, "{diagnosis}", "genotypes", "csv", "valid_ids_{fold}_{diagnosis}_{control}_{cohort}.ids"),
        global_configs = "{diagnosis}_{control}_{cohort}_{task_type}_global_config_{input_type}_{suffix_output_folder}.yaml",
        fusion_configs = "{diagnosis}_{control}_{cohort}_{task_type}_fusion_config_{input_type}_{suffix_output_folder}.yaml",
        output_configs = "{diagnosis}_{control}_{cohort}_{task_type}_output_config_{input_type}_{suffix_output_folder}.yaml",
        input_configs = lambda wildcards: \
        [ 
            "{diagnosis}_{control}_{cohort}_{task_type}_input_gln_config_{input_type}_{suffix_output_folder}.yaml", 
            "{diagnosis}_{control}_{cohort}_{task_type}_input_tabular_config_{input_type}_{suffix_output_folder}.yaml"
        ] if wildcards.input_type == "gencov" else \
        [ 
            "{diagnosis}_{control}_{cohort}_{task_type}_input_gln_config_{input_type}_{suffix_output_folder}.yaml" 
        ] if wildcards.input_type == "onlygen" else \
        [ 
            "{diagnosis}_{control}_{cohort}_{task_type}_input_tabular_config_{input_type}_{suffix_output_folder}.yaml"
        ]
    
    conda:
        "eir0.1.39"

    log:
        "logs/eirtrain/{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}.log"

    resources:
        project = config["project"],
        time = "12:00:00",
        mem = "150G" if config["gpu_or_cpu"] == "gpu" else "50G",
        threads = 10 if config["gpu_or_cpu"] == "gpu" else 4,
        # threads = "4"
        # partition = "normal"
        partition = "gpu" if config["gpu_or_cpu"] == "gpu" else "normal",
        gres="gpu:1" if config["gpu_or_cpu"] == "gpu" else ""
        # Resources to try out on CPU
        # project = config["project"],
        # time = "30:00",
        # mem = "1G",
        # nodes = "1",


        # TODO: don't forget to use these!
        # Resources for actual EIR runs

    params:
        batch_size=config["global_config"]["batch_size"],
        device="cuda:0" if config["gpu_or_cpu"] == "gpu" else "cpu"
    output:
        # "{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}/logging_history.log"
        "{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/results/{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}/{diagnosis}_{control}/validation_{diagnosis}_{control}_history.log"
        
    shell:
        """
        num_samples=$(tail -n +2 {input.train} | wc -l)
        # I set the sample_interval so that for each epoch EIR will have seen all the samples
        sample_interval=$(expr $num_samples / {params.batch_size})

        echo "Number of samples is:" $num_samples
        echo "Sample interval is:" $sample_interval

        # read hyperparameters from the hyperparameters file
        hyperparameters=$(cat {input.hyperparameters})

        # remove the last 5 characters from the input files
        input_gln=$(echo {input.input_configs[0]} | sed 's/\.yaml//')
        global_configs=$(echo {input.global_configs} | sed 's/\.yaml//')
        fusion_configs=$(echo {input.fusion_configs} | sed 's/\.yaml//')
        output_configs=$(echo {input.output_configs} | sed 's/\.yaml//')

        hyperparameters=$(echo $hyperparameters | awk -v input_toreplace="$input_gln" "{{gsub(/input_gln/, input_toreplace)}}1")
        hyperparameters=$(echo $hyperparameters | awk -v global_configs_toreplace="$global_configs" "{{gsub(/globals/, global_configs_toreplace)}}1")
        hyperparameters=$(echo $hyperparameters | awk -v fusion_configs_toreplace="$fusion_configs" "{{gsub(/fusion/, fusion_configs_toreplace)}}1")
        hyperparameters=$(echo $hyperparameters | awk -v output_configs_toreplace="$output_configs" "{{gsub(/output/, output_configs_toreplace)}}1")


        echo "Hyperparameters are:"
        echo $hyperparameters

        eirtrain \
        --global_configs {input.global_configs} \
        --$global_configs.sample_interval=$sample_interval \
        --$global_configs.checkpoint_interval=$sample_interval \
        --$global_configs.device={params.device} \
        --input_configs {input.input_configs} \
        --fusion_configs {input.fusion_configs} \
        --output_configs {input.output_configs} \
        --$input_gln.input_type_info.subset_snps_file={input.subset_snps_file} \
        --$output_configs.output_type_info.target_cat_columns="['{wildcards.diagnosis}_{wildcards.control}',]" \
        --$output_configs.output_info.output_source={input.train} \
        --$output_configs.output_info.output_name={wildcards.diagnosis}_{wildcards.control}_{wildcards.cohort}_{wildcards.task_type}_{wildcards.input_type}_{wildcards.suffix_output_folder}_hypercomb_{wildcards.hypercomb}_valid_{wildcards.fold} \
        --$global_configs.output_folder={wildcards.diagnosis}_{wildcards.control}_{wildcards.cohort}_{wildcards.task_type}_{wildcards.input_type}_{wildcards.suffix_output_folder}_hypercomb_{wildcards.hypercomb}_valid_{wildcards.fold}_thresh_{wildcards.threshold} \
        --$global_configs.manual_valid_ids_file={input.valid_ids} \
        $hyperparameters \
        &> {log}
        """



rule aggregate_validations:
    input:
        # logging_histories = rules.eirtrain.output,
        validation_files = expand("{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/results/{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}/{diagnosis}_{control}/validation_{diagnosis}_{control}_history.log",
            diagnosis=config["diagnoses_to_run"],
            control=config["control"],
            cohort=config["cohort"],
            task_type=config["task_type"],
            input_type=list(config["input_type"].keys())[0],
            hypercomb=range_hyper_combs,
            fold=range(1, k + 1),
            threshold=config["thresholds"],
            suffix_output_folder=config["suffix_output_folder"])
        
    output:
        predictions = "preds_kfold_valids_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}.csv"

    resources:
        partition = "normal",
        time = "30:00",
        mem = "1G",
        project = config["project"],
        nodes = "1",
        threads = 1
    

    log:
        "logs/aggregate_validations_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}.log"
    run:
        import os
        import re
        import pandas as pd
        import argparse


        main_directory = os.getcwd()

        # Initialize empty dataframe
        preds_combined = pd.DataFrame()
        for file in input.validation_files:
            run_name = file.split("/")[0]
            Diagnosis = run_name.split("_")[0]
            Control = run_name.split("_")[1]
            Cohort = run_name.split("_")[2]
            Task_type = run_name.split("_")[3]
            Input_type = run_name.split("_")[4]
            suffix_output_folder = run_name.split("_")[5]

            Hypercomb = run_name.split("_")[7]
            Fold = run_name.split("_")[9]
            Threshold = run_name.split("_")[11]
            SNPs = suffix_output_folder.split("~")[0]
            Model = suffix_output_folder.split("~")[1]
            Data = suffix_output_folder.split("~")[2]
            LoF = suffix_output_folder.split("~")[3]

            if Diagnosis in ["adhd", "skizo", "skizospek", "bipol", "autism", "majordd", "recurdd", "singledd", "diab1", "diab2"]:
                avg_metric = "roc-auc-macro"
            elif Diagnosis in ["height", "BMI", "blength", "bweight"]:
                avg_metric = "R2"


            df = pd.read_csv(file)

            # Find the column containing "roc-auc-macro"
            avg_metric_col = df.columns[df.columns.str.contains(avg_metric)][0]

            # Find the row index where this column is maximized
            max_avg_metric_row_index = df[avg_metric_col].idxmax()

            # Get the entire row with the maximum avg_metric value and extract the column "iteration"
            max_row = df.iloc[max_avg_metric_row_index]
            iteration_best_metric = int(max_row["iteration"])


            # go to folder "samples/iteration_best_metric"
            # Extract everything from the string file up to the last backslash
            folder = re.search(r'(.*)\/', file).group(1)
            os.chdir(os.path.join(folder, "samples", str(iteration_best_metric)))

            # Read the file "predictions.csv" as pandas dataframe
            preds = pd.read_csv("predictions.csv")

            # Add columns with the variables assigned above
            preds["Model"] = Model
            preds["Diagnosis"] = Diagnosis
            preds["Control"] = Control
            preds["Cohort"] = Cohort
            preds["SNPs"] = SNPs
            preds["Threshold"] = Threshold
            preds["Data"] = Data
            preds["Input_type"] = Input_type
            preds["Covariates"] = "None" if Input_type == "geno" else "sex + age + PC1:10"
            preds["LoF"] = LoF
            preds["Hypercomb"] = Hypercomb
            preds["Fold"] = Fold
            preds["Task_type"] = Task_type

            # Append the dataframe to the empty dataframe
            preds_combined = preds_combined.append(preds)

            os.chdir(main_directory)
        # Go back to the main directory
        os.chdir(main_directory)

        print("output file is: ", output)

        print("output[0] file is: ", output[0])
        # Save the dataframe as a CSV file
        preds_combined.to_csv(output.predictions, index=False)



rule create_combinations_input_configs:
    input:
        input_template = "template_input_config_tabs.yaml"

    output:
        tabdata_file = "{diagnosis}_{tabdata}_input_config.yaml"
    log:
        "logs/create_combinations_input_configs_{diagnosis}_{tabdata}.log"

    params:
        features_sets = config["tabdata"],
        diagnoses = config["diagnoses_to_run"],
        cohort = config["cohort"],
        metadata_path = config["metadata_path"]

    run:
        import os
        import yaml

        # Load the input YAML file
        with open(input.input_template, 'r') as yaml_file:
            input_config = yaml.safe_load(yaml_file)

        # diagnoses = params.diagnoses
        diagnosis = wildcards.diagnosis
        # features_sets = params.features_sets
        feature_set = wildcards.tabdata
        cohort = params.cohort
        metadata_path = params.metadata_path
        metadata_path = os.path.join(metadata_path, cohort)

        # Create a dictionary to map feature sets to corresponding columns
        feature_set_columns = {
            "bigstatsr": ["bigstatsr_adhd", "bigstatsr_autism", "bigstatsr_bipol", "bigstatsr_majordd", "bigstatsr_skizo"],
            "PGSp": ["PGS_adhd", "PGS_autism", "PGS_bipol", "PGS_majordd", "PGS_skizo"],
            "FGRSp": ["FGRS_adhd", "FGRS_autism", "FGRS_bipol", "FGRS_majordd", "FGRS_skizo"]
        }

        # Generate and write YAML files for each combination
        # for diagnosis in diagnoses:
            # print("########## {diagnosis} ##########".format(diagnosis=diagnosis))
            # print("Feature set: {feature_set}".format(feature_set=feature_set))
        columns_for_diagnosis = []  # Initialize an empty list for columns

        # Add columns based on feature_set
        for feature in feature_set.split('-'):
            # print("Feature: {feature}".format(feature=feature))
            features_to_add = feature_set_columns[feature]
            # Keep only the one matching the diagnosis
            features_to_add = [feature for feature in features_to_add if diagnosis in feature]
        
            
            columns_for_diagnosis.extend(features_to_add)



        input_config['input_info']['input_source'] = os.path.join(metadata_path, diagnosis, "genotypes", "csv", f"train_valid_{diagnosis}_nodiagkontrol2015I_{cohort}_prss_fgrss_bigstatsr.csv")
        new_yaml_data = {
            'input_info': input_config['input_info'],
            'input_type_info': {
                'input_con_columns': columns_for_diagnosis
            },
            'model_config': input_config['model_config']
        }

        # Create the filename based on diagnosis and features_set
        # filename = f"{diagnosis}_{feature_set}_input_config.yaml"

        # Write the new YAML file
        with open(output.tabdata_file, 'w') as new_yaml_file:
            yaml.dump(new_yaml_data, new_yaml_file, default_flow_style=False)


#TODO: put config files, hypercombs etc in different and tidy folders

# TODO: In the case of gencov, I don't know if the input_configs.yaml file should contain both the models. Need to test but I don't want to do it now 20/10/2023
# "/faststorage/jail/project/igpv/eirpsy/results/2023-12-04-EIR_GLN-best_5CV_all_Informed_genotypes_iPSYCH2012"
model_path = "/faststorage/jail/project/igpv/eirpsy/results/2023-12-04-EIR_GLN-best_5CV_all_Informed_genotypes_iPSYCH2012"
rule process_configs_for_eirtrain_pretrained:
    input:
        global_configs = os.path.join(model_path,"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/configs/global_config.yaml"),
        fusion_configs = os.path.join(model_path,"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/configs/fusion_config.yaml"),
        output_configs = os.path.join(model_path,"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/configs/output_configs.yaml"),
        input_configs = os.path.join(model_path,"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/configs/input_configs.yaml"),
        tabdata_file = "{diagnosis}_{tabdata}_input_config.yaml"

    output:
        global_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_global_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        fusion_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_fusion_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        output_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_output_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        input_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_input_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml"

    log:
        "logs/process_configs_for_eirtrain_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{tabdata}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.log"

    params:
        device = "cuda:0" if config["gpu_or_cpu"] == "gpu" else "cpu",
        old_suffix = config["suffix_output_folder"],
        new_suffix = config["suffix_output_folder_withtab"]

    run:
        import yaml

        with open(input.fusion_configs) as f:
            fusion_configs = yaml.load(f, Loader=yaml.FullLoader)

        with open(output.fusion_configs, 'w') as f:
            yaml.dump(fusion_configs, f)

        with open(input.global_configs) as f:
            global_configs = yaml.load(f, Loader=yaml.FullLoader)

        # Inject the device cpu or cuda:0 depending on if I want to run eirpredict on CPU or GPU. It's based on the config["gpu_or_cpu"] in the main_config.yaml file
        global_configs["device"] = params.device
        global_configs["output_folder"] = wildcards.tabdata + "_" + global_configs["output_folder"]
        global_configs["dataloader_workers"] = 4
        with open(output.global_configs, 'w') as f:
            yaml.dump(global_configs, f)

        # read input.output_configs
        with open(input.output_configs) as f:
            output_configs = yaml.load(f, Loader=yaml.FullLoader)

        output_configs = output_configs[0]
        input_output_output_source = output_configs["output_info"]["output_source"]


        # output_configs["output_info"]["output_name"] = "withtab_" + output_configs["output_info"]["output_name"]

        output_configs["output_info"]["output_name"] = wildcards.tabdata + "_" + output_configs["output_info"]["output_name"]



        # write output_configs as output.output_configs
        with open(output.output_configs, 'w') as f:
            yaml.dump(output_configs, f)

        # read input.input_configs
        with open(input.input_configs) as f:
            input_configs = yaml.load(f, Loader=yaml.FullLoader)


        if wildcards.input_type == "onlygen":
            with open(output.input_configs, 'w') as f:
                yaml.dump(input_configs[0], f)
        elif wildcards.input_type == "onlycov":
            with open(output.input_configs[0], 'w') as f:
                yaml.dump(input_configs[0], f)
        else:
            try:
                assert isinstance(input_configs, list), "input_configs should be a list"
                assert len(input_configs) == 2, "input_configs should have length 2 because it's gencov"
            except AssertionError as e:
                with open(log, 'w') as f:
                    f.write(e)
                    f.write("input_configs should be a list and have length 2 because it's gencov")
            with open(output.input_configs[0], 'w') as f:
                yaml.dump(input_configs[0], f)
            with open(output.input_configs[1], 'w') as f:    
                yaml.dump(input_configs[1], f)



"/faststorage/jail/project/igpv/eirpsy/results/2023-12-04-EIR_GLN-best_5CV_all_Informed_genotypes_iPSYCH2012"
rule eitrain_pretrained:
    input:
        global_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_global_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        fusion_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_fusion_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        output_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_output_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        input_configs = "withtab_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_input_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        input_configstabs = "{diagnosis}_{tabdata}_input_config.yaml"

    output:
        output_dir = directory("{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}"),
        output_touch = touch("{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/eitrain_pretrained.done")

    resources:
        project = config["project"],
        time = "60:00:00" if config["gpu_or_cpu"] == "cpu" else "12:00:00",
        mem = "150G" if config["gpu_or_cpu"] == "gpu" else "50G",
        threads = 10 if config["gpu_or_cpu"] == "gpu" else 4,
        # threads = "4"
        # partition = "normal"
        partition = "gpu" if config["gpu_or_cpu"] == "gpu" else "normal",
        gres="gpu:1" if config["gpu_or_cpu"] == "gpu" else ""
    
    params:
        model_file = os.path.join(model_path,"{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/saved_models/*.pt"),
        device = "cuda:0" if config["gpu_or_cpu"] == "gpu" else "cpu"

    log:
        "logs/eirtrain_pretrained_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.log"

    conda:
        "eir0.1.39"

    shell:
        """
        # remove the last 5 characters from the input files
        input_gln=$(echo {input.input_configs[0]} | sed 's/\.yaml//')
        global_configs=$(echo {input.global_configs} | sed 's/\.yaml//')
        fusion_configs=$(echo {input.fusion_configs} | sed 's/\.yaml//')
        output_configs=$(echo {input.output_configs} | sed 's/\.yaml//')

        mkdir -p {output.output_dir}
        eirtrain \
        --global_configs {input.global_configs} \
        --$global_configs.pretrained_checkpoint={params.model_file} \
        --$global_configs.strict_pretrained_loading=False \
        --input_configs {input.input_configs} {input.input_configstabs} \
        --fusion_configs {input.fusion_configs} \
        --output_configs {input.output_configs} \
        --$global_configs.device={params.device}  &> {log}
        """
# touch {output.output_touch}

rule process_configs_for_eirpredict:
    input:
        rules.eitrain_pretrained.output.output_touch

    output:
        global_configs = "test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_global_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        fusion_configs = "test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_fusion_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        output_configs = "test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_output_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        input_configs = "test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_input_config_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml",
        input_configstabs = "test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_input_configtabs_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.yaml"
    
    log:
        "logs/process_configs_for_eirpredict_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.log"

    params:
        device = "cuda:0" if config["gpu_or_cpu"] == "gpu" else "cpu",
        global_configs = rules.eitrain_pretrained.output.output_dir + "/configs/global_config.yaml",
        fusion_configs = rules.eitrain_pretrained.output.output_dir + "/configs/fusion_config.yaml",
        output_configs = rules.eitrain_pretrained.output.output_dir + "/configs/output_configs.yaml",
        input_configs = rules.eitrain_pretrained.output.output_dir + "/configs/input_configs.yaml"
        
    resources:
        partition = "normal",
        time = "30:00",
        mem = "1G",
        project = config["project"],
        threads = "1"
        # threads = 1,
        # gres = False
    run:
        import yaml

        with open(params.fusion_configs) as f:
            fusion_configs = yaml.load(f, Loader=yaml.FullLoader)

        with open(output.fusion_configs, 'w') as f:
            yaml.dump(fusion_configs, f)

        with open(params.global_configs) as f:
            global_configs = yaml.load(f, Loader=yaml.FullLoader)

        # Inject the device cpu or cuda:0 depending on if I want to run eirpredict on CPU or GPU. It's based on the config["gpu_or_cpu"] in the main_config.yaml file
        global_configs["device"] = params.device
        with open(output.global_configs, 'w') as f:
            yaml.dump(global_configs, f)

        # read input.output_configs
        with open(params.output_configs) as f:
            output_configs = yaml.load(f, Loader=yaml.FullLoader)

        output_configs = output_configs[0]
        train_valid_output_source = output_configs["output_info"]["output_source"]

        test_output_source = train_valid_output_source.replace("train_valid", "test")

        output_configs["output_info"]["output_source"] = test_output_source

        # write output_configs as output.output_configs
        with open(output.output_configs, 'w') as f:
            yaml.dump(output_configs, f)




        # read input.input_configs
        with open(params.input_configs) as f:
            input_configs = yaml.load(f, Loader=yaml.FullLoader)

        if wildcards.input_type == "onlygen" and wildcards.tabdata is None:
            with open(output.input_configs, 'w') as f:
                yaml.dump(input_configs[0], f)
        elif wildcards.input_type == "onlycov" and wildcards.tabdata is None:
            with open(output.input_configs[0], 'w') as f:
                yaml.dump(input_configs[0], f)
        elif wildcards.input_type == "onlygen" and not wildcards.tabdata is None:

            with open(output.input_configs, 'w') as f:
                yaml.dump(input_configs[0], f)

            input_configstab = input_configs[1]
            configtabs_input_info_input_source = input_configstab["input_info"]["input_source"]

            configtabs_input_info_input_source = configtabs_input_info_input_source.replace("train_valid", "test")
            input_configstab["input_info"]["input_source"] = configtabs_input_info_input_source

            with open(output.input_configstabs, 'w') as f:
                yaml.dump(input_configstab, f)



rule eirpredict:
    input:
        global_configs = rules.process_configs_for_eirpredict.output.global_configs,
        fusion_configs = rules.process_configs_for_eirpredict.output.fusion_configs,
        output_configs = rules.process_configs_for_eirpredict.output.output_configs,
        input_configs = rules.process_configs_for_eirpredict.output.input_configs,
        input_configtabs = rules.process_configs_for_eirpredict.output.input_configstabs
   
    output:
        output_dir = directory("test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}"),
        output_touch = touch("test_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/eirpredict.done")

    resources:
        project = config["project"],
        time = "12:00:00",
        mem = "150G" if config["gpu_or_cpu"] == "gpu" else "50G",
        threads = 10 if config["gpu_or_cpu"] == "gpu" else 4,
        # threads = "4"
        # partition = "normal"
        partition = "gpu" if config["gpu_or_cpu"] == "gpu" else "normal",
        gres="gpu:1" if config["gpu_or_cpu"] == "gpu" else ""
    
    params:
        model_file = "{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_valid_{fold}_thresh_{threshold}/saved_models/*.pt"

    log:
        "logs/eirpredict_{tabdata}_{diagnosis}_{control}_{cohort}_{task_type}_{input_type}_{suffix_output_folder}_hypercomb_{hypercomb}_thresh_{threshold}_valid_{fold}.log"

    conda:
        "eir0.1.39"

    shell:
        """
        mkdir -p {output.output_dir}
        eirpredict \
        --global_configs {input.global_configs} \
        --input_configs {input.input_configs} {input.input_configtabs} \
        --fusion_configs {input.fusion_configs} \
        --output_configs {input.output_configs} \
        --model_path {params.model_file} \
        --evaluate \
        --output_folder {output.output_dir} &> {log}
        """

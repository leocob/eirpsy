Code used in relation to the project presented in the preprint "Deep learning-based polygenic scores enhance generalizability of psychiatric disorders prediction, Cobuccio L. et al. 2025."


## 0_Set up

Softwares and versions used

* EIR version 0.1.39 https://eir.readthedocs.io/en/stable/
* Snakemake 7.18.2

Encode genotype data into numpy array with plink_pipelines https://github.com/arnor-sigurdsson/plink_pipelines

## 01_GLN

Run DL model architecture Genome-Local-Net (GLN) on genotype data

In folder 01_gln

```

snakemake --snakefile eir_pipeline.smk --configfilegln_config.yaml --latency-wait 120 --keep-incomplete --use-conda \
--cluster "sbatch -A igpv -t {resources.time} -p {resources.partition} --mem={resources.mem} --cpus-per-task {resources.threads} ${gres}" --jobs 50
```

# 02_Bigstatsr

Train bigstatsr on training of iPSYCH1 and predict on test iPSYCH1

In folder 02_bigstatsr

```

snakemake --snakefile bigstatsr_pipeline_test.smk --configfile config_bigstatsr_test.yaml --jobs 30 --latency-wait 120 --keep-incomplete --use-conda \
--cluster "sbatch -A igpv -t {resources.time} -p {resources.partition} --mem={resources.mem} --cpus-per-task {resources.threads} "
```

Use trained model to predict on iPSYCH2

```

snakemake --snakefile bigstatsr_pipeline_iPSYCH2.smk --configfile config_bigstatsr_iPSYCH2.yaml --jobs 30 --latency-wait 120 --keep-incomplete --use-conda \
--cluster "sbatch -A igpv -t {resources.time} -p {resources.partition} --mem={resources.mem} --cpus-per-task {resources.threads} "

```

## 03_Logistic regression integration

Follow code in 03_logistic_integration/logistic_integration.qmd


## 04_DL integration with EIR

In folder 04_DL_integration

```

snakemake --config gpu_or_cpu=$gpu_or_cpu --snakefile eir_pipeline_DL_integration.smk --configfile DL_integration_config.yaml --latency-wait 10 --keep-incomplete --use-conda \
--cluster "sbatch -A igpv -t {resources.time} -p {resources.partition} --mem={resources.mem} --cpus-per-task {resources.threads} ${gres}"
```

## 05_Plotting

Follow code in 05_plotting/plotting.qmd

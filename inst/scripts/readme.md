## Pre-processing steps

### 1. Pick a study to process and update HIPCCyto package

* Review [the status dashboard](https://hipccyto.s3.us-east-2.amazonaws.com/status.html) and pick a study to process.
* Make sure that the HIPCCyto package is up-to-date.

```sh
cd /fh/fast/gottardo_r/HIPCCyto/HIPCCyto
git pull
ml fhR
Rscript -e "devtools::install()"
```


### 2. Fetch files

```sh
cd /fh/fast/gottardo_r/HIPCCyto/data
ml fhR
Rscript -e "HIPCCyto::fetch_files('<STUDY_ACCESSION>')"
```

Replace `<STUDY_ACCESSION>` with the study acession you'd like to fetch and process.


### 3. Summarize files and review marker names and batch

```R
HIPCCyto:::summarize_study("<STUDY_ACCESSION>", input_dir = "/fh/fast/gottardo_r/HIPCCyto/data/<STUDY_ACCESSION>/ResultFiles/Flow_cytometry_result")
```

If you notice discrepancies in channel names or batch groups, update `data-raw/DATA.R` to specify custom processing parameters.

```sh
cd /fh/fast/gottardo_r/HIPCCyto/HIPCCyto
vim data-raw/DATA.R
Rscript -e 'source("data-raw/DATA.R"); devtools::install()'
```


### 4. Process study

```sh
inst/scripts/run.sh <STUDY_ACCESSION>
```

Replace `<STUDY_ACCESSION>` with the study acession you'd like to process. `run.sh` will launch a slurm job in gizmok node and store a log file in `/fh/fast/gottardo_r/HIPCCyto/logs` and gating sets in `/fh/fast/gottardo_r/HIPCCyto/data/<STUDY_ACCESSION>/GatingSets/<HIPCCYTO_VERSION>`. The intermediate files for debugging purposes will be stored in `/fh/scratch/delete10/gottardo_r/HIPCCyto/<STUDY_ACCESSION>`.


### 5. Review the QC reports and flag gates to impute

Review the study report (`/fh/fast/gottardo_r/HIPCCyto/data/<STUDY_ACCESSION>/GatingSets/<HIPCCYTO_VERSION>/study.html`) and the QC reports by gating set (`/fh/fast/gottardo_r/HIPCCyto/data/<STUDY_ACCESSION>/GatingSets/<HIPCCYTO_VERSION>/<GATINGSET_ACCESSION>/QC.html`). Flag any Lymphocytes gates that look incorrect and need to be imputed in "Gating/Lymphocytes gates" tab and "By Sample" tab.


### 6. Impute gates

In "Sample metadata" tab, click the "Make imputation code" button on top to create code that will impute Lymphocytes gates. Copy the code and paste it in the R console. The code will look something like this:

```R
gs = HIPCCyto:::impute_gates(
  gs_dir = '/fh/fast/gottardo_r/HIPCCyto/data/SDY301/GatingSets/v0.0.6/gs2',
  samples_to_impute = c(
    'AIRFV020V10_2.497320.fcs',
    'AIRFV020V11_2.497326.fcs',
    'AIRFV020V12_2.497349.fcs',
    'AIRFV030V11_2.497527.fcs',
    'AIRFV030V12_2.497533.fcs'
  )
)
```

`HIPCCyto:::impute_gates()` will backup the original gates and QC report for recovery, impute Lymphocytes gates in `sample_to_impute` using the nearest method by batch group by default, and save the imputed gating set, the new QC report, and the gate comparison report (`comparison.html`) in the gating set directory.


### 7. Review the QC reports

Review the QC reports once again.

If you don't like the imputed gates and would like to flag and impute differently, you can recover the original gates and the QC report:

```R
# in R console
HIPCCyto:::recover_original_gates(<GATINGSET_DIRECTORY_PATH>)
# for example
HIPCCyto:::recover_original_gates("/fh/fast/gottardo_r/HIPCCyto/data/SDY212/GatingSets/v0.0.6/gs2")
```

After recovering the original gates, repeat the steps 5 and 6.


### 8. Upload the QC reports to S3 bucket

```sh
/fh/fast/gottardo_r/HIPCCyto/HIPCCyto/inst/scripts/upload.R <STUDY_ACCESSION>
```


### 9. Update the status dashboard

```sh
/fh/fast/gottardo_r/HIPCCyto/HIPCCyto/inst/scripts/update_status.sh
```

# GeneLab bioinformatics processing pipeline for Agilent microarray data

> **This page holds an overview and instructions for how GeneLab processes Agilent microarray data.**  

---

**Date:**  July 5, 2022
**Revision:**   
**Document Number:** 

**Submitted by:**  


**Approved by:**  


---

- [GeneLab bioinformatics processing pipeline for Agilent microarray data](#genelab-bioinformatics-processing-pipeline-for-agilent-microarray-data)
- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [1a. Importing the Runsheets](#1a-importing-the-runsheets)
    - [1b. Importing the Raw Data](#1b-importing-the-raw-data)
  - [2. Trim/Filter Raw Data and Trimmed Data QC](#2-trimfilter-raw-data-and-trimmed-data-qc)
    - [2a. Trim/Filter Raw Data](#2a-trimfilter-raw-data)
---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|R|4.1.2|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.14.0|[https://bioconductor.org](https://bioconductor.org)|
|Limma|3.50.3|[https://bioconductor.org/packages/release/bioc/html/limma.html](https://bioconductor.org/packages/release/bioc/html/limma.html)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory.
> 
> All output files marked with a \# are published for each RNAseq processed dataset in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

---

## 1. Raw Data QC

<br>

### 1a. Importing the Runsheets

```r
dir = "<File_Path_Containing_Runsheets>"
GLDS_41_rs <- read.csv(file.path(dir, "Runsheet(GLDS-41).csv"), check.names = FALSE, fileEncoding = 'UTF-8-BOM') ## Outputs a dataframe
Factor_Value <- GLDS_41_rs[,"Factor Value[gravity]"]
```

**Parameter Definitions:**

- `dir` – variable that stores the path to the runsheet
- `Runsheet(GLDS-41).csv` – this is the specific runsheet for this data set
- `check.names = FALSE` - makes it so R doesn't check to see if the variable names are syntactically valid variable names. (If true, R could change the names of the variables)
- `fileEncoding = 'UTF-8-BOM'` - declares what encoding to be used on the file. This allows the character data to be re-encoded and removes the special characters that are inputted by R. 

**Input Data:**

- *.csv (runsheet data)

**Output Data:**

- No output files at this step

<br>

### 1b. Importing the Raw Data 

```r
datadir = "<File_Path_to_Input_Data>"
files = dir(path = file.path(datadir, "<Directory_containing_raw_data>"), pattern="*\\.txt$")

## Input: ".txt" files, Output: EListRaw (Raw expression levels)
raw_data <- limma::read.maimages(files, source="agilent", path = file.path(datadir, "<Directory_containing_raw_data>"), sep="\t", green.only = TRUE)
```

**Parameter Definitions:**

- `datadir` - variable containing file path to the raw data
- `pattern="*\\.txt$"` - allows for the detection of any file that ends with .txt within the directory
- `files` – variable that stores all of the ".txt" files at the given directory
- `source="agilent"` – this parameter specifies that we're reading in data that came from the Agilent platform. (There are other platforms you can specify with this flag)
- `path = file.path(datadir, "<Directory_containing_raw_data")` - this tells the function the file path to where the raw data is stored
- `sep = "\t"` - this tells the function what the data in the raw data files are separated by
- `green.only = TRUE` - this tells the function that we only need to read the green channel (Cys3) instead of both green and red. 

**Input Data:**

- *.txt (raw data files)

**Output Data:**

- No output files at this step

<br>

---

## 2. Trim/Filter Raw Data and Trimmed Data QC

<br>

### 2a. Trim/Filter Raw Data  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores NumberOfThreads \
  --phred33 \
  --illumina \ # if adapters are not illumina, replace with adapters used
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this paramater if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

```

**Parameter Definitions:**

- `--gzip` – compress the output files with `gzip`
- `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
- `--cores` - specify the number of threads available on the server node to perform trimming
- `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
- `--illumina` - defines the adapter sequence to be trimmed as the first 13bp of the Illumina universal adapter `AGATCGGAAGAGC`
- `--output_dir` - the output directory to store results
- `--paired` - indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
- `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastq.gz\# (trimmed reads)
- *trimming_report.txt\# (trimming report)

<br>
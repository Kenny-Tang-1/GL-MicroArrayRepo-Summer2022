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
  - [2. QC Raw Data](#2-qc-raw-data)
    - [2a. Plotting density of raw expression values](#2a-plotting-density-of-raw-expression-values)
    - [2b. Generating a pseudoimage](#2b-generating-a-pseudoimage)
    - [2c. Generating an MA plot](#2c-generating-an-ma-plot)
    - [2d. Generating a Foreground and Background plot](#2d-generating-a-foreground-and-background-plot)
  - [3. Normalization](#3-normalization)
    - [3a. Background correction](#3a-background-correction)
    - [3b. Normalizing between arrays](#3b-normalizing-between-arrays)
  - [4. QA of Normalized Data](#4-qa-of-normalized-data)
    - [4a. Density plot](#4a-density-plot)
    - [4b. MA Plot](#4b-ma-plot)
    - [4c. Boxplots](#4c-boxplots)
  - [5. Generate Probe Level and Gene Level Data](#5-generate-probe-level-and-gene-level-data)
    - [5a. Generating raw probe level data](#5a-generating-raw-probe-level-data)
    - [5b. Generate raw gene level data](#5b-generate-raw-gene-level-data)
    - [5c. Generate normalized probe level data](#5c-generate-normalized-probe-level-data)
    - [5d. Generate normalized gene level data](#5d-generate-normalized-gene-level-data)
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
levels = c("Control Groups" , "Weightlessness", "Gravity, Altered")
Treatment = factor(Factor_Value, levels = levels)
```

**Parameter Definitions:**

- `dir` – variable that stores the path to the runsheet
- `Runsheet(GLDS-41).csv` – this is the specific runsheet for this data set
- `check.names = FALSE` - makes it so R doesn't check to see if the variable names are syntactically valid variable names. (If true, R could change the names of the variables)
- `fileEncoding = 'UTF-8-BOM'` - declares what encoding to be used on the file. This allows the character data to be re-encoded and removes the special characters that are inputted by R. 
- `levels = levels` - designates the elements in the levels variable as the levels of the Treatment variable.
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

## 2. QC Raw Data

<br>

### 2a. Plotting density of raw expression values

```r
plotDensities(raw_data, log = TRUE, legend = "topright" , main = "Density of raw intensities for multiple arrays")
```

**Parameter Definitions:**

- `raw_data` – ElistRaw object that contains the raw expression values.
- `log = TRUE` - Plots the densities on the log2 scale.
- `legend = "topright"` - specifies where to place the legend on the graph
- `main = "Density of raw intensities for multiple arrays"` - specifies what the title of the graph should be

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

### 2b. Generating a pseudoimage
```r
num_rows <- nrow(raw_data)

find_factors <- function(num) {
    factors_list <- list()
    for (n in 1:num) {
        if((num %% n) == 0) {
            factors_list <- c(factors_list, n)
        }
    }
    return(factors_list)
}
factors <- find_factors(num_rows)

for (sample_name in colnames(raw_data$E)) {
    if (length(raw_data$printer) != 0) {
        imageplot(raw_data$E[,sample_name], layout=raw_data$printer, zlim = c((1.25*min(raw_data$E[,sample_name])),(0.75*max(raw_data$E[,sample_name]))), legend = TRUE, main = sample_name)
        }
    else {
        if (length(factors) %% 2 == 0) {
            rows <- factors[[length(factors)/2]]
            columns <- factors[[length(factors)/2+1]]
            }
        else {
            rows <- factors[[ceiling(length(factors)/2)]]
            columns <- factors[[ceiling(length(factors)/2)]]
            }
        imageplot(raw_data$E[,sample_name], layout=list(ngrid.r = 1, ngrid.c = 1, nspot.r = rows , nspot.c = columns), zlim = c((1.25*min(raw_data$E[,sample_name])),(0.75*max(raw_data$E[,sample_name]))), legend = TRUE, main = sample_name)
    }
} 
```

**Parameter Definitions:**
- `num_rows` - A variable that contains the number of rows (number of probesets) in the raw dataset.
- `find_factors` - A function that takes a number as an input and return the factors of that number. Returns a list of the factors.
- `factors` - A variable that stores the list of factors.
- `rows` - This is a variable that stores one of the 2 values from the middle of the factors list. This is used in the layout argument to specify the amount of rows in the array.
- `columns` - This is a variable that stores the other value from the middle of the factors list. This is used in the layout argument of imageplot to specify the amount of columns in the array.
- `raw_data$E[, sample_name]` - Uses the raw intensities (rows) for the specific column (sample_name) as the values for the plot.
- `layout = raw_data$printer` - The printer layout tells the function the arrangement of the spots on the microarray. This is used to tell the function the dimensions to use.
- `layout=list(ngrid.r = 1, ngrid.c = 1, nspot.r = rows , nspot.c = columns)` - If the raw data set does not contain a printer component, the code block will instead find the factors of the number of rows in the dataset and choose the factors that have the closest value to each other (middle values). Since any 2 factors that multiplied into the number of rows of the raw dataset would work, I chose to use the values in the middle to generate the most square arrangement.
- `zlim` - Specifying the extreme values of z to associate with colors low and high.
- `legend = TRUE` - Places the range of the z values at the bottom of the plot.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

### 2c. Generating an MA plot
```r 
count = 0
for (sample_name in colnames(raw_data$E)) {
    count = count + 1
    limma::plotMA(raw_data,array=count,xlab="Average log-expression",ylab="Expression log-ratio(this sample vs. others)", main = sample_name)
}
```
**Parameter Definitions:**
- `count` = Variable that will increase as each iteration of the forloop occurs. This allows me to change which array is being used to generate the FB plots.
- `sample_name` - A variable that holds the element of the list that we're iterating through. In this case, the list that we're iterating through are the column names of the raw_data$E component. These would be the names of the samples in our raw dataset.
- `raw_data` - EListRaw object that contains the raw intensities as well as information about each array.
- `array = count` - Array determines which column(sample) that is being fed into the function. A count was used to change the array in each iteration.
- `xlab = "Average log-expression"` - The x axis displays the average log-expression of the raw intensities of the sample.
- `ylab = "Expression log-ratio(this sample vs. others)"` - The y axis displays the Expression log-ratio of this sample versus the average across all other samples.
- `main = sample_name` - Determines the title of the plot. In thise case it will be the name of each sample.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

### 2d. Generating a Foreground and Background plot
```r 
count = 0
for (sample_name in colnames(raw_data$E)) {
    count = count + 1
    plotFB(raw_data, array = count, xlab = "log2 Background", ylab = "log2 Foreground", main = sample_name) 
}
```
**Parameter Definitions:**
- `count` = Variable that will increase as each iteration of the forloop occurs. This allows me to change which array is being used to generate the FB plots.
- `sample_name` - A variable that holds the element of the list that we're iterating through. In this case, the list that we're iterating through are the column names of the raw_data$E component. These would be the names of the samples in our raw dataset.
- `raw_data$E` - Foreground expression values held in the EListRaw object.
- `array = count` - Array determines which column(sample) that is being fed into the function. A count was used to change the array in each iteration.
- `xlab = "log2 Background"` - Label for the x axis.
- `ylab = "log2 Foreground"` - Label for the y axis.
- `main = sample_name` - Makes the title of each plot the name of the sample being used to generate the plot.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

---

## 3. Normalization 

<br>

### 3a. Background correction
```r
corrected_data <- limma::backgroundCorrect(raw_data, method = "normexp")
```

**Parameter Definitions:**
- `raw_data` - The dataset that the function will execute the function on
- `method = "normexp"` - Designates the method of background correction to be done on the data

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

### 3b. Normalizing between arrays
```r
norm_data <- normalizeBetweenArrays(corrected_data, method = "quantile")
```

**Parameter Definitions:**
- `corrected_data` - The background corrected data set
- `method = quantile` - Designates the method of normalization that will be done on the data. Quantile specifies that the empirical distribution of each column is identical.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

---

## 4. QA of Normalized Data

<br>

### 4a. Density plot
```r
plotDensities(norm_data, log = TRUE, legend = "topright" , main = "Density of raw intensities for multiple arrays")
```

**Parameter Definitions:**
- `norm_data` - Used the normalized data to generate the density plot.
- `log = TRUE` - Specifying that the plot will generated on the log2 scale.
- `legend = "topright"` - Specifying where the legend (containing the sample names) will be located on the plot.
- `main = "Density of raw intensities for multiple arrays"` - Title for the plot.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

### 4b. MA Plot
```r
count = 0 
for (sample_name in colnames(raw_data$E)) {
    count = count + 1
    limma::plotMA(norm_data, array = count, xlab = "Average Log-expression", ylab = "Expression Log-ratio (this sample vs. others)", main = sample_name)
}
```

**Parameter Definitions:**
- `count = 0` - Variable called count that is used to change which array generated the plot in each iteration.
- `norm_data` - Using the normalized data to generate the MA plot.
- `array = count` - Used the count variable to change which array generated the plot in each iteration.
- `xlab = "Average Log-expression"` - Title of the x-axis.
- `ylab = "Expression Log-ratio (this sample vs. others)"` - Title of the y-axis.
- `main = sample_name` - Title of the plot.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

### 4c. Boxplots
```r
boxplot(log2(raw_data$E)) # Comparing the raw data to the normalized data
suppressWarnings(boxplot(log2(norm_data$E)))
```

**Parameter Definitions:**
- `raw_data$E` - Raw intensities that are used to generate the boxplot.
- `norm_data$E` - Normalized intensities that are used to generate the boxplot.

**Input Data:**

- No input files needed at this step

**Output Data:**

- No output files at this step

<br>

## 5. Generate Probe Level and Gene Level Data

<br>

### 5a. Generating raw probe level data 
```r 
raw_intensities_df <- as.data.frame(raw_data$E, row.names = raw_data$genes$GeneName, col.names = colnames(raw_data$E))# Generates a dataframe containing the raw intensities with the genes as the rows and samples as column names.
raw_probe_level_data_df <- cbind(raw_data$genes, raw_intensities_df)
write.csv(raw_probe_level_data_df, file = file.path(dir, "/output_data/GLDS-41_output_data/Raw_probe_level_data-GLDS-41.csv"), row.names = FALSE)
```

**Parameter Definitions:**
- `raw_intensities_df` - Dataframe containing the raw intensities (raw_data$E).
- `row.names = raw_data$genes$GeneName` - Made each row name the gene name from the genes component of the raw data.
- `col.names = colnames(raw_data$E)` - Made each column the name of the sample.
- `raw_probe_level_data_df` - Dataframe containing the raw intensities as well information about each gene (GeneName, ProbeName, ProbeUID, Row, Col, etc.)
- `file = file.path(dir, "/output_data/GLDS-41_output_data/Raw_probe_level_data-GLDS-41.csv")` - Indicates where to output the file and what to name the csv file.
- `row.names = FALSE` - Tells R not to export the row names


**Input Data:**

- No input files needed at this step

**Output Data:**

- *.csv file (containing raw probe level data)

<br>

### 5b. Generate raw gene level data
```r
first_sample_raw_mean = paste0(colnames(raw_data$E)[1], "_mean")
sample_one_raw = colnames(raw_data$E)[1]
summarized_columns_raw <- raw_probe_level_data_df %>% group_by(GeneName) %>% summarize(!!first_sample_mean := mean(get(sample_one_raw))) # Using this data frame to merge with
for (sample_name_raw in colnames(raw_data$E)) {
    print(sample_name_raw)
    sample_name_mean = paste0(sample_name_raw,"_mean")
    sample_name_sd = paste0(sample_name_raw,"_sd")
    sample_name_median = paste0(sample_name_raw,"_median")
    summarized_columns_raw <- merge(summarized_columns_raw, raw_probe_level_data_df %>% group_by(GeneName) %>% summarize(!!sample_name_mean := mean(get(sample_name_raw)), !!sample_name_sd := sd(get(sample_name_raw)), !!sample_name_median := median(get(sample_name_raw)), .groups = "keep"))
    } # This is where I'm merging the other data frames to the first one, to generate 1 data frame with all of the columnes (Gene_name, sample_mean, sample_sd, sample_median)
write.csv(summarized_columns_raw, file = file.path(dir, paste0("/output_data/GLDS-41_output_data/Raw_gene_level_data-GLDS-41.csv")), row.names = FALSE)
```

**Parameter Definitions:**
- `first_sample_raw_mean` - Variable containing the name of first sample with _mean added. Used to name the column.
- `sample_one_raw` - Variable containing only the name of the first sample.
- `summarized_columns_raw` - A dataframe that will hold all of the summarized gene level information. 
- `sample_name_raw` - Name of the current sample.
- `sample_name_mean` - Variable that holds the current sample name followed by _mean. Used to name the column.
- `sample_name_sd` - Variable that holds the current sample name followed by _sd. Used to name the column.
- `sample_name_median` - Variable that holds the current sample name followed by _median. Used to name the column.
- `.groups = "keep"` - Tells the function to keep the grouping structure from raw_probe_level_data_df.
- `file = file.path(dir, paste0("/output_data/GLDS-41_output_data/Raw_gene_level_data-GLDS-41.csv"))` - Telling the function where to output the data and what to name the file.
- `row.names = FALSE` - Tells R not to import the row names into the csv files.

**Input Data:**

- No input files needed at this step

**Output Data:**

- *.csv file (containing raw gene level data)

<br>

### 5c. Generate normalized probe level data
```r
norm_intensities_df <- as.data.frame(norm_data$E, row.names = norm_data$genes$GeneName, col.names = colnames(norm_data$E))
norm_probe_level_data_df <- cbind(norm_data$genes, norm_intensities_df)
write.csv(norm_probe_level_data_df, file = file.path(dir, "/output_data/GLDS-41_output_data/Normalized_probe_level_data-GLDS-41.csv"), row.names = FALSE)
```

**Parameter Definitions:**
- `norm_intensities_df` - Variable containing the normalized intensities (norm_data$E).
- `row.names = norm_data$genes$GeneName` - Made each row name the gene name from the genes component of the normalized data.
- `col.names = colnames(norm_data$E)` - Made each column the name of the sample.
- `norm_probe_level_data_df` - Variable containing the normalized intensities as well information about each gene (GeneName, ProbeName, ProbeUID, Row, Col, etc.)
- `file = file.path(dir, "/output_data/GLDS-41_output_data/Normalized_probe_level_data-GLDS-41.csv")` - Indicates where to output the file and what to name the csv file.
- `row.names = FALSE` - Tells R not to export the row names


**Input Data:**

- No input files needed at this step

**Output Data:**

- *.csv file (containing normalized probe level data)

<br>

### 5d. Generate normalized gene level data
```r
first_sample_norm_mean = paste0(colnames(norm_data$E)[1], "_mean")
sample_one_norm = colnames(norm_data$E)[1]
summarized_columns_norm <- norm_probe_level_data_df %>% group_by(GeneName) %>% summarize(!!first_sample_norm_mean := mean(get(sample_one_norm)))
for (sample_name_norm in colnames(norm_data$E)) {
    print(sample_name_norm)
    sample_name_mean = paste0(sample_name_norm, "_mean")
    sample_name_sd = paste0(sample_name_norm, "_sd")
    sample_name_median = paste0(sample_name_norm, "_median")
    summarized_columns_norm <- merge(summarized_columns_norm, norm_probe_level_data_df %>% group_by(GeneName) %>% summarize(!!sample_name_mean := mean(get(sample_name_norm)), !!sample_name_sd := sd(get(sample_name_norm)), !!sample_name_median := median(get(sample_name_norm)), .groups = "keep"))
    }
write.csv(summarized_columns_norm, file = file.path(dir, "/output_data/GLDS-41_output_data/Normalized_gene_level_data-GLDS-41.csv"), row.names = FALSE)
```
**Parameter Definitions:**
- `first_sample_norm_mean` - Variable containing the name of first sample with _mean added. Used to name the column.
- `sample_one_norm` - Variable containing only the name of the first sample.
- `summarized_columns_norm` - A dataframe that will hold all of the summarized gene level information. 
- `sample_name_norm` - Name of the current sample.
- `sample_name_mean` - Variable that holds the current sample name followed by _mean. Used to name the column.
- `sample_name_sd` - Variable that holds the current sample name followed by _sd. Used to name the column.
- `sample_name_median` - Variable that holds the current sample name followed by _median. Used to name the column.
- `.groups = "keep"` - Tells the function to keep the grouping structure from raw_probe_level_data_df.
- `file = file.path(dir, paste0("/output_data/GLDS-41_output_data/Normalized_gene_level_data-GLDS-41.csv"))` - Telling the function where to output the data and what to name the file.
- `row.names = FALSE` - Tells R not to import the row names into the csv files.

**Input Data:**

- No input files needed at this step

**Output Data:**

- *.csv file (containing normalized gene level data)

<br>
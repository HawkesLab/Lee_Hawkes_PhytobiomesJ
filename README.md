# Lee_Hawkes_PhytobiomesJ

This repository includes data and code to supplement the manuscript by Lee & Hawkes titled "Plant and soil drivers of whole-plant microbiomes: variation in switchgrass fungi from coastal to mountain sites" for publication in Phytobiomes J.

## Data

There are four main data files (SAMmatrix.csv, ASVmatrix.csv, TAXmatrix.csv, and TAXtree.tre) and a metadata file (Metadata.xlsx) in the folder "data". 

- The *sample matrix ("SAMmatrix.csv")* includes information about each leaf, root, and soil sample collected from switchgrass plants in the field and associated environmental measurements collected at the plant-level and site-level. See the *metadata file "Metadata.xlsx"* for information about column headers.

- The *ASV matrix ("ASVmatrix.csv")* includes the number of amplicon sequence variant (ASV) reads for each sample (rows) and unique ASV (columns). The raw sequence data used to generate this matrix is available through the NCBI SRA database under accession number PRJNA648664 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA648664/).

- The *taxonomy matrix ("TAXmatrix.csv")* includes the taxonomic classification of each unique ASV as was used in the manuscript. See the *metadata file "Metadata.xlsx"* for information about column headers.

- The *taxonomy tree ("TAXtree.tre")* is a Newick-formated taxonomy-based hierarchy that was generated using the perl script "taxonomy_to_tree.perl" from Tedersoo et al (2008) and the taxonomy matrix from this study. Please the manuscript for details and full reference information.


## R Markdown scripts

Annotated R code to run the analyses can be found in R markdown scripts (*.Rmd) at the top of the repository directory. The easiest way to run the code is to clone or download the entire repository and use RStudio to interact wtih the files as a project by clicking on "Lee_Hawkes_PhytobiomesJ.Rmd". If working outside of the RStudio project environment, I would recommend setting your working directory to the root of this repository (Lee_Hawkes_PhytobiomesJ). Be aware that code within the R markdown files reference information in neighboring folders (i.e. code, data, data_intermediates), making the location of the Rmd in the working directory critical for the code to function properly.

Analyses in the R Markdown scripts are organized as follows:

### IllumFUN_makeASVmatrix.Rmd
Processes the raw sequence data to create the ASV matrix using the dada2 pipeline and lulu ASV curation.

### IllumFUN_Q0.Rmd
Summarizes fungal diversity based on the ASV matrix.

### IllumFUN_Q1.Rmd
Addresses the study question: Given within-plant habitat, how does fungal communities similarity vary with proximity in the landscape?

### IllumFUN_Q2*.Rmd
Addresses the study question: Which environmental variables underpin changes in fungal diversity?

#### ...a
Performs model selection to determine which environmental variables to include in structural equation model (SEM) for each within-plant habitat

#### ...bc
Constructs and evaluates SEMs

#### ...de
Examines ASV-environment relationships using Bayesian hierarchical modeling via the R package HMSC


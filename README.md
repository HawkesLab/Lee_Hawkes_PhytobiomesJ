# Lee_Hawkes_PhytobiomesJ

This repository includes four main data files and a metadata file to supplement the manuscript by Lee & Hawkes titled "Plant and soil drivers of whole-plant microbiomes: variation in switchgrass fungi from coastal to mountain sites" for publication in Phytobiomes J.

The sample matrix ("SAMmatrix.csv") includes information about each leaf, root, and soil sample collected from switchgrass plants in the field and associated environmental measurements collected at the plant-level and site-level. 

The ASV matrix ("ASVmatrix.csv") includes the number of amplicon sequence variant (ASV) reads for each sample (rows) and unique ASV (columns).

The taxonomy matrix ("TAXmatrix.csv") includes the taxonomic classification of each unique ASV as was used in the manuscript. 

The taxonomy tree ("TAXtree.tre") is a Newick-formated taxonomy-based hierarchy that was generated using the perl script "taxonomy_to_tree.perl" from Tedersoo et al (2008) and the taxonomy matrix from this study. Please the manuscript for details and full reference information.

Detailed information about column headers in the sample matrix and taxonomy matrix can be found in the metadata file "Metadata.xlsx". 

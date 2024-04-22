# MRA_cystic-fibrosis
Here is described the code used in the analysis and visualization of metaproteomics data of cystic fibrosis microbiota samples using Rstudio. [Phyloseq package](https://joey711.github.io/phyloseq/) was used as the main tool.

## Creating the physeq object
The raw MS proteomics data have been submitted to the [ProteomeXchange Consortium](http://www.proteomexchange.org) via the Proteomics Identifications Database (PRIDE) partner repository, with the database identifier [PXD029284](https://www.ebi.ac.uk/pride/archive/projects/PXD029284) [^1].
The peptide, protein and taxonomic identification and quantification was performed using Metalab software. Here, the **BuiltIn_taxa_refine** file from Metalab results is used to performed the analysis in Rstudio.
- Load the different packages you need to use
  ```
  library(readr)
  library(phyloseq)
  library(tidyverse)
  ```
- Load the file in Rstudio
  ```
  BuiltIn_taxa_refine <- read_delim("K:/BuiltIn.taxa.refine.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE)
  View(BuiltIn_taxa_refine)
  ```
  ![BuiltIn_taxa_refine](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/1d7ed589-e2d4-43b0-bd8b-e997ab85cc0a)

- Create the **tax, otu and sam-tables**
  1. tax_table. Data matrix with the taxonomic information.
  2. otu_table. Data matrix with the intensity information of each taxa (row) in each sample (column).
  3. We also need to add the metadata with the sample information as a new object (sampledata).
     It is created outside and loaded into the environment.

[^1]: Perez-Riverol Y, Csordas A, Bai J, Bernal-Llinares M, Hewapathirana S, Kundu DJ, Inuganti A, Griss 142 J, Mayer G, Eisenacher M, Pérez E, Uszkoreit J, Pfeuffer J, Sachsenberg T, Yılmaz Ş, Tiwary S, Cox 143 J, Audain E, Walzer M, Jarnuczak AF, Ternet T, Brazma A, Vizcaíno JA. 2019. The PRIDE database 144 and related tools and resources in 2019: improving support for quantification data. Nucleic Acids 145 Research 47:D442-D450.

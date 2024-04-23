# MRA_cystic-fibrosis
Here is described the code used in the analysis and visualization of metaproteomics data of cystic fibrosis microbiota samples using Rstudio. [Phyloseq package](https://joey711.github.io/phyloseq/) was used as the main tool.

## Creating the physeq object
The raw MS proteomics data have been submitted to the [ProteomeXchange Consortium](http://www.proteomexchange.org) via the Proteomics Identifications Database (PRIDE) partner repository, with the database identifier [PXD029284](https://www.ebi.ac.uk/pride/archive/projects/PXD029284) [^1].
The peptide, protein and taxonomic identification and quantification was performed using Metalab software [^2]. Here, the **BuiltIn_taxa_refine** file from Metalab results is used to performed the analysis in Rstudio. This file contains the identified taxa with at least 2 unique peptides. 
- Load the different packages you need to use
  ```
  library(readr)
  library(phyloseq)
  library(tidyverse)
  library(microbiome)
  library(tidyverse)
  library(paletteer)
  library(ggplot2)
  library(scales)
  ```
- Load the file in Rstudio
  ```
  BuiltIn_taxa_refine <- read_delim("K:/BuiltIn.taxa.refine.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE)
  View(BuiltIn_taxa_refine)
  ```
  ![BuiltIn_taxa_refine](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/1d7ed589-e2d4-43b0-bd8b-e997ab85cc0a)

The file contains the taxa identified (rows) in the experiment, with their taxonomic ranks and abundance in each samples (columns). To performed the analysis `phyloseq` package needs this information to be split in two different datasets, `tax_table` and `otu_table`.
- Create the **tax, otu and sam-tables**
  1. **tax_table**. Data matrix with the taxonomic information.
  ```
  TAX_table <- BuiltIn_taxa_refine[-c(1:2, 11:40)]   # Create the dataframe from the BuiltIn_taxa_refine file with only taxonomic information
  TAX_table[is.na(TAX_table)] <- ""                  # Replace the NA values with an empty cell.
  ```
  Convert the table into the right format and then into a phyloseq object (using `tax_table`).
  ```
  TAX_table <- as.matrix(TAX_table)
  TAX = tax_table(TAX_table)
  View(TAX)
  ```
  ![TAX_table](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/d763a28b-95db-4ba0-9522-f02e036945fc)

  2. **otu_table**. Data matrix with the intensity information of each taxa (row) in each sample (column).
  ```
  OTU_table <- BuiltIn_taxa_refine[-c(1:10)]         # Create the dataframe from the BuiltIn_taxa_refine file with only abundance information in the samples
  ```
  Convert the table into the right format and then into a phyloseq object (using `otu_table`).
  ```
  OTU_table <- as.matrix(OTU_table)
  OTU = otu_table(OTU_table, taxa_are_rows = TRUE)
  View(OTU)
  ```
  ![OTU table](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/44972cf2-4b4b-4b88-8833-50f3cc87daf6)

Combine both objects into a new phyloseq object named physeq
```
physeq = phyloseq(OTU, TAX)
```
  3. We also need to add the metadata with the sample information as a new object (`sampledata`).
     It is created outside and loaded into the environment.
  ```
  metadata<- read.csv("SAMPLE_DATA.csv", header = TRUE, sep = ";")
  View(metadata)
  ```
  ![metadata](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/4e48f2e2-ab57-4540-8d37-ccdc7de51271)

  Convert the table into a phyloseq object and give the same sample names as the ones in the object `physeq`.
  ```
  sampledata <- sample_data(metadata)
  row.names(sampledata) <- sample_names(physeq)
  ```
  Join the new object into physeq
  ```
  physeq = phyloseq(OTU, TAX, sampledata)
  physeq
  ```
  You can see the different objects merge in `physeq` an their characteristics:nº of taxa, samples, sample variables and the different taxonomic ranks.
## Preprocessing
In order to continue the analysis a preprocessing step is needed to prepare the data.
- Order the `Time` variable in the `sampledata` object.
  ```
  timeorder <- c("Initial", "Intermediate", "Early_CF")
  sample_data(physeq)$Time <- factor(sample_data(physeq)$Time, levels = timeorder)
  ```
- Convert `Patient` and `Sample` into factor
  ```
  sample_data(physeq)$Sample <- factor(sample_data(physeq)$Sample)
  sample_data(physeq)$Patient <- factor(sample_data(physeq)$Patient)
  ```
- Remove *human* taxa
  ```
  human <- c("Chordata")
  physeq = subset_taxa(physeq, !Phylum %in% human)
  ```
- Remove those taxa with a sum intensity equal 0 (because we remove *human* taxa). Just in case.
  ```
  physeq <- filter_taxa(physeq, function(x) sum(x) > 0, TRUE)
  ```
- Split taxonomic ranks
  In the *BuiltIn_taxa_refine* file all taxa of different taxonomic ranks are mixed together so we need to split the ranks into different phyloseq objects. In this case we are only going to need *Phylum* and *Species* levels.
  - Phylum
    We need to select only those taxa that remain at phylum level (the sum intensity of the rest of the taxa that belong to this phylum is included in those).
    ```
    classNA = subset_taxa(physeq, Class=="")                        # Only taxa without a name in Class rank is selected
    filterphyla = c("")
    phylumnotNA = subset_taxa(classNA, !Phylum %in% filterphyla)    # We remove the taxa with empty names at phylum level
    tax_table(phylumnotNA)
    ```
    ![misclassification](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/588ca680-e633-4ef8-a525-3108df6eee6c)

    We view the result table and we need to eliminate any taxa that is wrongly classified. In this case for example *Firmicutes bacterium ASF500* is classified as phylum but it's a species with only phylum and species rank names (thats the reason of the misclassification), so we need to remove it.
    ```
    phylumnotNA = subset_taxa(phylumnotNA, !Species %in% "Firmicutes bacterium ASF500")
    tax_table(phylumnotNA)
    phylumnotNA
    ```
  - Species
    ```
    speciesNA = subset_taxa(physeq, Species=="")
    speciesnotNA = subset_taxa(physeq, !Species %in% filterphyla)
    tax_table(speciesnotNA)
    speciesnotNA
    ```
- Calculate the relative abundances of this new objects
  ```
  phylumnotNA_relativ <- transform_sample_counts(phylumnotNA, function(x) x / sum(x))
  speciesnotNA_relativ <- transform_sample_counts(speciesnotNA, function(x) x / sum(x))
  ```
## Selection of taxa for plotting
We cant represent in the graph all the taxa present in the table, so we are going to select based on the relative abundance, and only those taxa with a high abundance will be selected for the plot. The rest are going to be merge into a new category named "Other". 
- Phylum: only those taxa with a relative abundance of 1% or more.
  ```
  abund_table_phylum <- taxa_sums(phylumnotNA_relativ)                       # Obtain the abundance table
  namesother_phylum <- names(abund_table_phylum[abund_table_phylum<0.01])    # Filter those with an abundance of <1%
  namesother_phylum
  ```
  Merge of these taxa in another category named "Other"
  ```
  phylum_other <- merge_taxa2(phylumnotNA_relativ, taxa=namesother_phylum, name = "Other")
  tax_table(phylum_other)
  ```
  
- Species: only the 20 more abundant species.
  ```
  abund_table_species <- taxa_sums(speciesnotNA_relativ)
  namesother_species <- abund_table_species[order(abund_table_species)]
  namesother_species <- head(namesother_species, 113)                        # Select the ones that are going to be merged. Total nº of species minus 20
  namesother_species <- names(namesother_species)
  namesother_species
  ```
  Merge of these taxa in another category named "Other"
  ```
  species_other <- merge_taxa2(speciesnotNA_relativ, taxa=namesother_species, name = "Other")
  tax_table(species_other)
  ```
  This is a note
> [!IMPORTANT]  
> Sometimes this code return a tax_table with some errors in tax rank names. For example:
> ![tax_table with errors](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/7c9d42f0-f44c-4f3c-9af1-fb784d43885a)
>
> Instead of rename the rank names of the new category `Other` it rename the species *`sp212`* tax names and the ones in the `Other`category are now *NA*. So we new to correct the names in both.
>```
> tax_table(speciesnotNA_relativ)["sp212"]                                    # Identify the species wrongly named
> tax_tab <- tax_table(species_other)                                         # Extract the tax_table
> row_idx <- which(row.names(tax_tab) == "sp212")                             # Select the row of the sp212 and of the Other category
> row_idother <- which(row.names(tax_tab) == "Other")
> tax_tab[row_idx, "Superkingdom"] <- "Bacteria"                              # Remane each of the ranks with the correct names
> tax_tab[row_idx, "Kingdom"] <- ""
> tax_tab[row_idx, "Phylum"] <- "Actinobacteria"
> tax_tab[row_idx, "Class"] <- "Actinobacteria"
> tax_tab[row_idx, "Order"] <- "Bifidobacteriales"
> tax_tab[row_idx, "Family"] <- "Bifidobacteriaceae"
> tax_tab[row_idx, "Genus"] <- "Bifidobacterium"
> tax_tab[row_idx, "Species"] <- "Bifidobacterium scardovii"
> tax_tab[row_idother, "Species"] <- "Other"                                  # Rename the Species rank of the Other category
> tax_table(species_other) <- tax_tab                                         # Sustitute the tax_table in the object with the corrected one
> tax_table(species_other)                                                    # Check the new table.
> ```

  
## Plotting
- Phylum
  We have to prepare the data from phyloseq object for plotting using the function `psmelt`. And we put the new categry "Other" at the end.
  ```
  dfphylum <- psmelt(phylum_other)
  dfphylum$Phylum <- fct_relevel(dfphylum$Phylum, "Other", after = Inf)
  ```
  We only want to plot the *Initial* and *Early_CF* stages, so we create a new object with these two times, and we filter the dataset with it.
  ```
  times_to_show <- c("Initial", "Early_CF")
  dfphylum_filtrate <- dfphylum[dfphylum$Time %in% times_to_show, ]
  ```
  Create a new palette with the disired colors and plot the taxa.
  ```
  my_colors <- paletteer_d("ggthemes::Tableau_20")
  ggplot(dfphylum_filtrate, aes(x = Patient, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Time, scales = "free_x") +
    theme_bw()+
    scale_y_continuous(labels = percent)+
    ylab("Relative abundance (%)") +
    scale_fill_manual(name = "Phylum", values = c(my_colors))
  ```
- Species
  ```
  dfspecies <- psmelt(species_other)
  dfspecies$Species <- fct_relevel(dfspecies$Species, "Other", after = Inf)
  dfspecies_filtrate <- dfspecies[dfspecies$Time %in% times_to_show, ]

  ggplot(dfspecies_filtrate, aes(x = Patient, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Time, scales = "free_x") +
  theme_bw()+
  scale_y_continuous(labels = percent)+
  ylab("Relative abundance (%)") +
  scale_fill_manual(name = "Species", values = c(mi_paleta,"#D7CE9F"))+                 # We need to add one more color because the palette only has 20.
  theme(legend.text.align = 0, legend.text = element_text(face = "italic"))             # The names in italic because are species names
  ```
The two plots will look like this:
![Taxa representation](https://github.com/Carmengar/MRA_cystic-fibrosis/assets/71711674/1d91181d-6e5f-47da-9415-1ed0eabc9311)


[^1]: Perez-Riverol Y, Csordas A, Bai J, Bernal-Llinares M, Hewapathirana S, Kundu DJ, Inuganti A, Griss 142 J, Mayer G, Eisenacher M, Pérez E, Uszkoreit J, Pfeuffer J, Sachsenberg T, Yılmaz Ş, Tiwary S, Cox 143 J, Audain E, Walzer M, Jarnuczak AF, Ternet T, Brazma A, Vizcaíno JA. 2019. The PRIDE database 144 and related tools and resources in 2019: improving support for quantification data. Nucleic Acids 145 Research 47:D442-D450.
[^2]: Cheng K, Ning Z, Zhang X, Li L, Liao B, Mayne J, Stintzi A, Figeys D. 2017. MetaLab: an automated 121 pipeline for metaproteomic data analysis. Microbiome 5:157.

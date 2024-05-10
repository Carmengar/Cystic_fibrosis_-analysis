# Metaproteomics analysis of microbiome cystic fibrosis fecal samples

# We are going to use phyloseq package for the visualization

# Creating the physeq object

  # First, the taxonomy table "BuiltIn_taxa_refine" from Metalab results is imported to Rstudio:
library(readr)
BuiltIn_taxa_refine <- read_csv("K:/BuiltIn.taxa.refine.csv")
View(BuiltIn_taxa_refine)

  # Load the packages
library(phyloseq)
library(tidyverse)
  # We need to create several phyloseq objects from this taxonomy table. 
    # 1: tax_table. Data matrix with the taxonomic information.
TAX_table <- BuiltIn_taxa_refine[-c(1:2, 11:40)]
View(TAX_table)
      # Replace the NA values with an empty cell.
TAX_table[is.na(TAX_table)] <- ""
      # Convert the table into the right format and then into a phyloseq object (using "tax_table").
TAX_table <- as.matrix(TAX_table)
TAX = tax_table(TAX_table)

    # 2: otu_table. Data matrix with the intensity information of each taxa (row) in each sample (column).
OTU_table <- BuiltIn_taxa_refine[-c(1:10)]
View(OTU_table)
      # Convert the table into the right format and then into a phyloseq object (using "otu_table").
OTU_table <- as.matrix(OTU_table)
OTU = otu_table(OTU_table, taxa_are_rows = TRUE)

    # Combine both objects into a new phyloseq object named physeq
physeq = phyloseq(OTU, TAX)

    # 3: We also need to add the metadata with the sample information as a new object (sampledata).
      # It is created outside and imported into the environment.

metadata<- read.csv("SAMPLE_DATA.csv", header = TRUE, sep = ";")
      # Convert the table into a phyloseq object and give the same sample names as the ones in the object physeq.
sampledata <- sample_data(metadata)
row.names(sampledata) <- sample_names(physeq)
    # Join the new object into physeq
physeq = phyloseq(OTU, TAX, sampledata)
physeq

#########################################
# Preprocesing

  # Order the "Time" variable in the sampledata object.
timeorder <- c("Initial", "Intermediate", "Early_CF")
sample_data(physeq)$Time <- factor(sample_data(physeq)$Time, levels = timeorder)

  # Convert Patient and sample into factor
sample_data(physeq)$Sample <- factor(sample_data(physeq)$Sample)
sample_data(physeq)$Patient <- factor(sample_data(physeq)$Patient)
  # Remove human taxa
human <- c("Chordata")
physeq = subset_taxa(physeq, !Phylum %in% human)

  # Remove those taxa with a sum intensity equal 0 (because we remove human taxa). Just in case.
physeq <- filter_taxa(physeq, function(x) sum(x) > 0, TRUE)

  # Because in metaproteomics data, all taxonomics ranks are mixed together, we need to
  # split them into different phyloseq objects (in this case in phylum and species level).

    # Phylo level:
classNA = subset_taxa(physeq, Class=="")
filterphyla = c("")
phylumnotNA = subset_taxa(classNA, !Phylum %in% filterphyla)
       # View the result and eliminate those wrong
tax_table(phylumnotNA)
phylumnotNA = subset_taxa(phylumnotNA, !Species %in% "Firmicutes bacterium ASF500")
tax_table(phylumnotNA)
phylumnotNA

    # Species level
speciesNA = subset_taxa(physeq, Species=="")
speciesnotNA = subset_taxa(physeq, !Species %in% filterphyla)
tax_table(speciesnotNA)
speciesnotNA

  # Relative abundances
phylumnotNA_relativ <- transform_sample_counts(phylumnotNA, function(x) x / sum(x))
speciesnotNA_relativ <- transform_sample_counts(speciesnotNA, function(x) x / sum(x))

#########################################

# Relative abundances visualization
# We are not going to represent in the graph all the taxa.

# Load packages
library(microbiome)
library(tidyverse)
  # Phylum: only those taxa with a relative abundance of 1% or more.
abund_table_phylum <- taxa_sums(phylumnotNA_relativ)
namesother_phylum <- names(abund_table_phylum[abund_table_phylum<0.01])
namesother_phylum
      # Merge of these taxa in another category named "Other"
phylum_other <- merge_taxa2(phylumnotNA_relativ, taxa=namesother_phylum, name = "Other")
tax_table(phylum_other)

  # Species: only the 20 more abundant species.
abund_table_species <- taxa_sums(speciesnotNA_relativ)
namesother_species <- abund_table_species[order(abund_table_species)]
namesother_species <- head(namesother_species, 113) # Total nº of species minus 20
namesother_species <- names(namesother_species)
namesother_species
      # Merge of these taxa in another category named "Other"
species_other <- merge_taxa2(speciesnotNA_relativ, taxa=namesother_species, name = "Other")
tax_table(species_other)

# Visualization

  # Phylum
dfphylum <- psmelt(phylum_other)

    # "Other" at the end
dfphylum$Phylum <- fct_relevel(dfphylum$Phylum, "Other", after = Inf)

    # Construct the plot
library(paletteer)
library(ggplot2)
library(scales)

times_to_show <- c("Initial", "Early_CF")
dfphylum_filtrate <- dfphylum[dfphylum$Time %in% times_to_show, ]

my_colors <- paletteer_d("ggthemes::Tableau_20")

ggplot(dfphylum_filtrate, aes(x = Patient, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Time, scales = "free_x") +
  theme_bw()+
  scale_y_continuous(labels = percent)+
  ylab("Relative abundance (%)") +
  scale_fill_manual(name = "Phylum", values = c(my_colors))

  # Species
dfspecies <- psmelt(species_other)

      # "Other" at the end
dfspecies$Species <- fct_relevel(dfspecies$Species, "Other", after = Inf)

      # Plot
dfspecies_filtrate <- dfspecies[dfspecies$Time %in% times_to_show, ]

ggplot(dfspecies_filtrate, aes(x = Patient, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Time, scales = "free_x") +
  theme_bw()+
  scale_y_continuous(labels = percent)+
  ylab("Relative abundance (%)") +
  scale_fill_manual(name = "Species", values = c(mi_paleta,"#D7CE9F"))+
  theme(legend.text.align = 0, legend.text = element_text(face = "italic"))

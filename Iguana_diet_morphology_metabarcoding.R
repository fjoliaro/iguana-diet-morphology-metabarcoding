##################################### Venn Diagrams #######################################

#Make venn diagrams- iguana scat samples eliminated - each method

#import data from text file
iguana_scat_samples_eliminated_morph_meta <- read.delim("~/Desktop/Iguana_Diet/Iguana_diet_molec_visual/iguana_scat_samples_eliminated_morph_meta.txt")
View(iguana_scat_samples_eliminated_morph_meta)

#install ggvenn
install.packages("ggvenn")
library(ggvenn)
library(ggplot2)

#convert data to list of sets for ggvenn
list_of_sets_eliminated <- list(
  Metabarcoding = na.omit(unique(iguana_scat_samples_eliminated_morph_meta$Metabarcoding)),
  Morphology = na.omit(unique(iguana_scat_samples_eliminated_morph_meta$Morphology))
)


#make venn diagram
ggvenn(list_of_sets_eliminated, 
       stroke_size = 1, 
       set_name_size = 7,
       text_size = 10,
       show_percentage = FALSE, auto_scale = FALSE, fill_color = brewer.pal(n = 3, name = "Dark2")[1:2])


#Make venn diagram of plant species overlap between morphology and metabarcoding
#import plant lists for each method
iguana_scat_plant_list_venn <- read.delim("~/Desktop/Iguana_Diet/Iguana_diet_molec_visual/iguana_scat_plant_list_venn.txt")
View(iguana_scat_plant_list_venn)

#convert to plant list sets for ggvenn
plant_list_sets <- list(
  Metabarcoding = na.omit(unique(iguana_scat_plant_list_venn$Metabarcoding)),
  Morphology = na.omit(unique(iguana_scat_plant_list_venn$Morphology))
)

#make venn diagram
ggvenn(plant_list_sets, 
       stroke_size = 1, 
       set_name_size = 7,
       text_size = 10,
       show_percentage = FALSE, 
       auto_scale = TRUE) + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(expand = expansion(mult = 0.2))

##################### Paired sample analysis between morphology and metabarcoding #############################

#install packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)


#set wd
setwd("~/Desktop/Iguana_Diet/Iguana_diet_molec_visual_paired/")


#import visual and DNA data
iguana_scat_DNA_data <- t(read.delim("iguana_scat_DNA_pres_abs_genus_collapsed_numeric_paired.txt", row.names = 1)) #bringing in with the first column as the row name, then transforming so it has samples across the top and rownames are the taxonomy
iguana_scat_DNA_data # add a 'n' for nucleic acid
colnames(iguana_scat_DNA_data) <- paste0(colnames(iguana_scat_DNA_data),"n") # tack the n on the end with paste0
colnames(iguana_scat_DNA_data)

iguana_scat_visual_data <- t(read.delim("iguana_scat_visual_pres_abs_collapsed_numeric_paired.txt", row.names = 1))
iguana_scat_visual_data # sample names are the same so add a 'v' to the end so that we can bring in metadata to distinguish visual from DNA samples
colnames(iguana_scat_visual_data) <- paste0(colnames(iguana_scat_visual_data),"v") # tack the v on the end with paste0
colnames(iguana_scat_visual_data) 



## convert to phyloseq objects
DNA_otu_table <- otu_table(iguana_scat_DNA_data, taxa_are_rows = T)
DNA_otu_table
any(rowSums(iguana_scat_DNA_data == 0)) == ncol(iguana_scat_DNA_data)

visual_otu_table <- otu_table(iguana_scat_visual_data, taxa_are_rows = T)
visual_otu_table
any(rowSums(iguana_scat_visual_data == 0)) == ncol(iguana_scat_visual_data)

all_otu_table <- merge_phyloseq(DNA_otu_table, visual_otu_table) # combine them into one county matrix

# bring in metadata text file
metadata <- import_qiime_sample_data("IguanaScat_metadata_paired.txt") 

data <- merge_phyloseq(all_otu_table, metadata) # combine otu table and metadata
data
data@otu_table


# Calculate Jaccard distance matrix
jaccard.dist <- vegdist(t(data@otu_table), method = "jaccard")
jaccard.ord <- metaMDS(jaccard.dist) # ordinate with NMDS

############################## ANALYSIS METHOD JACCARD DISTANCE PERMANOVA ########################################

# create DF for plotting "Analysis" (morphology/metabarcoding)
jaccard.scrs <- scores(jaccard.ord$points) # extracting co-ords for each point (community)
jaccard.scrs_Analysis <- cbind(as.data.frame(jaccard.scrs), Analysis = data@sam_data$Analysis) # combine metadata with locations
jaccard.scrs_Analysis # now we've got the co-ords and metadata info
jaccard.cent_Analysis <- aggregate(cbind(MDS1,MDS2) ~ Analysis, data = jaccard.scrs_Analysis, FUN = mean) # get the co-ords for the centroids
jaccard.segs_Analysis <- merge(jaccard.scrs_Analysis, setNames(jaccard.cent_Analysis, c("Analysis","cMDS1","cMDS2")), by = "Analysis", sort = F)

data.df <- as(data@sam_data, "data.frame") # convert metadata to a DF for stats

#run permanova with adonis2 on jaccard distance by "Analysis" with strata for paired samples
adonis2(jaccard.dist ~ Analysis, data.df, strata = data.df$Sample) #significant, R2=0.01292, F=3.4807, p=0.001

#permanova including covariate "Island"
adonis2(jaccard.dist ~ Island + Analysis, data.df, strata = data.df$Sample) # R2=0.21481, F=7.8424, p=0.001

#including "Year"
adonis2(jaccard.dist ~ Year + Analysis, data.df, strata = data.df$Sample) #R2=0.05057, F=7.0579, p=0.001

#including "subspecies"
adonis2(jaccard.dist ~ Subspecies + Analysis, data.df, strata = data.df$Sample) #R2=0.04258, F=5.8934, p=0.001

#including "visitation"
adonis2(jaccard.dist ~ Visitation + Analysis, data.df, strata = data.df$Sample) #R2=0.04789, F=6.6652, p=0.001

#including all covariates
adonis2(jaccard.dist ~ Subspecies + Visitation + Year + Island + Analysis, data.df, strata = data.df$Sample) #R2= 0.23564, F=7.9229, p=0.001


#store and export adonis2 results
PERMANOVA_results_morph_meta_paired_covariates <- adonis2(jaccard.dist ~ Subspecies + Visitation + Year + Island + Analysis, data.df, strata = data.df$Sample)
PERMANOVA_results_df_morph_meta_paired_covariates <- as.data.frame(PERMANOVA_results_morph_meta_paired_covariates)
write.csv(PERMANOVA_results_df_morph_meta_paired_covariates, "PERMANOVA_results_morph_meta_paired_covariates.csv", row.names = TRUE)

#define grouping factor by Analysis method (morph/meta) and block by Sample
analysis_grouping <- data.df$Analysis
paired_block <- data.df$Sample

#run disperson with betadisper on jaccard distance by "Analysis"
permdisp_Analysis <- betadisper(jaccard.dist, group = analysis_grouping)

#run permutation test on dispersion
permutest(permdisp_Analysis, pairwise = TRUE, permutations = 999)


#store the permdisp result 
PERMDISP_results_morph_meta_paired <- permutest(permdisp_Analysis, pairwise = TRUE, permutations = 999)
PERMDISP_results_df_morph_meta_paired <- as.data.frame(PERMDISP_results_morph_meta_paired$tab)
write.csv(PERMDISP_results_df_morph_meta_paired, "PERMDISP_results_morph_meta_paired.csv", row.names = TRUE)

#boxplots of permdisp Analysis (morph/meta)
boxplot(permdisp_Analysis)

#PCoA of permdisp Analysis (morph/meta)
plot(permdisp_Analysis)

# plot Analysis (morph/meta) NMDS
ggplot(jaccard.scrs_Analysis) + theme_bw() +
  labs(y= "NMDS2", x= "NMDS1") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, colour = Analysis), size = 3) +
  scale_colour_brewer(palette = "Dark2")

# plot Analysis (morph/meta) NMDS with centroids

ggplot(jaccard.segs_Analysis) + theme_bw() +
  labs(y= "NMDS2", x= "NMDS1") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  scale_x_continuous(limits = c(-0.75,0.75)) +
  scale_y_continuous(limits = c(-0.75,0.75)) +
  geom_point(aes(x=MDS1, y=MDS2, colour = Analysis), size = 3) +
  geom_point(aes(x=cMDS1, y=cMDS2, colour = Analysis), size = 10) +
  scale_colour_brewer(palette = "Dark2") + theme(text = element_text(size = 20))


############################ ALPHA DIVERSITY (richness, Shannon)#############################

#calculate per-sample (alpha) diversity for both methods
alpha_diversity <- estimate_richness(data, measures = c("Observed", "Shannon"))
alpha_diversity_df <- cbind(alpha_diversity, data.df) #combine with metadata into dataframe
#summarize the data
alpha_div_summary <- summary(alpha_diversity_df)
alpha_div_summary
write.csv(alpha_div_summary, "alpha_diversity_summary.csv", row.names = TRUE)

# bring in metadata text file for DNA (metabarcoding) only paired samples
metadata_DNA <- import_qiime_sample_data("IguanaScat_metadata_DNA_paired.txt") 

data_DNA <- merge_phyloseq(DNA_otu_table, metadata_DNA) # combine DNA otu table and DNA metadata
data_DNA
data_DNA@otu_table

data.df_DNA <- as(data_DNA@sam_data, "data.frame") # convert DNA metadata to a DF for stats

#calculate alpha diversity for metabarcoding
alpha_diversity_DNA <- estimate_richness(data_DNA, measures = c("Observed", "Shannon"))
alpha_diversity_DNA_df <- cbind(alpha_diversity_DNA, data.df_DNA) #combine with metadata into dataframe
#summarize the metabarcoding data
alpha_div_summary_DNA <- summary(alpha_diversity_DNA_df)
alpha_div_summary_DNA
write.csv(alpha_div_summary_DNA, "alpha_diversity_summary_DNA_paired.csv", row.names = FALSE)

# bring in metadata text file for visual(morphology) only
metadata_visual <- import_qiime_sample_data("IguanaScat_metadata_visual_paired.txt") 

data_visual <- merge_phyloseq(visual_otu_table, metadata_visual) # combine visual otu table and visual metadata
data_visual
data_visual@otu_table

data.df_visual <- as(data_visual@sam_data, "data.frame") # convert DNA metadata to a DF for stats

#calculate alpha diversity for morphology
alpha_diversity_visual <- estimate_richness(data_visual, measures = c("Observed", "Shannon"))
alpha_diversity_visual_df <- cbind(alpha_diversity_visual, data.df_visual) #combine with metadata into dataframe
#summarize the morphology data
alpha_div_summary_visual <- summary(alpha_diversity_visual_df)
alpha_div_summary_visual
write.csv(alpha_div_summary_visual, "alpha_diversity_summary_visual_paired.csv", row.names = FALSE)

#run wilcoxon signed-rank test and paired t-test for observed richness between analysis methods
observed_metabarcoding <- alpha_diversity_df$Observed[alpha_diversity_df$Analysis == "Metabarcoding"]
observed_morphology <- alpha_diversity_df$Observed[alpha_diversity_df$Analysis == "Morphology"]
wilcox.test(observed_metabarcoding, observed_morphology, paired = TRUE)
t.test(observed_metabarcoding, observed_morphology, paired = TRUE)

#export wilcox signed-rank test results for observed richness
wilcox_observed <- wilcox.test(observed_metabarcoding, observed_morphology, paired = TRUE)
wilcox_observed_df <- data.frame(
  Statistic = wilcox_observed$statistic,
  P_value = wilcox_observed$p.value,
  Alternative_Hypothesis = wilcox_observed$alternative
)
write.csv(wilcox_observed_df, "wilcox_signed_rank_test_results_observed.csv", row.names = FALSE)

#run wilcoxon signed-rank test and t-test for shannon diversity between analysis methods
shannon_metabarcoding <- alpha_diversity_df$Shannon[alpha_diversity_df$Analysis == "Metabarcoding"]
shannon_morphology <- alpha_diversity_df$Shannon[alpha_diversity_df$Analysis == "Morphology"]
wilcox.test(shannon_metabarcoding, shannon_morphology, paired = TRUE)
t.test(shannon_metabarcoding,shannon_morphology, paired = TRUE)

#export wilcox test results for shannon
wilcox_shannon <- wilcox.test(shannon_metabarcoding, shannon_morphology, paired = TRUE)
wilcox_shannon_df <- data.frame(
  Statistic = wilcox_shannon$statistic,
  P_value = wilcox_shannon$p.value,
  Alternative_Hypothesis = wilcox_shannon$alternative
)


#################################TEST SIGNIFICANCE BETWEEN ANALYSIS METHOD FOR TOP TAXA####################################

# Melt the phyloseq object
melted_data <- psmelt(data)

#Group by Analysis method
melted_data$Group <- as.factor(melted_data$Analysis)


#get unique names of the taxa
taxa_names <- unique(melted_data$OTU)

#subset data for each top 12 taxon
Casasia.clusiifolia_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Gentianales_4055.f__Rubiaceae_24966.g__Casasia_58499.s__Casasia.clusiifolia_58500", ]
Rachicallis.americana_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Gentianales_4055.f__Rubiaceae_24966.g__Rachicallis_171717.s__Rachicallis.americana_1717177", ]
Manilkara.bahamensis_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Ericales_41945.f__Sapotaceae_3737.g__Manilkara_3740.s__Manilkara.bahamensis_37400", ]
Guaiacum.sanctum_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Zygophyllales_403666.f__Zygophyllaceae_43873.g__Guaiacum_45188.s__Guaiacum.sanctum_45189", ]
Guapira.discolor_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__undef__0.o__Caryophyllales_3524.f__Nyctaginaceae_3536.g__Guapira_427802.s__Guapira.discolor_427803", ]
Vitis.vinifera_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Vitales_403667.f__Vitaceae_3602.g__Vitis_3603.s__Vitis.vinifera_29760", ]
Reynosia.septentrionalis_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Rosales_3744.f__Rhamnaceae_3608.g__Reynosia_106701.s__Reynosia.septentrionalis_1226765", ]
Croton.linearis_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Malpighiales_3646.f__Euphorbiaceae_3977.g__Croton_100370.s__Croton.linearis_1003700", ]
Lactuca.sativa_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Asterales_4209.f__Asteraceae_4210.g__Lactuca_4235.s__Lactuca.sativa_4236", ]
Thalassia_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__Liliopsida_4447.o__Alismatales_16360.f__Hydrocharitaceae_26319.g__Thalassia_55495.__", ]
Conocarpus.sericeus_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Myrtales_41944.f__Combretaceae_3954.g__Conocarpus_134916.s__Conocarpus.sericeus_578696", ]
Cucumis.melo_data <- melted_data[melted_data$OTU == "Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Cucurbitales_71239.f__Cucurbitaceae_3650.g__Cucumis_3655.s__Cucumis.melo_3656", ]

#perform wilcox signed-rank test for top taxa abundance between analysis methods
Casasia_metabarcoding <- Casasia.clusiifolia_data$Abundance[Casasia.clusiifolia_data$Analysis == "Metabarcoding"] 
Casasia_morphology <- Casasia.clusiifolia_data$Abundance[Casasia.clusiifolia_data$Analysis == "Morphology"]
wilcox_casasia_paired <- wilcox.test(Casasia_metabarcoding, Casasia_morphology, paired = TRUE, exact = FALSE)

Rachicallis_metabarcoding <- Rachicallis.americana_data$Abundance[Rachicallis.americana_data$Analysis == "Metabarcoding"] 
Rachicallis_morphology <- Rachicallis.americana_data$Abundance[Rachicallis.americana_data$Analysis == "Morphology"]
wilcox_rachicallis_paired <- wilcox.test(Rachicallis_metabarcoding, Rachicallis_morphology, paired = TRUE, exact = FALSE)

Manilkara_metabarcoding <- Manilkara.bahamensis_data$Abundance[Manilkara.bahamensis_data$Analysis == "Metabarcoding"] 
Manilkara_morphology <- Manilkara.bahamensis_data$Abundance[Manilkara.bahamensis_data$Analysis == "Morphology"]
wilcox_manilkara_paired <- wilcox.test(Manilkara_metabarcoding, Manilkara_morphology, paired = TRUE, exact = FALSE)

Guaiacum_metabarcoding <- Guaiacum.sanctum_data$Abundance[Guaiacum.sanctum_data$Analysis == "Metabarcoding"] 
Guaiacum_morphology <- Guaiacum.sanctum_data$Abundance[Guaiacum.sanctum_data$Analysis == "Morphology"]
wilcox_guaiacum_paired <- wilcox.test(Guaiacum_metabarcoding, Guaiacum_morphology, paired = TRUE, exact = FALSE)

Guapira_metabarcoding <- Guapira.discolor_data$Abundance[Guapira.discolor_data$Analysis == "Metabarcoding"] 
Guapira_morphology <- Guapira.discolor_data$Abundance[Guapira.discolor_data$Analysis == "Morphology"]
wilcox_guapira_paired <- wilcox.test(Guapira_metabarcoding, Guapira_morphology, paired = TRUE, exact = FALSE)

Vitis_metabarcoding <- Vitis.vinifera_data$Abundance[Vitis.vinifera_data$Analysis == "Metabarcoding"] 
Vitis_morphology <- Vitis.vinifera_data$Abundance[Vitis.vinifera_data$Analysis == "Morphology"]
wilcox_vitis_paired <- wilcox.test(Vitis_metabarcoding, Vitis_morphology, paired = TRUE, exact = FALSE)

Reynosia_metabarcoding <- Reynosia.septentrionalis_data$Abundance[Reynosia.septentrionalis_data$Analysis == "Metabarcoding"] 
Reynosia_morphology <- Reynosia.septentrionalis_data$Abundance[Reynosia.septentrionalis_data$Analysis == "Morphology"]
wilcox_reynosia_paired <- wilcox.test(Reynosia_metabarcoding, Reynosia_morphology, paired = TRUE, exact = FALSE)

Croton_metabarcoding <- Croton.linearis_data$Abundance[Croton.linearis_data$Analysis == "Metabarcoding"] 
Croton_morphology <- Croton.linearis_data$Abundance[Croton.linearis_data$Analysis == "Morphology"]
wilcox_croton_paired <- wilcox.test(Croton_metabarcoding, Croton_morphology, paired = TRUE, exact = FALSE)

Lactuca_metabarcoding <- Lactuca.sativa_data$Abundance[Lactuca.sativa_data$Analysis == "Metabarcoding"] 
Lactuca_morphology <- Lactuca.sativa_data$Abundance[Lactuca.sativa_data$Analysis == "Morphology"]
wilcox_lactuca_paired <- wilcox.test(Lactuca_metabarcoding, Lactuca_morphology, paired = TRUE, exact = FALSE)

Thalassia_metabarcoding <- Thalassia_data$Abundance[Thalassia_data$Analysis == "Metabarcoding"] 
Thalassia_morphology <- Thalassia_data$Abundance[Thalassia_data$Analysis == "Morphology"]
wilcox_thalassia_paired <- wilcox.test(Thalassia_metabarcoding, Thalassia_morphology, paired = TRUE, exact = FALSE)

Conocarpus_metabarcoding <- Conocarpus.sericeus_data$Abundance[Conocarpus.sericeus_data$Analysis == "Metabarcoding"] 
Conocarpus_morphology <- Conocarpus.sericeus_data$Abundance[Conocarpus.sericeus_data$Analysis == "Morphology"]
wilcox_conocarpus_paired <- wilcox.test(Conocarpus_metabarcoding, Conocarpus_morphology, paired = TRUE, exact = FALSE)

Cucumis_metabarcoding <- Cucumis.melo_data$Abundance[Cucumis.melo_data$Analysis == "Metabarcoding"] 
Cucumis_morphology <- Cucumis.melo_data$Abundance[Cucumis.melo_data$Analysis == "Morphology"]
wilcox_cucumis_paired <- wilcox.test(Cucumis_metabarcoding, Cucumis_morphology, paired = TRUE, exact = FALSE)


#perform wilcox rank sum test for top taxa between analysis methods
wilcox_casasia <- wilcox.test(Abundance ~ Analysis, data = Casasia.clusiifolia_data)
wilcox_rachicallis <- wilcox.test(Abundance ~ Analysis, data = Rachicallis.americana_data)
wilcox_manilkara <- wilcox.test(Abundance ~ Analysis, data = Manilkara.bahamensis_data)
wilcox_guaiacum <- wilcox.test(Abundance ~ Analysis, data = Guaiacum.sanctum_data)
wilcox_guapira <- wilcox.test(Abundance ~ Analysis, data = Guapira.discolor_data)
wilcox_vitis <- wilcox.test(Abundance ~ Analysis, data = Vitis.vinifera_data)
wilcox_reynosia <- wilcox.test(Abundance ~ Analysis, data = Reynosia.septentrionalis_data)
wilcox_croton <- wilcox.test(Abundance ~ Analysis, data = Croton.linearis_data)
wilcox_lactuca <- wilcox.test(Abundance ~ Analysis, data = Lactuca.sativa_data)
wilcox_vitis <- wilcox.test(Abundance ~ Analysis, data = Vitis.vinifera_data)
wilcox_thalassia <- wilcox.test(Abundance ~ Analysis, data = Thalassia_data)
wilcox_conocarpus <- wilcox.test(Abundance ~ Analysis, data = Conocarpus.sericeus_data)
wilcox_cucumis <- wilcox.test(Abundance ~ Analysis, data = Cucumis.melo_data)

#save raw p values for paired signed rank test to list
#make empty list for p values first
p_values <- list()
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Gentianales_4055.f__Rubiaceae_24966.g__Casasia_58499.s__Casasia.clusiifolia_58500"]] <- wilcox_casasia_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Gentianales_4055.f__Rubiaceae_24966.g__Rachicallis_171717.s__Rachicallis.americana_1717177"]] <- wilcox_rachicallis_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Ericales_41945.f__Sapotaceae_3737.g__Manilkara_3740.s__Manilkara.bahamensis_37400"]] <- wilcox_manilkara_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Zygophyllales_403666.f__Zygophyllaceae_43873.g__Guaiacum_45188.s__Guaiacum.sanctum_45189"]] <- wilcox_guaiacum_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__undef__0.o__Caryophyllales_3524.f__Nyctaginaceae_3536.g__Guapira_427802.s__Guapira.discolor_427803"]] <- wilcox_guapira_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Vitales_403667.f__Vitaceae_3602.g__Vitis_3603.s__Vitis.vinifera_29760"]] <- wilcox_vitis_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Rosales_3744.f__Rhamnaceae_3608.g__Reynosia_106701.s__Reynosia.septentrionalis_1226765"]] <- wilcox_reynosia_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Malpighiales_3646.f__Euphorbiaceae_3977.g__Croton_100370.s__Croton.linearis_1003700"]] <- wilcox_croton_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__asterids_71274.o__Asterales_4209.f__Asteraceae_4210.g__Lactuca_4235.s__Lactuca.sativa_4236"]] <- wilcox_lactuca_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__Liliopsida_4447.o__Alismatales_16360.f__Hydrocharitaceae_26319.g__Thalassia_55495.__"]] <- wilcox_thalassia_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Myrtales_41944.f__Combretaceae_3954.g__Conocarpus_134916.s__Conocarpus.sericeus_578696"]] <- wilcox_conocarpus_paired$p.value
p_values[["Root.k__Viridiplantae_33090.p__Streptophyta_35493.c__sub__rosids_71275.o__Cucurbitales_71239.f__Cucurbitaceae_3650.g__Cucumis_3655.s__Cucumis.melo_3656"]] <- wilcox_cucumis_paired$p.value

#convert to data frame
wilcox_paired_top_taxa_df <- data.frame(Taxon = names(p_values),p_raw = unlist(p_values))


#adjust p-values for false discovery rate with Benjamini-Hochberg (FDR)
wilcox_paired_top_taxa_df$p_adjusted <- p.adjust(wilcox_paired_top_taxa_df$p_raw, method = "fdr")
#export summary as csv file
wilcox_paired_top_taxa_summary <- summary(wilcox_paired_top_taxa_df)
wilcox_paired_top_taxa_summary
write.csv(wilcox_paired_top_taxa_summary, "wilcox_paired_top_taxa_summary.csv", row.names = FALSE)
#export data frame as csv.
write.csv(wilcox_paired_top_taxa_df, "wilcox_paired_top_taxa_df.csv", row.names = FALSE)


#filter significant taxa by ajusted p-values
significant_top_taxa_paired <- wilcox_paired_top_taxa_df %>% filter(p_adjusted < 0.05) %>% arrange(p_adjusted)
print(significant_top_taxa_paired)

################################## PLOT TOP 12 TAXA AS BARCHART #####################################

#load packages
library(ggplot2)
library(tidyr)

#import excel data top taxa list for each method
iguana_diet_top_taxa <- read.delim("~/Desktop/Iguana_Diet/Iguana_diet_molec_visual_paired/iguana_diet_top_taxa_observed_paired.txt")

#sort by descending order not alphabetical
iguana_diet_top_taxa$Species<-factor(iguana_diet_top_taxa$Species, levels = c(iguana_diet_top_taxa$Species[c(1:12)]))

#pivot to create method column
iguana_diet_top_taxa<-pivot_longer(iguana_diet_top_taxa, cols = c("Metabarcoding", "Morphology"), names_to = "Analysis", values_to = "Observed")

#make barchart
iguana_diet_barchart<-ggplot()
iguana_diet_barchart<-iguana_diet_barchart + geom_col(data = iguana_diet_top_taxa, aes(x=Species, y=Observed, fill=Analysis), position = "dodge") + 
  labs(y="Observed", x="") + theme_bw() + scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(face = "italic")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18)) + theme(axis.title.y = element_text(size = 18),legend.title = element_text(size = 18))
iguana_diet_barchart

iguana_diet_barchart<-ggplot()
iguana_diet_barchart<-iguana_diet_barchart + geom_col(data = iguana_diet_top_taxa, aes(x=Species, y=Observed, fill=Analysis), position = "dodge") + 
  labs(y="Observed", x="") + theme_bw() + scale_fill_brewer(palette = "Dark2") + theme(axis.text.x = element_text(face = "italic")) + theme(text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
iguana_diet_barchart
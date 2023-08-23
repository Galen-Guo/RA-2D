# Load and prep files ####
###load libraries ####
library(tidyverse)
library(vegan)
###load output files ####

####load taxa files ####
genome_out_raw <- read_tsv(file = "MetaLabMAG_output/taxonomy_analysis/Genome.tsv") %>%
  type.convert(as.is = TRUE)

taxa_out_raw <- read_tsv(file = "MetaLabMAG_output/taxonomy_analysis/Taxa.tsv") %>%
  type.convert(as.is = TRUE) 
colnames(taxa_out_raw) <- str_replace(names(taxa_out_raw), "Intensity Guo_", "Guo_")

####load function files ####
func_out_raw <- read_tsv(file = "MetaLabMAG_output/functions.tsv") %>%
  type.convert(as.is = TRUE) 
colnames(func_out_raw) <- str_replace(names(func_out_raw), "Intensity Guo_", "Guo_")

####load metadata files ####
metadata <- read.delim2(file = "MetaLabMAG_output/metainfo.tsv") %>%
  mutate(Raw = str_replace(Raw, "Intensity Guo_", "Guo_"))

###### Change column name to whatever fit the metadata
colnames(metadata) <- c("Raw", # keep it as raw to reconcile with the rest of the code
                        "ExperimentName", # this one too
                        "Treatment",
                        "Dilution",
                        "Centrifuge")

### work with functional table with protein annotation group ####
func_out <- func_out_raw
colname_func_out <- as.data.frame(colnames(func_out))
colnames(colname_func_out) <- c("Raw") 
colname_func_out <- colname_func_out %>%
  left_join(metadata, by = "Raw") 
col <- colname_func_out$Raw
colnames(func_out) <-col

##### colnames of extra info and category to choose from ####
# [183] "Protein name", [184] "Seed ortholog", [185] "Description" ,[186] "Taxonomy Id"                                          
# [187] "Taxonomy name"[188] "Preferred name", 
# [189] "Gene_Ontology_id", [190] "Gene_Ontology_name" [191] "Gene_Ontology_namespace",
# [192] "EC_id", [193] "EC_de"[194] "EC_an", [195] "EC_ca"[196] 
# "KEGG_ko" [197] "KEGG_Pathway_Entry" [198] "KEGG_Pathway_Name"[199] "KEGG_Module"                                          
# [200] "KEGG_Reaction" [201] "KEGG_rclass" [202] "BRITE"[203] "KEGG_TC" 
# [204] "CAZy" [205] "BiGG_Reaction"[206] "PFAMs"                                                
# [207] "COG accession" [208] "COG category" [209] "COG name"                                             
# [210] "NOG accession"[211] "NOG category"[212] "NOG name"  

### combine functional table with metadata ####
func_metadata <- func_out %>%
  select(Name, contains("Guo"))%>%
  column_to_rownames("Name") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Raw") %>%
  inner_join(metadata)

### substract treatment intensity by control intensity ####

# ##### Calculate average control intensity  ####
# func_water_only <- func_metadata %>%
#   filter(Treatment == "ddH2O") %>%
#   group_by(Individual) %>%
#   summarise(across(where(is.numeric), ~ mean(.))) %>%
#   column_to_rownames("Individual") %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column("protein") 
# 
# func_DMSO_only <- func_metadata %>%
#   filter(Treatment == "DMSO") %>%
#   group_by(Individual) %>%
#   summarise(across(where(is.numeric), ~ mean(.))) %>%
#   column_to_rownames("Individual") %>%
#   t() %>% as.data.frame()%>%
#   rownames_to_column("protein") 
# 
# # rbind(func_water_only, func_DMSO_only) %>%
# #   write.table("func_ctrl.txt")
# 
# #### substract intensity  ####
# 
# substr_func_water <-  func_metadata %>%
#   column_to_rownames("Label") %>%
#   filter(Treatment != "Erythrosine_B_Free_Acid") %>%
#   select(-c(Treatment, Concentration, ExperimentName, Volunteer)) %>%
#   t() %>% as.data.frame() %>%
#   mutate(across(.cols = everything(), ~ . - func_water_only)) 
# 
# substr_func_water_V47 <-  func_metadata %>%
#   column_to_rownames("Label") %>%
#   filter(Treatment != "Erythrosine_B_Free_Acid") %>%
#   select(-c(Treatment, Concentration, ExperimentName, Volunteer)) %>%
#   t() %>% as.data.frame() %>%
#   mutate(across(.cols = everything(), ~ . - func_water_only$V47)) 
# 
# substr_func_DMSO_IBD <-  func_metadata %>%
#   column_to_rownames("Raw") %>%
#   filter(Individual == "IBD1") %>%
#   filter(Treatment == "Erythrosine_B_Free_Acid") %>%
#   select(-c(Treatment, Concentration, Individual, Replicate)) %>%
#   t() %>% as.data.frame() %>%
#   mutate(across(.cols = everything(), ~ . - func_water_only$IBD1))  
# 
# substr_func_DMSO_V47 <-  func_metadata %>%
#   column_to_rownames("Raw") %>%
#   filter(Individual == "V47") %>%
#   filter(Treatment == "Erythrosine_B_Free_Acid") %>%
#   select(-c(Treatment, Concentration, Individual, Replicate)) %>%
#   t() %>% as.data.frame() %>%
#   mutate(across(.cols = everything(), ~ . - func_water_only$V47)) 

### add column that contains protein group info ####
## prep table
func_group_info <- func_out %>%
  select(-c(contains("Guo"), Group_ID, Protein_ID))
# 
# substr_func_data_raw <- cbind(substr_func_water_IBD, substr_func_DMSO_IBD, 
#                           substr_func_water_V47, substr_func_DMSO_V47) %>%
#   rownames_to_column("Name") %>%
#   left_join(func_group_info, by = join_by(Name))
# 
# func_out = substr_func_data_raw
### separate table into protein CAT (KEGG vs COG, NOG.. etc) ####
for (cat in c("COG", "KEGG", "NOG")) {
  obj_name <- paste(cat,"_table", sep = "")
  assign(obj_name, 
         func_out %>%
           select(Name, contains("Guo"), starts_with(cat)))
  done
}





# Create hetmap with each COG category ####

### combine cog C with metadata ####
COG_C_metadata <- COG_table %>%
  filter(`COG category` == "C") %>% 
  select(`COG accession`, contains("Guo")) %>%
  # mutate(`COG accession` = if_else(`COG accession` == "NA", "unknown", `COG accession`)) %>%
  group_by(`COG accession`) %>%
  summarise(across(where(is.numeric), ~ mean(.))) %>%
  ungroup() %>%
  column_to_rownames("COG accession") %>%
  t() %>%  as.data.frame() %>%
  rownames_to_column("Raw") %>%
  inner_join(metadata) 


### if you want to filter by whatever factor change the appropriate value ####

COG_C_metadata_input <- COG_C_metadata %>%
  filter(Treatment == "KES")
### if you want to filter keep default table  ####

COG_C_metadata_input <- COG_C_metadata

### create annotation color table ####

annotation_row <- COG_C_metadata_input %>%
  column_to_rownames("Raw") %>%
  # select(-c(starts_with("COG"), Replicate)) %>%
  dplyr::select(Dilution, Treatment, Centrifuge)

### create input table for pheatmap annotation ####


pheatmap_input <- COG_C_metadata_input %>%
  select(Raw, starts_with("COG")) %>%
  column_to_rownames("Raw") %>%
  select(-which(colSums(.,na.rm = TRUE) < 100000))

drows <- vegan::vegdist(pheatmap_input,
                        method = "robust.aitchison")
dcols <- vegan::vegdist(t(pheatmap_input),
                        method = "robust.aitchison")

pheatmap::pheatmap(pheatmap_input,
                   scale = "row",
                   annotation_row = annotation_row,
                   # clustering_distance_rows = drows, 
                   clustering_distance_cols = dcols,
                   show_rownames = FALSE,
                   main = "clustering by row and column using Robust Aitchison distance, scaled by row")



# Alpha diversity ####
shanTD_map <- func_out %>%
  select(contains("Guo"))%>%
  # column_to_rownames("Name") %>%
  vegan::decostand(method ="hellinger") %>%
  t() %>% as.data.frame() %>%
  mutate(Shannon = vegan::diversity(., "shan"),
         ShanTD = Shannon/log(vegan::specnumber(.))) %>%
  rownames_to_column("Raw") %>%
  inner_join(metadata) %>%
  select(Shannon, ShanTD, colnames(metadata))

shanTD_map %>%
  # filter(Treatment != "EMPTY",
  #        Treatment != "REF") %>%
  ggplot(aes(x = Centrifuge, y = Shannon, fill = Dilution,
             color = Centrifuge)) +
  # geom_line(size = 1.5, show.legend = FALSE) +
  geom_boxplot(color = "black", show.legend = FALSE, alpha = .5) +
  geom_point(shape = 21, size = 5, color = "black" , 
             position = position_dodge(width = .75)) +
  facet_grid(
    cols = vars(Treatment),
    scales = "free", space = "free") +
  theme_bw() + 
  theme(line = element_line(size = 1), 
        text = element_text(size = 18),
        # axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(color = "black"),
        # axis.title.y = element_blank(),
        legend.text = element_text(size =10)) +
  labs(xlab = "Treatment", ylab = "ShanTD") +
  labs(title = "Hellinger transformed-functional quantification using TMT-labelled peptides")

# ordination ####

### data transformation using robust centre-log ratio ####

#### prep table to have unique identifier as rownames and no chr in the table. (remove excess info) ####
func_table <- COG_table %>%
  select(-starts_with("COG")) %>%
  column_to_rownames("Name") 

#### transform data using rclr (assuming compositional dataset) ####
func_table_transformed <- func_table %>%
  vegan::decostand(method ="hellinger") %>%
  
  na.omit()

### PCA ####
mds_func_raw <- vegan::pcnm(vegan::vegdist(func_table))
PCA_func_raw <- princomp(func_table)

PCA_func <- metadata[metadata$Raw %in% colnames(func_table_transformed),] %>%
  mutate(PCA1 = vegan::scores(PCA_func_raw,display =  "species")[,c(1)]) %>%
  mutate(PCA2 = vegan::scores(PCA_func_raw,display =  "species")[,c(2)]) %>%
  # mutate(PCA1 = PCA_func_raw$loadings[,1]) %>%
  # mutate(PCA2 = PCA_func_raw$scores[,2]) %>%
  filter(Treatment != "EMPTY") %>%
  ggplot(aes(
    x = PCA1,
    y = PCA2,
    shape = Treatment,
    color = Dilution)) +
  geom_point(
    # shape = 21,
    aes(shape = Treatment),
    show.legend = TRUE
  ) +
  # ggrepel::geom_label_repel(aes( 
  #   x = PCA1,
  #   y = PCA2,
  #   label =  Concentration)) + #to check if labels are correct
  # scale_fill_manual(
  #   # values = c("#a6cee3", #high
  #   #                            # "#1f78b4", 
  #   #                            # "#b2df8a", #med
  #   #                            "#33a02c" #low
  #   #                            ),
  #                   # labels = c("HIGH", 
#                   #            "MED", 
#                   #            "LOW"
#                   #            )
#                   )+
# scale_shape_manual(values = c(21,
#                               22,
#                               23,
#                               24,
#                               25),
#                    # labels = c("IBD1", 
#                    #            "V47")
#                               ) +
labs(x = "PCA1 (25.8%)",
     y = "PCA3 (14.8%)",
) +
facet_wrap(vars(Treatment),
           ncol = 3)
# facet_grid(cols = vars(meta3),
#            rows = vars(meta1)) 

PCA_func


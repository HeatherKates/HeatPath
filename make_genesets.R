library(msigdbr)
library(dplyr)

dir.create("data", showWarnings = FALSE)

organisms <- c("Homo sapiens" = "human", "Mus musculus" = "mouse")

for (org in names(organisms)) {
  tag <- organisms[[org]]
  message("Processing gene sets for: ", tag)
  
  gene_sets <- list()
  
  # ---- GO (BP, MF, CC) ----
  go_df <- msigdbr(species = org, category = "C5")
  go_list <- split(go_df$gene_symbol, go_df$gs_name)
  go_list <- go_list[lengths(go_list) >= 10]
  gene_sets$GO <- go_list
  
  # ---- KEGG ----
  kegg_df <- msigdbr(species = org, category = "C2") %>%
    filter(gs_subcat %in% c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"))
  
  if (nrow(kegg_df) > 0) {
    kegg_list <- split(kegg_df$gene_symbol, kegg_df$gs_name)
    kegg_list <- kegg_list[lengths(kegg_list) >= 10]
    gene_sets$KEGG <- kegg_list
  }
  
  # ---- Reactome ----
  reactome_df <- msigdbr(species = org, category = "C2") %>%
    filter(gs_subcat == "CP:REACTOME")
  if (nrow(reactome_df) > 0) {
    reactome_list <- split(reactome_df$gene_symbol, reactome_df$gs_name)
    reactome_list <- reactome_list[lengths(reactome_list) >= 10]
    gene_sets$Reactome <- reactome_list
  }
  
  # ---- BioCarta ----
  biocarta_df <- msigdbr(species = org, category = "C2") %>%
    filter(gs_subcat == "CP:BIOCARTA")
  if (nrow(biocarta_df) > 0) {
    biocarta_list <- split(biocarta_df$gene_symbol, biocarta_df$gs_name)
    biocarta_list <- biocarta_list[lengths(biocarta_list) >= 10]
    gene_sets$BioCarta <- biocarta_list
  }
  
  # ---- Save ----
  saveRDS(gene_sets, file = paste0("data/gene_sets_", tag, ".rds"))
}



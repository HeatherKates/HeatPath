library(msigdbr)
library(dplyr)
library(stringr)

dir.create("data", showWarnings = FALSE)

organisms <- c("Homo sapiens" = "human", "Mus musculus" = "mouse")

for (org in names(organisms)) {
  tag <- organisms[[org]]
  message("Processing gene sets for: ", tag)
  
  gene_sets <- list()
  
  # Helper function to build a gene set list and metadata
  build_gene_set <- function(df, name = "unknown") {
    gene_list <- split(df$gene_symbol, df$gs_name)
    gene_list <- gene_list[lengths(gene_list) >= 10]
    
    meta <- df %>%
      distinct(gs_name, gs_description) %>%
      filter(gs_name %in% names(gene_list)) %>%
      mutate(
        set_id = gs_name,
        name_clean = str_replace_all(gs_name, "_", " "),
        searchable_text = paste(gs_name, gs_description, sep = " | ")
      )
    
    return(list(
      genes = gene_list,
      metadata = meta
    ))
  }
  
  # ---- GO (BP, MF, CC) ----
  go_df <- msigdbr(species = org, category = "C5")
  gene_sets$GO <- build_gene_set(go_df, "GO")
  
  # ---- KEGG ----
  kegg_df <- msigdbr(species = org, category = "C2") %>%
    filter(gs_subcat %in% c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"))
  if (nrow(kegg_df) > 0) {
    gene_sets$KEGG <- build_gene_set(kegg_df, "KEGG")
  }
  
  # ---- Reactome ----
  reactome_df <- msigdbr(species = org, category = "C2") %>%
    filter(gs_subcat == "CP:REACTOME")
  if (nrow(reactome_df) > 0) {
    gene_sets$Reactome <- build_gene_set(reactome_df, "Reactome")
  }
  
  # ---- BioCarta ----
  biocarta_df <- msigdbr(species = org, category = "C2") %>%
    filter(gs_subcat == "CP:BIOCARTA")
  if (nrow(biocarta_df) > 0) {
    gene_sets$BioCarta <- build_gene_set(biocarta_df, "BioCarta")
  }
  
  # ---- Save ----
  saveRDS(gene_sets, file = paste0("data/gene_sets_", tag, ".rds"))
}


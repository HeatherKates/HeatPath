library(msigdbr)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(readr)

dir.create("data", showWarnings = FALSE)

org <- "Homo sapiens"

# Step 1: Lookup table for GO term → GO:ID
go_term_lookup <- AnnotationDbi::select(GO.db, keys = keys(GO.db),
                                        columns = c("GOID", "TERM"), keytype = "GOID") %>%
  mutate(term = tolower(TERM)) %>%
  dplyr::select(go_id = GOID, term)

# Step 2: Retrieve MSigDB gene sets
databases <- list(
  GO = msigdbr(species = org, category = "C5"),
  KEGG = msigdbr(species = org, category = "C2") %>% filter(gs_subcat %in% c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS")),
  Reactome = msigdbr(species = org, category = "C2") %>% filter(gs_subcat == "CP:REACTOME"),
  BioCarta = msigdbr(species = org, category = "C2") %>% filter(gs_subcat == "CP:BIOCARTA")
)

# Step 3: Mapping function for GO
map_gs_name_to_go <- function(gs_name_vector) {
  cleaned <- gs_name_vector %>%
    str_remove("^GOBP_|^GOMF_|^GOCC_") %>%
    str_replace_all("_", " ") %>%
    tolower()
  
  matched_go_ids <- sapply(cleaned, function(term) {
    idx <- match(term, go_term_lookup$term)
    if (!is.na(idx)) go_term_lookup$go_id[idx] else NA
  })
  
  matched_go_ids
}

# Step 4: ID extractors for KEGG, Reactome
extract_external_id <- function(df, db) {
  case_when(
    db == "KEGG"     ~ str_extract(df$gs_url, "hsa\\d{5}"),
    db == "Reactome" ~ str_extract(df$gs_url, "R-HSA-\\d+"),
    db == "BioCarta" ~ NA_character_,
    TRUE             ~ NA_character_
  )
}

# Step 5: Build and bind all
converter_list <- lapply(names(databases), function(db) {
  df <- databases[[db]] %>%
    distinct(msigdbr_id = gs_id, gs_name, gs_description, gs_url)
  
  if (db == "GO") {
    df <- df %>%
      mutate(
        external_id = map_gs_name_to_go(gs_name),
        db = db
      )
  } else {
    df <- df %>%
      mutate(
        external_id = extract_external_id(., db),
        db = db
      )
  }
  
  df %>% filter(!is.na(external_id))
})

# Step 6: Combine and write
converter_df <- bind_rows(converter_list) %>%
  distinct(external_id, msigdbr_id, db, gs_name, gs_description)

write_tsv(converter_df, "data/msigdbr_id_converter.txt")

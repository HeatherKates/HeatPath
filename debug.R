GO <- readRDS("data/gene_sets_human.rds")

logcpm_data <- "test_data/Sharma_LCPM.csv"
  df <- read.csv(logcpm_data)
  gene_ids <- df[[1]]
  keep <- !is.na(gene_ids) & gene_ids != ""
  df <- df[keep, , drop = FALSE]
  rownames(df) <- make.unique(df[[1]])
  df[[1]] <- NULL
  logcpm_data= as.matrix(df)
  
  
sample_info <- "test_data/Sharma_sampleinfo.csv"
    
    df <- read.csv(sample_info)
    colnames(df)[1:2] <- c("sample_name", "group")
    df <- df %>% filter(sample_name %in% colnames(logcpm_data))
    sample_info=df

    # Pick the same DB and pathway you'd select in the app
    db <- "GO"
    pathway <- "GOBP_2FE_2S_CLUSTER_ASSEMBLY"
    
    genes_in_pathway <- GO[[db]][[pathway]]
    length(genes_in_pathway)  # how many genes are in the pathway
    
    matched_genes <- intersect(rownames(logcpm_data), genes_in_pathway)
    length(matched_genes)      # how many matched
    matched_genes              # which ones matched
    
# Graphics and analysis selection fungi BLW2
# Author: Anja Logo
# Started: 16.09.24
# Last changes: 16.01.24

# Loading environment
source(file = "setup_F24.R")

# ATTENTION load one or the other
path <- "10000_reads/"
#path <- "20000_reads/"
output_path <- "../../../04_Dissemination/01_Publikationen/figures_pdf/"

# Functions-----

# Functions from external file
source("../../20230227_Bioassays_funct.R") # Loading Functions

# Function: Transform characters into single letters with [] and seperated by |
transform_characters <- function(text) {
  unique_chars <- unique(strsplit(text, "")[[1]])
  transformed_text <- paste0("[", paste(unique_chars, collapse = "]|["), "]")
  return(transformed_text)
}

# Function: boxplot with asterix
plot_box_ast = function(data_sum_all, data.test, color.code, variable, y.text) {
  
  # Create column with stars for composts that are significantly different from the peat substrate
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, "*", ""))
  
  data.test$variable = variable
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  labelling <- list(
    '1'="May 2022",
    '2'="July 2022",
    '3'="September 2022",
    '4'="May 2023"
  )
  batch_labeller <- function(variable,value){
    return(labelling[value])
  }
  plot = data.test %>%
    ggplot( aes(x = treatment, y = variable)) +
    geom_boxplot(position = position_dodge(width = 0.9), width=0.8, fill = color.code, outlier.size = 0.8) +
    background_pr+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, size =10), axis.title = element_text(size =12),
          legend.position = "none", ) +labs ( x = "", y = y.text )+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    facet_grid(.~batch, scales= "free", switch = "x", space ="free_x", labeller = batch_labeller)+
    theme(strip.placement = "outside", strip.text.x = element_text(size=11))+
    geom_text(data = value_max, aes(x=treatment, y = 1 + quantile, label = data_sum_all$sig), vjust=0, size= 5, colour ="red")
}

# Scatter plot
plot_scatter <-function(data, variable1, variable2, method = "spearman", xaxis = "Relative biomass [%]", yaxis = "Rating [5-25]"){
  data$x <-data[,variable1];data$y <- data[, variable2]
  cor_res  <- cor.test(data$x, data$y, method = method)
  rho <- round(cor_res$estimate, 2)
  p_val <- ifelse(cor_res$p.value < 0.001, "< 0.001", paste("=", round(cor_res$p.value,3)))
  plot <- ggplot(data = data, aes(x = x, y = y)) +
    geom_point() +
    geom_text(x = min(data$x), y = max(data$y), 
              label = paste0("rho = ", rho, ", p ", p_val),
              hjust = 0, vjust = 1, size = 4) +
    # Label the axes with the variable names
    labs(x = xaxis, y = yaxis) +
    theme_bw()
  return(plot)
}

# Function: select batch concentration
batch.select.fun = function(data, batch.select){
  data %>%
    filter(paste(batch, conc) %in% paste(batch.select$batch, batch.select$conc))
}

sig.differences = function(x){
  significance <- ifelse(x < 0.001, "***",
                         ifelse(x < 0.01, "**",
                                ifelse(x < 0.05, "*",
                                       ifelse(x < 0.1, ".", " "))))
  print(significance)
}

# Discard missing ASVs in a dataframe
discard_missing_asvs <- function(asv_table) {
  present_columns <- apply(asv_table, 2, function(col) any(col > 0))
  return(asv_table[, present_columns])
}

# Random color generator
generate_random_color <- function() {
  rgb(runif(1), runif(1), runif(1))
} # Function to generate a random RGB color

# data.asv|ISS: ASV/genera table with in rows and samples in columns
# top & flop: meta_file for the selected composts
# path: location where the output files should be saved
# ncompres: number of composts the ASV/genera must be present, default 6
ISA_PALMA = function(top, flop, data.ISS, data.asv, path, fct ="all", path.plant, ncompres = 6, nperm =9999, cutoff = TRUE, correction = FALSE){
  
  # Check if 'method' is one of the supported methods
  if (!(fct %in% c("r.g", "IndVal.g", "aldex", "ancom", "maaslin", "all"))) {
    stop("Invalid method. Supported methods: 'r.g', 'IndVal.g', 'aldex', 'ancom', 'maaslin', 'all'")
  }
  
  # Check if 'path.plant is in the supported systems
  if (!(path.plant %in% c("cress", "cugu", "curs", "4v4", "14v6"))) {
    stop("Invalid method. Supported systems: 'cress', 'cugu', 'curs', '4v4', '14vs6' ")
  }
  
  # Subset data to only get ASV/genera table
  top$group = "top" ; flop$group = "flop"
  selected_composts <- rbind(top, flop)
  
  data = data.ISS[, selected_composts$treatment]
  if (cutoff == FALSE) {
    # rarefied
    data_01 <-ifelse(data > 0,1,0) %>% rowSums() %>% as.data.frame()
  } else {
    data_01 <-ifelse(data >= 1,1,0) %>% rowSums() %>% as.data.frame()
  }
  data.ISS.red = data[data_01 > (ncompres-1),]
  rm(data, data_01)
  # unrarefied
  data.asv.red = data.asv[rownames(data.ISS.red),]
  list[["nASV"]] = nrow(data.asv.red)
  #list[["nBLW1"]] =  rownames(data.asv.red) %in% asv_names_BLW1_99 %>% sum()
  
  # Point Biserial correlation based on rarefied data
  if(fct %in% c("r.g", "all")){
    # ASVs in columns and samples in row needed
    res.assoc <- multipatt(x=t(data.ISS.red), cluster=selected_composts$group, func= "r.g",
                           control = how(nperm=nperm), duleg=T, print.perm = T)
    res <-res.assoc$sign # save all the outputs/comparisons
    if (correction == TRUE) {
      res$qval = p.adjust(res$p.value, method ="fdr")
      res.ind = res[res$s.top ==1 & res$qval < 0.1,]
    } else{
      res.ind = res[res$s.top ==1 & res$p.value < 0.05,]}
    PBC = na.omit(res.ind) %>% dplyr::arrange(-stat)
    write.csv(PBC, file = paste0(path, path.plant, ifelse(cutoff == TRUE, "_cut",""), "_PBC.csv"))
    print("PBC DONE")
  }
  
  # Indicator species analysis
  if(fct %in% c("IndVal.g", "all")){
    # ASVs in columns and samples in row needed
    res.assoc <- multipatt(x=t(data.ISS.red), cluster=selected_composts$group, func= "IndVal.g",
                           control = how(nperm=nperm), duleg=T, print.perm = T)
    res <-res.assoc$sign # save all the outputs/comparisons
    if (correction == TRUE) {
      res$qval = p.adjust(res$p.value, method ="fdr")
      res.ind = res[res$s.top ==1 & res$qval < 0.1,]
    } else{
      res.ind = res[res$s.top ==1 & res$p.value < 0.05,]}
    ISA = na.omit(res.ind) %>% dplyr::arrange(-stat)
    write.csv(ISA, file = paste0(path ,path.plant, ifelse(cutoff == TRUE, "_cut",""), "_ISA.csv"))
    print("ISA DONE")
  }
  
  # dissimiliarity based MANOVA
  
  #if(fct %in% c("dbMANOVA", "all")){
  
  #manova.indic <- dbMANOVAspecies(t(data.ISS.red), selected_composts$group, nrep = 999)
  #results.manova <- as.data.frame(cbind(manova.indic$test$names, manova.indic$test$pvalue))
  #colnames(results.manova) = c("otu", "pvalue")
  #res <- results.manova %>% dplyr::filter(pvalue < 0.05) # Only select the ones which are enriched in the suppressive composts
  
  # Select only the genera which are enriched in the top composts
  #selected_composts = as.data.frame(selected_composts)
  #rownames(selected_composts) = selected_composts$treatment
  
  #res = res[!res$otu == "GLOBAL",]
  #data_merged =merge(selected_composts[,c("treatment", "group")], t(data.ISS.red[res$otu,]), by =0)
  #rownames(data_merged) <- data_merged$Row.names
  #data_merged$Row.names <- NULL; data_merged$treatment <-NULL
  
  # Calculate column sums for top and flop separately
  #top_sums <- colSums(data_merged[data_merged$group =="top", 2:ncol(data_merged)])
  #flop_sums <- colSums(data_merged[data_merged$group =="flop", 2:ncol(data_merged)])
  
  # Create a data frame indicating whether top or flop was higher for each column
  #comparison_df <- data.frame(
  #  otu = names(top_sums),
  #  Top = as.integer(top_sums > flop_sums),
  #  Flop = as.integer(top_sums < flop_sums)
  #)
  #res.0.1 =merge(comparison_df, res, by = "otu")
  #tax.ind =res.0.1[res.0.1$Top  ==1,]
  #MANOVA = tax.ind[order(tax.ind$pvalue, decreasing = F),] # order by p value
  #write.csv(MANOVA, file = paste0(path ,path.plant, "_MANOVA.csv"))}
  
  # aldex based on unrarefied data
  if(fct %in% c("aldex", "all")){
    conds <- ifelse(colnames(data.asv.red) %in% top$treatment, "top", "flop")
    data_int= as.data.frame(data.asv.red)
    data_int = round(data_int, 0) # Round to full numbers
    x <- aldex.clr(reads = data_int, conds = conds, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex_res <- aldex.ttest(x, paired.test = FALSE, verbose =FALSE)
    aldex_effect <-aldex.effect(x, CI =TRUE, verbose =FALSE)
    aldex_out = data.frame(aldex_res, aldex_effect)
    if (correction ==TRUE) { 
      ALDEX = aldex_out[aldex_out$wi.eBH < 0.05,]
    } else{
      ALDEX = aldex_out[aldex_out$wi.ep < 0.05,]}
    ALDEX = ALDEX[ALDEX$effect > 0, ] # Positive effect means higher in the first group which is top!
    write.csv(ALDEX, file = paste0(path , path.plant, ifelse(cutoff == TRUE, "_cut",""), "_ALDEX.csv"))
    print("ALDEX DONE")
  }
  # maaslin rarefied data
  if(fct %in% c("maaslin", "all")){
    metadata = selected_composts
    rownames(metadata)= metadata$treatment 
    data = data.ISS.red %>% t() %>% as.data.frame()
    #rownames(data) == rownames(metadata) # rownames of metadata and data must be the same!
    fit_data = Maaslin2(
      input_data = data, 
      input_metadata = metadata, 
      output = "demo_output", 
      transform = "AST",
      normalization ="TSS",
      standardize = FALSE,
      fixed_effects = "group")
    res= fit_data$results
    res = res[res$coef > 0,] # subsetting only the ones that are associated with the top composts
    if(correction == TRUE){
      MAASLIN = res[res$qval<0.05,]
    } else{
      MAASLIN = res[res$pval<0.05,]
    }
    write.csv(MAASLIN, file = paste0(path ,path.plant, ifelse(cutoff == TRUE, "_cut",""), "_MAASLIN.csv"))
    print("MAASLIN DONE")
    
  }
  #if(fct %in% c("presab", "all")){
  #  metadata = selected_composts
  #  rownames(metadata)= metadata$treatment 
  #  metadata = merge(metadata %>% select(c("group")), data.ISS.red %>% t(), by =0) # Attention merge always changes the order!
  #  row.names(metadata) = metadata$Row.names
  #  metadata$Row.names <-NULL
  #  meta_aggregate = aggregate(metadata[,-1], by = list(metadata$group), mean)
  #  meta_aggregate = meta_aggregate[meta_aggregate$Group.1 == "flop",] 
  #  meta_aggregate = meta_aggregate[,-1]
  #  presab = data.ISS.red[meta_aggregate == 0, ]}
  
  if(fct == "all"){
    # Venndiagram for the four methods
    x <- list(PBC = rownames(PBC), ALDEx2 = rownames(ALDEX), MaAsLin2 = MAASLIN$feature, ISA = rownames(ISA))
    g = ggVennDiagram(x)
    ggsave(g, file =paste0(path ,path.plant, ifelse(cutoff == TRUE, "_cut",""), "_venndiagram.png"), height= 5, width =5)
    
    # Count occurrences of each feature
    all_features <- c(rownames(PBC), rownames(ALDEX),MAASLIN$feature, rownames(ISA))
    feature_counts <- table(all_features)
    
    # Select features that occur in at least three groups
    list[["min2"]] <- names(feature_counts[feature_counts >= 2])
    list[["min3"]] <- names(feature_counts[feature_counts >= 3])
    list[["all"]] <- names(feature_counts[feature_counts >= 4])
  }
  return(list)}

# Indicator saving results helper funciton. Attention  needs tax.calss.red, path2, path.plant loaded
comparison_Indicator_analysis <-function(output, level, cutoff =TRUE){
  all = c(output$all, output$min2, output$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  
  if (level == "ASV") {
    tax_select =tax.class.red[unique(all),]
  }
  if (level == "genus") {
    tax_select = data.asv.genus[unique(all),]
  }
  if (level == "family") {
    tax_select = data.asv.family[unique(all),] 
  }
  if (level == "phyla") {
    tax_select = data.asv.phyla[unique(all),] 
  }
  
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  if (level == "ASV"){
    order =ISS.rob.comp.avg[rownames(tax_select_count), top$treatment] %>% rowSums() %>%
      sort(decreasing = TRUE) %>% as.data.frame()
  }
  
  if (level == "genus") {
    order =data.ISS.genus[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
      sort(decreasing = TRUE) %>% as.data.frame()
  }
  
  if (level == "family") {
    order =data.ISS.family[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
      sort(decreasing = TRUE) %>% as.data.frame()
  }
  
  if(level =="phyla"){
    order =data.ISS.phyla[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
      sort(decreasing = TRUE) %>% as.data.frame()
  }
  tax_select_count_order =tax_select_count[rownames(order),]
  print(output$nASV)
  write.csv(tax_select_count_order, file = paste0(path2, path.plant,ifelse(cutoff == TRUE, "_cut",""), "_method_comparison.csv"))
}

# Phylogenetical tree based on indicartive AVSvs, round
phylo_tree_ASV <- function(ASV, sequences, clade_colors, all =FALSE){
  if(all == FALSE){ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species)}
  else{ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species, assay_combination)}
  
  ASV <-ASV %>%
    mutate(lowest = case_when(
      Species != "unclassified" ~ Species,
      Genus != "unclassified" ~ Genus,
      Family != "unclassified" ~ Family,
      Order != "unclassified" ~ Order,
      Class != "unclassified" ~ Class,
      Phyla != "unclassified" ~ Phyla,
      TRUE ~ "Unclassified"
    ))
  
  ASV$lowest <- sapply(ASV$lowest, function(x) unlist(strsplit(x, split='_sp', fixed=TRUE))[1]) %>% as.vector()
  
  # read FASTA file, select only the FASTA files that are needed
  required_ids <- ASV$ASV
  subset_sequences <- sequences[names(sequences) %in% required_ids]
  
  # Convert seqinr sequences to a format compatible with msa (if aligning)
  dna_sequences <- sapply(subset_sequences, function(x) toupper(paste(x, collapse = "")))
  dna_string_set <- DNAStringSet(dna_sequences)
  
  # Create an allignment
  alignment <- msa(dna_string_set)
  aligned_sequences <- as.DNAbin(as(alignment, "DNAMultipleAlignment"))
  
  # Create a distance matrix and construct a phylogenetic tree
  dist_matrix <- dist.dna(aligned_sequences, model = "JC69")
  phylo_tree <- nj(dist_matrix)
  
  # Metadata
  if(all == FALSE){
    metadata <- data.frame(
      ASV = names(subset_sequences),  # Ensure this matches the sequence IDs
      Phylum = ASV[match(names(subset_sequences), ASV$ASV), "Phyla"],
      lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"])}
  else{
    metadata <- data.frame(
      ASV = names(subset_sequences),  # Ensure this matches the sequence IDs
      Phylum = ASV[match(names(subset_sequences), ASV$ASV), "Phyla"],
      lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"],
      BLW1 = ASV[match(names(subset_sequences), ASV$ASV), "BLW1_all"],
      combination = ASV[match(names(subset_sequences), ASV$ASV), "assay_combination"])}
  
  # Draw the tree
  
  # Create a ggtree plot
  p <- ggtree(phylo_tree, layout="circular", branch.length="none") %<+% metadata +
    geom_tippoint(aes(color = Phylum), size=4, shape =16, position =position_nudge(x = -1)) +  theme_tree2()+
    scale_color_manual(values = clade_colors$color)
  
  
  # Extract data from ggtree object
  p_data <- p$data
  
  # Calculate angles and adjust text alignment
  p_data <- p_data %>%
    mutate(angle2 = (ifelse(angle < 90 | angle > 270 , angle , angle+ 180)),
           hjust = ifelse(angle < 90 | angle > 270, -0.1, 1.1))
  
  # Add labels with rotation and alignment adjustments
  
  if(all == TRUE){
    p <- p + ggnewscale::new_scale_color() +
      geom_tippoint(aes(color = combination), size=3, shape=17) +
      scale_color_manual(values = c("cress_cugu" =   "#469990", "cugu_curs" =  "#F58231",
                                    "cress_cugu_curs" = "#9A6324")) 
  }
  legend1 = ggpubr::get_legend(p+ theme(legend.position = "right"))
  
  p1 =p +ggnewscale::new_scale_colour()  + geom_text2(aes(label = lowest, angle = angle2, hjust = hjust), data = p_data, size = 3)+
    theme(plot.margin = margin(2, 2, 2, 2, "cm"),
          axis.line = element_blank(),  # Remove axis lines
          axis.text = element_blank(),  # Remove axis text
          axis.ticks = element_blank(),  # Remove axis ticks
          panel.grid = element_blank(),  # Remove grid lines
          legend.position = "none")
  
  return(list(phylo_tree, p1, legend1))
}

# Adjusted version flat for publication
phylo_tree_ASV_flat <- function(ASV, sequences, clade_colors, all =FALSE, size =2, xlim.factor = 1.35, offset = 0.5){
  if(all == FALSE){ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species)}
  else{ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species,assay_combination)}
  
  ASV <-ASV %>%
    mutate(lowest = case_when(
      Species != "unclassified" ~ Species,
      Genus != "unclassified" ~ Genus,
      Family != "unclassified" ~ Family,
      Order != "unclassified" ~ Order,
      Class != "unclassified" ~ Class,
      Phyla != "unclassified" ~ Phyla,
      TRUE ~ "Unclassified"
    ))
  
  ASV$lowest <- sapply(ASV$lowest, function(x) unlist(strsplit(x, split='_sp', fixed=TRUE))[1]) %>% as.vector()
  
  # read FASTA file, select only the FASTA files that are needed
  required_ids <- ASV$ASV
  subset_sequences <- sequences[names(sequences) %in% required_ids]
  
  # Convert seqinr sequences to a format compatible with msa (if aligning)
  dna_sequences <- sapply(subset_sequences, function(x) toupper(paste(x, collapse = "")))
  dna_string_set <- DNAStringSet(dna_sequences)
  
  # Create an allignment
  alignment <- msa(dna_string_set)
  aligned_sequences <- as.DNAbin(as(alignment, "DNAMultipleAlignment"))
  
  # Create a distance matrix and construct a phylogenetic tree
  dist_matrix <- dist.dna(aligned_sequences, model = "JC69")
  phylo_tree <- nj(dist_matrix)
  
  # Metadata
  if(all == FALSE){
    metadata <- data.frame(
      ASV = names(subset_sequences),  # Ensure this matches the sequence IDs
      Phylum = ASV[match(names(subset_sequences), ASV$ASV), "Phyla"],
      lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"])}
  else{
    metadata <- data.frame(
      ASV = names(subset_sequences),  # Ensure this matches the sequence IDs
      Phylum = ASV[match(names(subset_sequences), ASV$ASV), "Phyla"],
      lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"],
      combination = ASV[match(names(subset_sequences), ASV$ASV), "assay_combination"])}
  
  # Draw the tree
  
  # Create a ggtree plot
  p <- ggtree(phylo_tree, branch.length="none") %<+% metadata +
    geom_tippoint(aes(color = Phylum), size=3, shape =16, position =position_nudge(x = 0)) +  theme_tree2()+
    scale_color_manual(values = clade_colors$color)
  
  
  # Extract data from ggtree object
  p_data <- p$data
  
  # Add labels with rotation and alignment adjustments
  
  if(all == TRUE){
    p <- p + ggnewscale::new_scale_color() +
      geom_tippoint(aes(color = combination), size=3, shape=17, position =position_nudge(x = -1)) +
      scale_color_manual(values = c("cress_cugu" =   "#469990", "cugu_curs" =  "#F58231",
                                    "cress_cugu_curs" = "#9A6324")) 
  }
  legend1 = ggpubr::get_legend(p+ theme(legend.position = "right",
                                        legend.text = element_text(size =9),
                                        legend.title = element_text(size =10, face ="bold")))
  
  p1 =p +ggnewscale::new_scale_colour()  + geom_tiplab(aes(label = lowest), data = p_data, size = size, offset =offset, hjust =0, align =TRUE, linetype =NA)+
    theme(axis.line = element_blank(),  # Remove axis lines
          axis.text = element_blank(),  # Remove axis text
          axis.ticks = element_blank(),  # Remove axis ticks
          panel.grid = element_blank(),  # Remove grid lines
          legend.position = "none")+
    xlim(NA, max(p$data$x) * xlim.factor)
  
  return(list(phylo_tree, p1, legend1))
}

# Compare directly among the three systems and draw the phylotrees

comparison_pathogen_plant <-function(level, robdetect, cutoff, nmethods =3) {
  path2 <- paste0(path1, level, robdetect, "/")
  
  cress <- read.csv(file =paste0(path2, "cress", ifelse(cutoff ==TRUE, "_cut","" ), "_method_comparison.csv"))
  rownames(cress) <-cress[,1]; cress[,1] <-NULL
  
  cugu <- read.csv(file =paste0(path2, "cugu", ifelse(cutoff ==TRUE, "_cut","" ),"_method_comparison.csv"))
  rownames(cugu) <-cugu[,1]; cugu[,1] <-NULL
  
  curs <- read.csv(file =paste0(path2, "curs", ifelse(cutoff ==TRUE, "_cut","" ),"_method_comparison.csv"))
  rownames(curs) <-curs[,1]; curs[,1] <-NULL
  
  
  # Initialize the list to store the results
  results_list <- list()
  
  # Define the comparison types
  comparison_types <- c(1, 2, 3)
  
  # Define assays and their combinations for comparison
  assays <- c("cress", "cugu", "curs")
  assay_combinations <- list(
    c("cress", "cugu"), 
    c("cress", "curs"), 
    c("cugu", "curs"),
    assays # Include full combination
  )
  
  # Function to create comparison list
  create_comparison_list <- function(include_assays, comparison) {
    comparison_list <- setNames(lapply(include_assays, function(assay) {
      # Fetch the correct dataset based on the assay name
      if(assay == "cress") {
        cress %>% filter(Freq >= comparison) %>% rownames()
      } else if(assay == "cugu") {
        cugu %>% filter(Freq >= comparison)  %>% rownames()
      } else if(assay <= "curs") {
        curs %>% filter(Freq >= comparison)  %>% rownames()
      }
    }), include_assays)
    
    return(comparison_list)
  }
  
  # Function to process overlap
  process_overlap <- function(x, comparison, assay_combination) {
    
    if (level == "ASV") {
      overlapp <- tax.class.red[Reduce(intersect, x),]
      overlapp$ASV <- rownames(overlapp)
    } else{
      overlapp <- Reduce(intersect, x) %>% as.data.frame()
    }
    
    if (nrow(overlapp) == 0)  {
      overlapp <-NULL
    } else{
      overlapp$comparison <- comparison
      
      # Add the assay combination as a string (e.g., "cress_cugu")
      overlapp$assay_combination <- paste(assay_combination, collapse = "_")
      
      # Generate and save the Venn diagram
      venn_plot <- ggVennDiagram(x) # Modify this based on your actual Venn diagram function
      filename <- paste0(path2, "Venndiagram_", overlapp$comparison, "_", overlapp$assay_combination, ".png")
      ggsave(filename = filename, plot = venn_plot, height = 5, width = 5)
    }
    return(overlapp)
  }
  
  # Loop over each comparison type
  for (comparison in comparison_types) {
    # Perform comparisons for all assay combinations
    for (assay_pair in assay_combinations) {
      comparison_list <- create_comparison_list(assay_pair, comparison)
      results_list <- c(results_list, list(process_overlap(comparison_list, comparison, assay_pair)))
    }
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, results_list)
  write.csv(final_results, file =paste0(path2,"comparisons", ifelse(cutoff ==TRUE, "_cut","" ), "_path_plants_systems.csv"), row.names = FALSE)
  
  if (level == "ASV") {
    
  if (nmethods == 3) {
    cress <- cress %>% filter(Freq > 1)
    cugu <- cugu %>%  filter(Freq > 1)
    curs <- curs %>%  filter(Freq > 1)
  }
  cress$ASV <- rownames(cress)
  cugu$ASV <- rownames(cugu)
  curs$ASV <- rownames(curs)
  
  names_phyla <- c(unique(cress$Phyla), unique(cugu$Phyla), unique(curs$Phyla)) %>% unique() %>% sort()
  phylacol <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "lightblue", "violet", "orange")
  colors <- data.frame(color = phylacol[1:length(names_phyla)],
                       phyla = names_phyla)
  
  clade_colors1 = colors %>% filter(phyla %in% unique(cress$Phyla)) %>% dplyr::select(color)
  p1 <-phylo_tree_ASV_flat(cress, sequences, clade_colors1, all=FALSE, size =4, xlim = 2.5, offset = 0.3)
  save_plot(
    filename = paste0("figures/", path, "indicator_species_analys/", level, robdetect, "/F_BLW2_GU_cress",
                      ifelse(cutoff == TRUE, "_cut", ""), "_phylogenetical_tree_rob", nmethods,".png"),
    plot = p1[[2]],
    base_height = 5, # Adjust these parameters for dynamic scaling
    base_width = 5
  )
  
  clade_colors2 = colors %>% filter(phyla %in% unique(cugu$Phyla)) %>% dplyr::select(color)
  p2 <-phylo_tree_ASV_flat(cugu, sequences, clade_colors2, all=FALSE, size =4, xlim = 2.5, offset = 0.3)
  save_plot(
    filename = paste0("figures/", path, "indicator_species_analys/", level, robdetect, "/F_BLW2_GU_cuc",
                      ifelse(cutoff == TRUE, "_cut", ""), "_phylogenetical_tree_rob", nmethods,".png"),
    plot = p2[[2]],
    base_height = 5, # Adjust these parameters for dynamic scaling
    base_width = 5
  )
  
  clade_colors3 = colors %>% filter(phyla %in% unique(curs$Phyla)) %>% dplyr::select(color)
  p3 <-phylo_tree_ASV_flat(curs, sequences, clade_colors3, all=FALSE, size =4, xlim = 2.5, offset = 0.3)
  save_plot(
    filename = paste0("figures/", path, "indicator_species_analys/", level, robdetect, "/F_BLW2_RS_cuc",
                      ifelse(cutoff == TRUE, "_cut", ""), "_phylogenetical_tree_rob", nmethods, ".png"),
    plot = p3[[2]],
    base_height = 5, # Adjust these parameters for dynamic scaling
    base_width = 5
  )
  combinedplot <-ggarrange(p1[[2]], p2[[2]], p3[[2]], p1[[3]], p2[[3]], p3[[3]], ncol =3, nrow =2, heights = c(1,0.2), widths = c(0.65,1,1)) #cut
  
  save_plot(
    filename = paste0("figures/", path, "indicator_species_analys/", level, robdetect, "/F_BLW2",
                      ifelse(cutoff == TRUE, "_cut", ""),"_combined_plot_rob", nmethods, ".png"),
    plot = combinedplot,
    base_height = 12,  # Adjust the height according to your combined plot
    base_width = 12     # Adjust the width as needed
  )
  }
}

# Heat map for indicative ASVs
heatmap_all <- function(topflop, asv_names, ds, tip_order, presentation = FALSE, size = 8, adjustcor = FALSE){
  data.asv <-asv.rob.comp.avg %>% as.matrix() # Based on unrarefied asv table
  data.asv.TSS <- prop.table(data.asv, margin = 1) %>% as.data.frame() # Calculate relative abundance
  
  ind.ASV.sum <- matrix(NA, nrow = length(asv_names), ncol = 25) %>% as.data.frame() # create empty data.frame
  colnames(ind.ASV.sum)<- c("ASV", "tax", "n",
                            "T9_cress", "F9_cress", "ratio.9.9_cress",
                            "T9_cugu", "F9_cugu", "ratio.9.9_cugu",
                            "T9_curs", "F9_curs", "ratio.9.9_curs",
                            "T3", "F4", "ratio.3.4",
                            "rho_cress", "p_cress", "padj_cress", "rho_cugu", "p_cugu", "padj_cugu", "rho_curs", "p_curs",
                            "padj_curs", "sites")
  ind.ASV.sum$ASV <- asv_names
  min.value <- min(data.asv.TSS[data.asv.TSS >0])/2 # Half of the minimum abundance (to avoid zeros in the data set)
  data.asv.TSS[data.asv.TSS == 0] <- min.value
  
  ind.ASV.sum$tax <-tax.class.red[asv_names,] %>%
    mutate(tax = case_when(
      Species != "unclassified" ~ Species,
      Genus != "unclassified" ~ Genus,
      Family != "unclassified" ~ Family,
      Order != "unclassified" ~ Order,
      Class != "unclassified" ~ Class,
      Phyla != "unclassified" ~ Phyla,
      TRUE ~ "Unclassified"
    )) %>% select("tax")
  
  # Number of composts present
  for (i in 1:nrow(ind.ASV.sum)) {
    ind.ASV.sum[i, "n"] <-(data.asv.TSS[, asv_names[i]] >min.value) %>% sum()
  }
  
  # Abundance in composts
  factors <- c("T9_cress", "F9_cress","T9_cugu", "F9_cugu","T9_curs", "F9_curs","T3", "F4")
  
  for(j in 1:length(factors)){
    for (i in 1:nrow(ind.ASV.sum)) {
      ind.ASV.sum[i,factors[j]] <-data.asv.TSS[topflop[[factors[j]]]$treatment, asv_names[i]]  %>% mean() %>% log10() %>% round(2)
    }
  }
  
  # Calculates ratios
  ind.ASV.sum$ratio.9.9_cress <-((10^ind.ASV.sum$T9_cress)/(10^ind.ASV.sum$F9_cress)) %>% round(2)
  ind.ASV.sum$ratio.9.9_cugu <-((10^ind.ASV.sum$T9_cugu)/(10^ind.ASV.sum$F9_cugu)) %>% round(2)
  ind.ASV.sum$ratio.9.9_curs <-((10^ind.ASV.sum$T9_curs)/(10^ind.ASV.sum$F9_curs)) %>% round(2)
  ind.ASV.sum$ratio.3.4 <-((10^ind.ASV.sum$T3)/(10^ind.ASV.sum$F4)) %>% round(2)
  
  # present in which sites?
  for (i in 1:nrow(ind.ASV.sum)) {
    composts <-data.asv.TSS[data.asv.TSS[, asv_names[i]] >min.value,] %>% rownames()
    sites <-df_BLW2.div %>% filter(treatment %in% composts) %>% select(site_ID) %>% unique()
    ind.ASV.sum[i, "sites"] <- sites$site_ID %>% sort() %>% paste0(collapse = " ")
  }
  
  # Correlation test
  rho.name <- c("rho_cress", "rho_cugu", "rho_curs")
  p.name <- c("p_cress", "p_cugu", "p_curs")
  p.adj <- c("padj_cress", "padj_cugu", "padj_curs")
  rho.sig <- c("rho.sig.cress", "rho.sig.cugu", "rho.sig.curs")
  
  for (j in 1:3) {
    for (i in 1:nrow(ind.ASV.sum)) {
      cor <-cor.test(data.asv.TSS[ , asv_names[i]], df_BLW2.div[,ds[j]], method ="spearman")
      ind.ASV.sum[i, rho.name[j]] <- cor$estimate %>% round(2)
      ind.ASV.sum[i, p.name[j]] <- cor$p.value %>% round(3)
    }
    
    # Adjust for multiple testing
    if (adjustcor == TRUE) {
      ind.ASV.sum[, p.adj[j]] <- p.adjust(ind.ASV.sum[, p.name[j]], method = "BH")
      ind.ASV.sum[, rho.sig[j]] <-ifelse(ind.ASV.sum[, p.adj[j]] < 0.05 , ind.ASV.sum[, rho.name[j]], NA)
      
    } else {
      ind.ASV.sum[, rho.sig[j]] <-ifelse(ind.ASV.sum[, p.name[j]] < 0.05 , ind.ASV.sum[, rho.name[j]], NA)
    }}
  
  
  rownames(ind.ASV.sum) <- ind.ASV.sum$ASV; ind.ASV.sum$ASV <- NULL
  
  # Sort by the order in the phylogenetical tree
  ind.ASV.sum <- ind.ASV.sum[rev(tip_order),]
  
  # Create heat maps
  column_width <- unit(0.5, "cm")
  heatmap_data <- ind.ASV.sum[, factors]
  heatmap_data[heatmap_data==log10(min.value) %>% round(2)] <-NA
  
  color_ramp <- colorRampPalette(c("steelblue", "darkblue"))(100)
  
  H1 <-Heatmap(
    as.matrix(ind.ASV.sum %>% select(starts_with("rho.sig"))),
    na_col = "lightgrey",
    name = "Sperman's\nRho",
    col = color_ramp,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = gpar(fontsize = size, fontface = "bold"),
    column_names_rot = 0,
    row_names_gp = gpar(fontsize = size),
    row_names_side = "left",
    width = column_width*3, 
    heatmap_legend_param = list(title_gp = gpar(fontsize = size, fontface = "bold", just = "center")),
   column_labels = c("I", "II", "III")
  )
  
  if (presentation == TRUE) {
    color_ramp <- colorRampPalette(c("lightcoral", "darkred"))(100)
    H2 <- Heatmap(
      as.matrix(heatmap_data),
      name = "Log10",
      na_col = "grey",
      col = color_ramp,
      show_row_names = TRUE,
      show_column_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_gp = gpar(fontsize = size, fontface = "bold"),
      row_names_gp = gpar(fontsize = size, fontface = "bold"),
      width = column_width*ncol(heatmap_data), 
      heatmap_legend_param = list(title_gp = gpar(fontsize = size, fontface = "bold")) )
  }
  
  
  color_ramp <-colorRampPalette(c("lightgreen", "darkgreen"))(100)
  
  H3 <- Heatmap(
    as.matrix(ind.ASV.sum %>% select(n)),
    name = "Number of\ncomposts",
    col = color_ramp,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = gpar(fontsize = size, fontface = "bold"),
    row_names_gp = gpar(fontsize = size, fontface = "bold"),
    width = column_width, 
    heatmap_legend_param = list(title_gp = gpar(fontsize = size, fontface = "bold"))
  )
  
  # Heat maps for sites presence
  ind.ASV.sum$A <-grepl("A", ind.ASV.sum$sites, fixed = TRUE)
  ind.ASV.sum$B <-grepl("B", ind.ASV.sum$sites, fixed = TRUE)
  ind.ASV.sum$C <-grepl("C", ind.ASV.sum$sites, fixed = TRUE)
  ind.ASV.sum$D <-grepl("D", ind.ASV.sum$sites, fixed = TRUE)
  ind.ASV.sum$E <-grepl("E", ind.ASV.sum$sites, fixed = TRUE)
  ind.ASV.sum$F <-grepl("F", ind.ASV.sum$sites, fixed = TRUE)
  ind.ASV.sum$G <-grepl("G", ind.ASV.sum$sites, fixed = TRUE)
  
  # Create a logical matrix from the presence/absence data
  presence_absence_df <- ind.ASV.sum %>% select(A,B,C,D,E,F,G)
  # Define color mapping for presence/absence
  
  # Create the left-side annotation for presence/absence
  if (presentation == FALSE) {
    left_annotation <- rowAnnotation(
      df = presence_absence_df,
      col = list(
        A = c("TRUE" = site.colors[1], "FALSE" = "white"),
        B = c("TRUE" = site.colors[2], "FALSE" = "white"),
        C = c("TRUE" = site.colors[3], "FALSE" = "white"),
        D = c("TRUE" = site.colors[4], "FALSE" = "white"),
        E = c("TRUE" = site.colors[5], "FALSE" = "white"),
        F = c("TRUE" = site.colors[6], "FALSE" = "white"),
        G = c("TRUE" = site.colors[7], "FALSE" = "white")),
      annotation_name_gp = gpar(fontsize = size, fontface = "bold"),
      annotation_name_rot = 0, 
      show_annotation_name = TRUE,
      show_legend = c("bar" = FALSE))
    
    heat <- H3 + left_annotation + H1 +H2
    heatmap_ggplot <- draw(heat, heatmap_legend_side = "right")
    heatmap_grob <- grid::grid.grabExpr(draw(heatmap_ggplot))
    return(heatmap_grob)  
    
  }else{
    left_annotation <- rowAnnotation(
      df = presence_absence_df,
      col = list(
        A = c("TRUE" = "lightgrey", "FALSE" = "darkred"),
        B = c("TRUE" = "lightgrey", "FALSE" = "darkred"),
        C = c("TRUE" = "lightgrey", "FALSE" = "darkred"),
        D = c("TRUE" = "lightgrey", "FALSE" = "darkred"),
        E = c("TRUE" = "lightgrey", "FALSE" = "darkred"),
        F = c("TRUE" = "lightgrey", "FALSE" = "darkred"),
        G = c("TRUE" = "lightgrey", "FALSE" = "darkred")),
      annotation_name_gp = gpar(fontsize = size, fontface = "bold"),
      annotation_name_rot = 0, 
      show_annotation_name = TRUE,
      show_legend = FALSE
    )
    
    combined_legend <- Legend(
      labels = c("Present", "Absent"),
      title = "Composting\nSite",
      legend_gp = gpar(fill = c("lightgrey", "darkred")),
      title_gp = gpar(fontsize = size, fontface = "bold"),
      labels_gp = gpar(fontsize = size)
    )
    
    heat <- H3 + left_annotation + H1 
    heatmap_ggplot <- draw(
      heat, 
      heatmap_legend_side = "right",
      heatmap_legend_list = list(combined_legend)  # Add the custom legend to the main legend column
    )
    heatmap_grob <- grid::grid.grabExpr(draw(heatmap_ggplot))
    return(list(heatmap_grob, ind.ASV.sum))  
    
  }
  
}

# Remove all unnecessary factors etc.
rm(size, background_pr, mypalette)


# Data------
## General------

set.seed(100)  # set random seed

# compost to exclude and factors
compost.exclude = c("K6","K10","K17","K27","K37")
factor.category <- c("site", "company", "comp.system")  # categorical factors
factor.continous <-c( "DS",  "max_WHC", "pH", "sal", "OD550", "NO2", "NO3","NH4",
                      "Nmin","NO3.Nmin","PO4","Ntot", "Ctot","Corg","Corg.N", "basal","FDA")
factor.bioassay <- c("rel.cp","gp.cp","biomass.cp", "bca.cp","ratio.cp","rank.cp","rank.overall.cp",
                     "biomass.cu.gu","rel.cu.gu", "surv.cu.gu", "bca.cu.gu", "ratio.cu.gu", "avg.cu.gu",
                     "gp.cu.gu", "rank.cu.gu", "rank.overall.cu.gu",
                     "biomass.cu.rs", "rel.cu.rs", "surv.cu.rs", "bca.cu.rs", "ratio.cu.rs", "avg.cu.rs",
                     "gp.cu.rs", "rank.cu.rs", "rank.overall.cu.rs")
factor.site <- c("age",  "plant.content", "soil.content", "max.temp")

# Meta data
df_BLW2 = read.csv(file ="data/meta_data_compost_BLW_without_errors.csv", sep =";")
factors = c("batch", "site", "company", "comp.system", "site_ID") # Change to factor
for (i in factors) {
  df_BLW2[, i] <- as.factor(df_BLW2[,i])
}
df_BLW2.r = df_BLW2 %>% filter(!treatment %in% compost.exclude) # Reduce normal distribution to 37 composts
rm(factors)


# Compost_ID for changing the labeling of the data frame
compost_ID <- read.csv(file ="data/Key_treatment_compostID.csv", sep= ";")
compost_ID$treatment = as.factor(compost_ID$treatment)
compost_ID$compost_ID = as.factor(compost_ID$compost_ID)
compost_ID$site_ID = as.factor(compost_ID$site_ID)
compost_ID <- compost_ID %>% filter(!treatment %in% c("K27", "K37")) # Exclude two error composts
compost_ID$treatment <-droplevels(compost_ID$treatment)

# Key for sequencing samples
design.BLW2 <- read.csv(file ="data/design_BLW2.csv")
design.BLW2$site <- as.factor(design.BLW2$site)

# Add the new numbers to the design.BLW2
design.BLW2 = merge(design.BLW2, compost_ID, by ="treatment", ) # Add new labeling
rownames(design.BLW2) <- design.BLW2$ID
design.BLW2$ID <- NULL
design.BLW2 = design.BLW2[order(rownames(design.BLW2)),]
design.BLW2$batch.y <-NULL
setnames(design.BLW2, "batch.x", "batch")
design.BLW2 = design.BLW2 %>% filter(!treatment %in% compost.exclude) # 37 compost

# For plots

# For the bioassay plots
graphic.style = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, size =14),
                      legend.position = "right",
                      axis.text.y = element_text(size =16),
                      axis.title.y = element_text(size=16),
                      axis.ticks.y = element_line(colour = "lightgrey"),
                      axis.line = element_line(colour = "black"),
                      panel.grid.major.y = element_line(colour ="lightgrey"),
                      panel.grid.major.x= element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_rect(fill ="white", colour ="white"))

bg_theme <- theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.3, vjust = 0, size = 14, color ="black"),
    axis.text.y = element_text(size = 14, color ="black"),
    text = element_text(size = 14, color ="black"),
    axis.title = element_text(size=14),
    legend.text = element_text(size=14)
  )

site.colors =  c("#81d4fa", "#01579b", "#ffecb3", "#e65100","#b39ddb","#ff9800", "#0288d1", "lightgrey")
site.colors.dark =  c("#81d4fa","#01579b","#9D8604FF","#e65100","#b39ddb", "#ff9800", "#0288d1")  

# Color for bioassays
color.cress.dark <- rgb(143, 170, 220,maxColorValue = 255)
color.cress.light <- rgb(189, 215, 238, maxColorValue =255)
color.cugu.dark <- rgb(169, 209, 142, maxColorValue =255)
color.cugu.light <- rgb(197, 224, 180, maxColorValue =255)
color.curs.dark <- rgb(244, 177, 131, maxColorValue =255)
color.curs.light <- rgb(248, 203, 173, maxColorValue =255)
color.dark <- c(color.cress.dark, color.cugu.dark, color.curs.dark)
color.light <-c(color.cress.light, color.cugu.light, color.curs.light)

rm(color.cress.dark, color.cress.light, color.cugu.dark, color.cugu.light, color.curs.dark,
   color.curs.light)


## Sequencing----
# asv file without rarefying
asv = read.table("data/asv.F.BLW2.txt")
ISS = read.table(paste0("data/", path,"ISS.F.BLW2.txt"))
tax = read.table("data/tax.F.BLW2.txt")
tax.r = read.table("data/tax.F.r.txt") # Without compost exclude, but with peat!
asv.rob.comp.avg = read.table(paste0("data/",path,"asv.rob.comp.avg.txt"))
asv.rob.comp.avg2 = read.table(paste0("data/",path,"asv.rob.comp.avg_2.txt"))

# Preparation tax file
tax.class.red <- lapply(tax.r, function(x) gsub("_fam_Incertae_sedis", "_IS", x))
tax.class.red <- lapply(tax.class.red, function(x) sub(".*__", "", x)) %>% as.data.frame()
rownames(tax.class.red) <- rownames(tax.r)

# rarefied
ISS.rob = read.table(paste0("data/",path,"ISS.rob.txt"))
ISS.rob.avg = read.table(paste0("data/",path, "ISS.rob.avg.txt"))
ISS.rob.comp = read.table(paste0("data/",path,"ISS.rob.comp.txt"))
ISS.rob.comp.avg = read.table(paste0("data/",path,"ISS.rob.comp.avg.txt"))

# rarefied 2/4 replicates
ISS.rob2 = read.table(paste0("data/",path,"ISS.rob_2.txt"))
ISS.rob.avg2 = read.table(paste0("data/",path, "ISS.rob.avg_2.txt"))
ISS.rob.comp2 = read.table(paste0("data/",path,"ISS.rob.comp_2.txt"))
ISS.rob.comp.avg2 = read.table(paste0("data/",path,"ISS.rob.comp.avg_2.txt"))

ISS.comp.avg <- read.table(paste0("data/",path, "ISS.F.BLW2.comp.txt")) #Rarefied but not only the robust ones
design.BLW2.r <- design.BLW2[rownames(design.BLW2) %in% colnames(ISS.comp.avg),]
ISS.comp.avg = ISS.comp.avg %>% t() %>% aggregate(list(design.BLW2.r$treatment), mean) # calculate mean for the four replicates
rownames(ISS.comp.avg) <- ISS.comp.avg$Group.1; ISS.comp.avg$Group.1 <-NULL
ISS.comp.avg <- ISS.comp.avg %>% t() %>% as.data.frame()
(rowSums(ISS.comp.avg) ==0) %>% sum() # Check if all the ASVs without appreance are gone

ISS.comp.avg = as.matrix(ISS.comp.avg) %>% t()
ISS.prop = prop.table(ISS.comp.avg, margin = 1) * 100  # get proportions
ISS.prop = as.data.frame(ISS.prop)
ISS.comp.avg = as.data.frame(ISS.comp.avg) %>% t()
ISS.comp.avg.t <- ISS.comp.avg %>% t()

# Alpha diversity for 37 compost + peat substrate

# Meta file with averages
df_BLW2.div = read.csv(paste0("data/",path,"df_BLW2.div.csv"))
df_BLW2.div$X <-NULL
df_BLW2.div$batch <- as.factor(df_BLW2.div$batch)

# Replicates separately with peat
ISS.alpha = read.csv(paste0("data/",path,"ISS.alpha.F.BLW2.csv"))
rownames(ISS.alpha) <- ISS.alpha$X; ISS.alpha$X <- NULL

# Beta diversity 37 composts (using 100 rarefied abundance table)
ISS.bray <- read.table(paste0("data/",path,"ISS.bray.txt"))
ISS.bray.avg <- read.table(paste0("data/",path,"ISS.bray.avg.txt")) %>% as.matrix() #BC-dis mean of four replicates 
ISS.bray.avg.comp <-ISS.bray.avg[!grepl("Std", rownames(ISS.bray.avg)),!grepl("Std", colnames(ISS.bray.avg))]

ISS.jac <- read.table(paste0("data/",path, "ISS.jac.txt"))
ISS.jac.avg <- read.table(paste0("data/",path, "ISS.jac.avg.txt"))
ISS.jac.avg.comp <- ISS.jac.avg[!grepl("Std", rownames(ISS.jac.avg)),!grepl("Std", colnames(ISS.jac.avg))]

# Fasta file ASVs
sequences <- read.fasta(file = "Sequencing_data/8.all.ASV_ITSx.ITS2.fasta", seqtype = "DNA")

# Analysis------
## Sequencing overview----
# rarecurve(t(asv), step = 100, col = "blue", cex = 0.6) # ASVs

tax.comp.avg <- tax.r[ISS.comp.avg.t %>% colnames(),]

#sink(file =paste0("output/",path, "Taxonomy_summary.txt"))
paste("Phyla")
tax.comp.avg$Phyla %>% table() %>% sort() %>% nrow()
tax.comp.avg$Phyla %>% table() %>% sort() %>% as.data.frame()
1/nrow(tax.comp.avg)* (nrow(tax.comp.avg)-(16+273))*100

paste("Families")
tax.comp.avg$Family %>% table() %>% sort() %>% nrow()
tax.comp.avg$Family %>% table() %>% sort() %>% as.data.frame()
1/nrow(tax.comp.avg)* (nrow(tax.comp.avg)-(191+559))*100

paste("Genera")
tax.comp.avg$Genus %>% table() %>% sort() %>% nrow()
tax.comp.avg$Genus %>% table() %>% sort() %>% as.data.frame()
1/nrow(tax.comp.avg)* (nrow(tax.comp.avg)-(243+732))*100
#sink()

# Present in all 37 composts

ISS.comp.avg.all <- ISS.comp.avg.t[, apply(ISS.comp.avg.t, 2, function(col) !any(col == 0))]
#sink(file =paste0("output/",path, "Paritioning_of_ASVs_among_composts.txt"))
paste("In total", ncol(ISS.comp.avg.all), "ASVs were present in all 36 composts. That is",
      round(ncol(ISS.comp.avg.all)/ ncol(ISS.comp.avg.t) *100, 2), "percent of all", 
      ncol(ISS.comp.avg.t), "ASVs")
paste("Classifcation on family level")
tax.r[colnames(ISS.comp.avg.all),"Family"] %>% table() %>% sort() %>% as.data.frame()
paste("Classifcation on genera level")
tax.r[colnames(ISS.comp.avg.all),"Genus"] %>% table() %>% sort() %>% as.data.frame()

paste("Present in 9 out of 10 composts")
threshold <- 33
ISS.comp.avg.all <- ISS.comp.avg.t[, apply(ISS.comp.avg.t, 2, function(col) sum(col != 0) >= threshold)]
ncol(ISS.comp.avg.all)
ncol(ISS.comp.avg.all)/ ncol(ISS.comp.avg.t) *100
tax.r[colnames(ISS.comp.avg.all),"Family"] %>% table() %>% sort() %>% as.data.frame()
tax.r[colnames(ISS.comp.avg.all),"Genus"] %>% table() %>% sort() %>% as.data.frame()

paste("Present in at least 10 composts")
threshold <- 10
ISS.comp.avg.all <- ISS.comp.avg.t[, apply(ISS.comp.avg.t, 2, function(col) sum(col != 0) >= threshold)]
ncol(ISS.comp.avg.all)
ncol(ISS.comp.avg.all)/ ncol(ISS.comp.avg.t) *100
tax.r[colnames(ISS.comp.avg.all),"Family"] %>% table() %>% sort() %>% as.data.frame()
tax.r[colnames(ISS.comp.avg.all),"Genus"] %>% table() %>% sort() %>% as.data.frame()

#sink()

rm(ISS.comp.avg.all, threshold)

# Venndiagram preperation

ISS.comp.avg.t <- as.data.frame(ISS.comp.avg.t)
ISS.comp.avg.t$treatment <- rownames(ISS.comp.avg.t)
merged_data =merge(df_BLW2[,c("treatment", "batch", "site")],ISS.comp.avg.t , by = "treatment")
rownames(merged_data) <- merged_data$treatment
merged_data$treatment <- NULL
ISS.comp.avg.t$treatment <-NULL
## Composting site----

composting_sites <- c("Riehen", "Leibstadt", "Uster", "Fehraltdorf", "Bergdietikon")
asv_tables_by_site <- lapply(composting_sites, function(site) {
  subset_data <- merged_data[merged_data$site == site, ]
  asv_table_subset <- subset_data[, -c(1, 2)]  # Remove metadata columns
  return(asv_table_subset)
})
names(asv_tables_by_site) <- composting_sites

# Apply the function to each site's ASV table
asv_tables_by_site_filtered <- lapply(asv_tables_by_site, discard_missing_asvs)

Reduce(intersect,lapply(asv_tables_by_site_filtered, colnames)) %>% length() # Number of ASVs shared in all composting sites

# Extract ASV names from each site
# Increase margin or position around the plot

asv_names_by_site <- lapply(asv_tables_by_site_filtered, colnames)
custom_colors <- c("#81d4fa", "#01579b", "#ffecb3", "#e65100", "#b39ddb")
# Create a list of sets for the Venn diagram
sets <- lapply(asv_names_by_site, function(x) as.character(x))
venn1 <- venn.diagram(
  x = sets,
  category.names = names(asv_tables_by_site_filtered),
  filename = NULL,
  output = TRUE,
  fill = custom_colors, 
  margin = 0.1,
  cat.dist = c(0.2, 0.25, 0.18, 0.2, 0.25)
)


# Can be used to check whether certain microbes belong to a specific group
partion =get.venn.partitions(sets, force.unique = TRUE, keep.elements = TRUE,
                             hierarchical = FALSE)
partion$..count.. # Number of ASVS in the different groups
ASV_names_different_groups =partion$..values.. # For the different groups the ASV names
ASV_names_groups =partion$..set.. # Groups
#ASV_names_different_groups[[2]] # How to access the different groups


## Batch-----

batch <- c(1:4)
asv_tables_by_batch <- lapply(batch, function(batch) {
  subset_data <- merged_data[merged_data$batch == batch, ]
  asv_table_subset <- subset_data[, -c(1, 2)]  # Remove metadata columns
  return(asv_table_subset)
})
names(asv_tables_by_batch) <- batch

# Apply the function to each batch's ASV table
asv_tables_by_batch_filtered <- lapply(asv_tables_by_batch, discard_missing_asvs)

Reduce(intersect,lapply(asv_tables_by_batch_filtered, colnames)) %>% length() # Number of ASVs shared in all composting batchs

# Extract ASV names from each batch
# Increase margin or position around the plot

asv_names_by_batch <- lapply(asv_tables_by_batch_filtered, colnames)
custom_colors <- c("gray80","gray60", "gray40", "gray20" )
# Create a list of sets for the Venn diagram

sets <- lapply(asv_names_by_batch, function(x) as.character(x))
venn2 <- venn.diagram(
  x = sets,
  category.names = names(asv_tables_by_batch_filtered),
  filename = NULL,
  output = TRUE,
  fill = custom_colors,
  margin = 0.1
)

# Can be used to check whether certain microbes belong to a specific group
partion =get.venn.partitions(sets, force.unique = TRUE, keep.elements = TRUE,
                             hierarchical = FALSE)
partion$..count.. # Number of ASVS in the different groups
ASV_names_different_groups =partion$..values.. # For the different groups the ASV names
ASV_names_groups =partion$..set.. # Groups
#ASV_names_different_groups[[2]] # How to access the different groups


#png(paste0("figures/",path, "Sup4_Venn_site_batch.png"), width = 4000, height = 2000 , res=300)

# Arrange and save the plots
grid.arrange(
  grobs = list(venn1, venn2),
  ncol = 2  # Number of columns
)

# Close the PNG device
#dev.off()

rm(custom_colors, asv_tables_by_site, asv_tables_by_site_filtered, sets, composting_sites, partion,
   ASV_names_different_groups,ASV_names_groups, venn1, venn2, asv_names_by_batch, batch,
   merged_data, asv_names_by_site)

# Figure 4 Alpha diversity-----

# Good's coverage
data <- asv
data <- ISS

Goods <- rep(NA, ncol(data))
for (i in 1:ncol(data)) {
  Goods[i] = 1 - (sum(data[,i] ==1) / sum(data[,i]))
}
summary(Goods)
cbind(Goods, colnames(data)) %>% View()
sd(Goods)*100

data <- design.BLW2 %>% filter(!treatment %in% compost.exclude)
data <- data[ISS %>% colnames(),]
data1 <- ISS.alpha[rownames(ISS.alpha) %in% rownames(data),]
rownames(data) == rownames(data1)

#sink(file = paste0("output/", path, "Summary_ASV_richness.txt"))
paste("Peat substrate")
data2[data1 %>% rownames() %>% startsWith("Std"),] %>% select(sobs, evenness, shannon, invsimpson) %>% summary()
paste("all 37 composts without K44_2, K41_2 and K19_1-4")
data2[data1 %>% rownames() %>% startsWith("K"),] %>% select(sobs, evenness, shannon, invsimpson) %>% summary()
paste("Standard deviation for sobs, evenness, shannon and inversed simpson")
select <- data1[data1 %>% rownames() %>% startsWith("K"),] %>% select(sobs, evenness, shannon, invsimpson) 
sd(select$sobs); sd(select$evenness); sd(select$shannon); sd(select$invsimpson)
#sink()


# Figure 3
g1 =ggplot(data1 , aes(x = reorder(data$compost_ID, evenness, na.rm =T), y = evenness, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
  labs(x = "", y ="Fungal evenness")+
  scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
  scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0,0.8, by = 0.1))+
  bg_theme +
  theme(legend.position = "bottom", legend.key = element_blank(), legend.title = element_text(size =14),
        legend.text = element_text(size=14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = "lightgrey"),
        axis.title = element_text(size=14))

#ggsave(file =paste0("figures/", path,"F_Evenness_colored_by_site.png"), height =6, width = 14)


# Supplementary Alpha diversity

s1 <-ggplot(data1 , aes(x = reorder(data$compost_ID, sobs, na.rm =T), y = sobs, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
  labs(x = "", y ="Fungal richness")+
  scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
  scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 50))+
  bg_theme +
  theme(legend.position = "none", legend.key = element_blank(), legend.title = element_text(size =14),
        legend.text = element_text(size=14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = "lightgrey"))

s2 <- ggplot(data1 , aes(x = reorder(data$compost_ID, shannon, na.rm =T), y = shannon, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
  labs(x = "", y ="Fungal Shannon diversity")+
  scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
  scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5))+
  bg_theme +
  theme(legend.position = "none", legend.key = element_blank(), legend.title = element_text(size =14),
        legend.text = element_text(size=14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = "lightgrey"))


legend <- get_legend(s1+ theme(legend.position = "bottom"))
plot_grid(s1, s2, legend, nrow =3, axix ="none", rel_heights = c(1,1,0.15),
          ncol =1,labels =c("A", "B")) 

#ggsave(filename =paste0("Figures/",path, "/F_Sup4_Alpha_diversity_adapted.png"), height =10, width =14)

#Statistics for Supplementary file

data1$treatment <-gsub("_.*", "",rownames(data1) )

ANOVA <-aov(evenness ~ treatment, data= data1)
summary(ANOVA)

#sink(paste0("output/", path, "Compsoting_site_batch_effec_alpha_div.txt"))

paste("Differences in alpha diversity metrics")

paste("Composting site")
data.test <- df_BLW2.div %>% filter(!site_ID %in% c("F", "G"))

paste("sobs")
kruskal.test(mean.sobs ~ site_ID, data = data.test)
pairwise.wilcox.test(data.test$mean.even, data.test$site_ID, p.adj ="BH") 

paste("Evenness")
kruskal.test(mean.even ~ site_ID, data = data.test)
pairwise.wilcox.test(data.test$mean.even, data.test$site_ID, p.adj ="BH") 

paste("Shannon diversity")
kruskal.test(mean.shannon ~ site_ID, data = data.test)
pairwise.wilcox.test(data.test$mean.shannon, data.test$site_ID, p.adj ="BH")

paste("Inversed Simpsons")
kruskal.test(mean.ivsimp ~ site_ID, data = data.test)
pairwise.wilcox.test(data.test$mean.ivsimp, data.test$site_ID, p.adj ="BH")

paste("batch")

paste("sobs")
kruskal.test(mean.sobs ~ batch, data = df_BLW2.div)
pairwise.wilcox.test(df_BLW2.div$mean.sobs, df_BLW2.div$batch, p.adj ="BH") # Pairwise differences not significant

paste("evenness")
kruskal.test(mean.even ~ batch, data = df_BLW2.div)
pairwise.wilcox.test(df_BLW2.div$mean.even, df_BLW2.div$batch, p.adj ="BH") 

paste("Shannon diversity")
kruskal.test(mean.shannon ~ batch, data = df_BLW2.div)
pairwise.wilcox.test(df_BLW2.div$mean.shannon, df_BLW2.div$batch, p.adj ="BH") # Pairwise differences not significant

paste("Iversed simpson")
kruskal.test(mean.ivsimp ~ batch, data = df_BLW2.div)
pairwise.wilcox.test(df_BLW2.div$mean.ivsimp, df_BLW2.div$batch, p.adj ="BH") # Pairwise differences not significant
#sink()

# Homogeneity not given for all of them -> spearman's rank sum correlation
shapiro.test(df_BLW2.div$mean.sobs)
shapiro.test(df_BLW2.div$mean.shannon)
shapiro.test(df_BLW2.div$mean.even)
shapiro.test(df_BLW2.div$mean.ivsimp)
shapiro.test(df_BLW2.div$rel.cp)
shapiro.test(df_BLW2.div$rel.cu.gu)
shapiro.test(df_BLW2.div$rel.cu.rs)

list <- list()
path.plant <- c("rel.cp", "rel.cu.gu", "rel.cu.rs")
color.dark.dark <- c("#64769A", "#769263", "#AA7B5B")
alpha.div <- c("mean.sobs", "mean.even", "mean.shannon", "mean.ivsimp")
yaxis <- c("ASV richness", "Evenness", "Shannon diversity", "Inversed Simpson")
for(j in 1:4){
for (i in 1:3) {
  list[[i]] <- plot_scatter(df_BLW2.div, variable1 = path.plant[i], variable2 =  alpha.div[j], xaxis = "Disease suppression [%]", yaxis =yaxis[j] )+
    geom_text_repel(aes(label = Number, color = as.factor(path.plant[i]) ), show.legend =F,  max.overlaps =20)+
    scale_color_manual(values = color.dark[i])
}
ggarrange(list[[1]], list[[2]], list[[3]], ncol =3, labels = c("A", "B", "C"))
#ggsave(filename = paste0("figures/", path, "F_BLW2_cor_ds_",alpha.div[j], ".png"), height =4, width =12, units= "in")
}

g5 <- plot_scatter(df_BLW2.div, variable1 = "rel.cp", variable2 = "mean.even", xaxis = "Disease suppression GU-cress [%] ", yaxis = "Fungal evenness" )+
  geom_text_repel(aes(label = Number, show.legend =F,  max.overlaps =20))+
  geom_smooth(se =F, method ="lm", colour = "#64769A" )+ bg_theme


g2 = ggplot(df_BLW2.div) +
  geom_point(aes(y = mean.sobs, x = rel.cp, color = "rel.cp", shape ="rel.cp"), size =2.5) +
  geom_smooth(aes(y = mean.sobs, x = rel.cp), method = "lm", se = F, color= color.dark.dark[1], linewidth =1) +
  geom_point(aes(y = mean.sobs, x = rel.cu.gu, color = "rel.cu.gu", shape = "rel.cu.gu"), size=2.5) +
  geom_smooth(aes(y = mean.sobs, x = rel.cu.gu), method = "lm", se = F, color = color.dark.dark[2], linewidth =1) +
  geom_point(aes(y = mean.sobs, x = rel.cu.rs, color = "rel.cu.rs", shape ="rel.cu.rs"), size =2.5) +
  geom_smooth(aes(y = mean.sobs, x = rel.cu.rs), method = "lm", se= F, color = color.dark.dark[3], linewidth =1) +
  labs(x = "Disease suppression [%]", y = "") +
  scale_color_manual(name = "", values = c("rel.cp" = color.dark[1], "rel.cu.gu" = color.dark[2], "rel.cu.rs" = color.dark[3]),
                     labels = c("rel.cp" = "GU-cress", "rel.cu.gu" = "GU-cucumber", "rel.cu.rs" = "RS-cucumber")) +
  scale_shape_manual(name = "", values = c("rel.cp" = 19, "rel.cu.gu" = 15, "rel.cu.rs" = 17),
                     labels = c("rel.cp" = "GU-cress", "rel.cu.gu" = "GU-cucumber", "rel.cu.rs" = "RS-cucumber")) +
  guides(color = guide_legend(title = "")) + bg_theme+
  theme(legend.position = "top", legend.key.size = unit(0, "cm"), 
        #legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.box.margin = margin(t = 0, r = -0, b = -0.5, l = 0, unit = "cm")
  )

#ggsave(g2, filename = paste0("figures/",path, "F_BLW2_sobs_disease_suppression.png"), height =5, width = 6)


# Composting site

# Statistics
kruskal_result <- kruskal.test(mean.even ~ site_ID, data = df_BLW2.div %>% filter(!site_ID %in% c("F", "G")))
pairwise_results <-pairwise.wilcox.test(data.test$mean.even, data.test$site_ID, p.adj ="BH") 
KRUSKAL =tri.to.squ(pairwise_results$p.value)
compact_letters_df =multcompLetters(KRUSKAL)

value_max = df_BLW2.div %>% filter(!site_ID %in% c("F", "G")) %>% group_by(site_ID) %>% dplyr::summarize(max_value = max(mean.even),
                                                                 quantile = quantile(mean.even, probs = 0.75))
value_max$site_ID <- as.factor(value_max$site_ID)
value_max$letter <- c("a", "ab", "ab", "ab", "b")

# Plot
g3 = ggplot(df_BLW2.div %>% filter(!site_ID %in% c("F", "G")), aes(x= site_ID, y= mean.even, fill = site_ID, color =site_ID))+
  geom_boxplot()+ theme(legend.position = "none")+
  scale_color_manual(values= c(site.colors.dark, "grey22"))+
  scale_fill_manual(values= adjustcolor(c(site.colors.dark[1:5],"grey22"), alpha.f = 0.2))+
  ylab("Fungal evenness") + xlab("Composting site")+
  bg_theme +theme(axis.text.x = element_text(angle = 0, hjust=0.3, vjust=0))+
  geom_text(data =value_max, aes(x=site_ID, y = quantile, label = letter), hjust=-0.2, vjust =-0.5,  size= 6, colour ="black")

# Figure 3D
# Statistics
kruskal_result <- kruskal.test(mean.sobs ~ batch, data = df_BLW2.div)
pairwise_results = pairwise.wilcox.test(df_BLW2.div$mean.sob, df_BLW2.div$batch, p.adj ="BH")
#comparison_letters <- multcompLetters(as.vector(pairwise_results$p.value))
KRUSKAL =tri.to.squ(pairwise_results$p.value)
compact_letters_df =multcompLetters(KRUSKAL)

value_max = df_BLW2.div %>% group_by(batch) %>% dplyr::summarize(max_value = max(mean.sobs),
                                                                 quantile = quantile(mean.sobs, probs = 0.75))
value_max$batch <- as.factor(value_max$batch)
value_max$letter <- c("a", "b", "b", "ab")

# Plot wiht letters
g4 = ggplot(df_BLW2.div , aes(x= batch, y= mean.sobs, fill = batch))+
  geom_boxplot()+
  scale_fill_manual(values =c (rep("white",4)))+
  scale_x_discrete(labels= c("I", "II", "III", "IV"))+
  ylab("") + xlab("")+
  bg_theme + theme(axis.text.x = element_text(angle = 0, hjust=0.3, vjust=0))+
  geom_text(data =value_max, aes(x=batch, y = 20 + quantile, label = letter), hjust=-0.2, vjust =-0.5,  size= 6, colour ="black")

g = plot_grid(g2, g3, g4, nrow =1, axix ="none", ncol =3, rel_widths = c(2,2,2),labels =c("B", "C", "D")) 
plot_grid(g1, g, nrow=2, ncol =1, rel_heights = c(1.1,1), labels =c("A", ""))

#ggsave(filename =paste0("figures/",path, "F_Alpha_diversity_plot_A_D.png"), height =9, width =14)
g = plot_grid(g5, g3, nrow =1, axis ="non", ncol =2, rel_widths = c(1,0.4), labels = c("E", "F"))
ggsave(g, filename ="Figures/F_Alpha_diversity_plot_E_F.png", height =4.5, width =7)

rm(g1, g2,g3, g4, data, data1, g, cor1, cor2, cor3, kruskal_result, pairwise_results, KRUSKAL, compact_letters,
   compact_letters_df, color1, color2, color3)

### Tbl S7_CP_BC------
# Correlation analysis compost properties + compost age and alpha diversity metrics

data <-df_BLW2.div
factor.comp = c( "mean.sobs", "mean.even", "mean.shannon", "mean.ivsimp")
factor.list = c("age",factor.continous)
correlation_estimates <- data.frame(Variable1 = character(0), Variable2 = character(0),
                                    Estimate = numeric(0),
                                    P_Value = numeric(0),
                                    P_adj = numeric(0),
                                    n = numeric(0))
correlation_matrix <- matrix(NA, nrow = length(factor.list), ncol = length(factor.comp))

for (i in 1:length(factor.list)){
  for (j in 1:length(factor.comp)) {
    var1 <- factor.list[i]
    var2 <- factor.comp[j]
    
    cor <-cor.test(data[[var1]], data[[var2]], method = "spearman", use= "p")
    n <- data[c(var1, var2)] %>% na.omit() %>% nrow()
    estimate <- cor$estimate
    p_value <- cor$p.value
    correlation_estimates <- rbind(correlation_estimates, data.frame(Variable1 = var1, Variable2 = var2,
                                                                     Estimate = estimate,
                                                                     P_value = p_value, n = n))
    if (p_value < 0.05) {
      # Calculate the correlation and store it in the appropriate position
      correlation <- cor(data[[var1]], data[[var2]], method = "spearman", use= "p")
      correlation_matrix[i, j] <- correlation
    }
  }
}
rownames(correlation_matrix) = factor.list; colnames(correlation_matrix) = factor.comp
correlation_matrix[is.na(correlation_matrix)] <-0
corrplot(correlation_matrix, method = "number", diag= TRUE)

rownames(correlation_estimates) = NULL
correlation_estimates[,3] = round(correlation_estimates[,3], 2)
correlation_estimates[,4] = round(correlation_estimates[,4], 3)
correlation_estimates <-correlation_estimates %>% arrange(Variable2)

# Correction for multiple testing

list <-list()
for (i in factor.comp) {
  list[[i]] <-correlation_estimates %>% filter(Variable2 == i)  %>%
    mutate(p_adj = p.adjust(P_value, "BH") %>% round(3))
}
correlation_estimates <-do.call(rbind, list)

rownames(correlation_estimates) <- NULL
#write.csv(correlation_estimates, file = paste0("Output/", path, "F_Alpha_diversity_compost_prop_multiple_testing_filled.csv"), row.names = FALSE)

# Beta diversity

data <- df_BLW2.r %>% filter(treatment %in% df_BLW2.div$treatment)
factor.list = c("age",factor.continous)
correlation_estimates <- matrix(nrow = length(factor.list), ncol =4) %>% as.data.frame()
colnames(correlation_estimates) <- c("variable", "R2", "Fvalue", "P_Value")

for (i in 1:length(factor.list)) {
  data$x <-data[,factor.list[i]]
  data_red <- data[!is.na(data$x), ] 
  ISS.bray.avg.comp.red <-  ISS.bray.avg.comp[rownames(ISS.bray.avg.comp) %in% data_red$treatment ,
                                              colnames(ISS.bray.avg.comp) %in% data_red$treatment ]
  adonis_results <-adonis2(ISS.bray.avg.comp.red ~ x, data = data_red, permutations = 999)
  correlation_estimates$variable[i] <- factor.list[i]
  correlation_estimates$R2[i] <-adonis_results$R2[1] %>% round(3)
  correlation_estimates$Fvalue[i] <- adonis_results$F[1] %>% round(2)
  correlation_estimates$P_Value[i] <- adonis_results$`Pr(>F)`[1] %>% round(3)
}
# Correct for multiple testing

correlation_estimates <-correlation_estimates %>% mutate(
  p_adj = p.adjust(P_Value, method ="BH") %>% round(3)
)
correlation_estimates$R2 <- correlation_estimates$R2*100

write.csv(correlation_estimates, file =paste0("Output/", path, "F_Beta_diversity_compost_prop_multiple_testing_filled.csv"), row.names = FALSE)

## Beta diversity-----

data <-ISS.bray.avg.comp
#data <- ISS.bray.avg.comp[!rownames(ISS.bray.avg.comp) %in% c("K30"), !colnames(ISS.bray.avg.comp) %in% c("K30")]
nmds.BLW2 <- metaMDS(data)
nmds.BLW2.points = nmds.BLW2$points %>% as.data.frame()
colnames(nmds.BLW2.points) = c("nmds1", "nmds2")
nmds.BLW2.points$treatment = rownames(nmds.BLW2.points)
df_BLW2.beta <- merge(df_BLW2.r, nmds.BLW2.points , by ="treatment")
nmds.BLW2$stress # 0.17 10000, 0.14 20000

# composting site

find_hull <- function(df_BLW2.beta) df_BLW2.beta[chull(df_BLW2.beta$nmds1, df_BLW2.beta$nmds2), ]
hulls <- plyr::ddply(df_BLW2.beta, "site_ID", find_hull)

g1 =ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2))+
  geom_point(aes(color = site_ID))+
  theme_classic(base_size = 14) + theme(legend.position = "none") +
  scale_fill_manual(values = site.colors.dark, name = "Composting\nsite")+
  scale_color_manual(values = site.colors.dark, name = "Composting\nsite")+
  geom_text(x=0.28,y=-0.28,label='stress = 0.18', size=4)+
  geom_polygon(data=hulls,alpha=0.2,aes(fill=site_ID))+
  ylab("NMDS2")+xlab("NMDS1")+
  coord_fixed(ratio =1)+
  geom_text_repel(data = df_BLW2.beta, aes(x = nmds1 , y = nmds2 , label = Number, color = site_ID), show.legend = FALSE, max.overlaps =20)

#ggsave(g1, filename = paste0("figures/", path, "F_NMDs_Betadiv_composting_site_with_20.png"), height = 3, width =8)
#ggsave(g1, filename = paste0("figures/", path, "F_NMDs_Betadiv_composting_site.png"), height = 8, width =8)


# Disease suppression
color = c("orange", "firebrick","orchid4","blue3","cornflowerblue")
color = c( "#F1F1F1","#FFC59E" ,"#E1BB4E", "yellow", "#9BB306", "#26A63A", "darkgreen")
color = c("darkred", "yellow")
 

g2 = ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cp))+
  theme_classic(base_size = 14) +
  geom_point(aes(fill = rel.cp), size =4, pch=21, colour ="black")+
  coord_fixed(ratio =1) +
  theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
  scale_fill_gradientn(colors= color, limits = c(0,110), name= "Disease\n suppression [%]")

g3= ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cu.gu))+
  theme_classic(base_size = 14) +
  geom_point(aes(fill = rel.cu.gu), size =4, pch=21, colour="black")+
  coord_fixed(ratio =1) +
  theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
  scale_fill_gradientn(colours = color,
                       limits =c(0,110), name= "Disease\nsuppression [%]")

g4= ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cu.rs))+
  theme_classic(base_size = 14) +
  geom_point(aes(fill = rel.cu.rs), size =4, pch =21, colour="black")+
  coord_fixed(ratio =1) + 
  theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
  scale_fill_gradientn(colours = color,
                       limits = c(0,110),name= "Disease\nsuppression [%]")

legend1 = ggpubr::get_legend(g2+ theme(legend.position = "bottom"))
legend2 = ggpubr::get_legend(g1+ theme(legend.position = "bottom"))

plots_without_legends <- plot_grid(g2, g3, g4, g1, ncol = 2, align = 'hv', labels = c("A", "B", "C", "D"))
legends <- plot_grid(legend1, legend2, ncol = 2, align ='v')
final_plot <- plot_grid(plots_without_legends, legends, ncol = 1, align = 'v', rel_heights = c(1, 0.1))
ggsave(final_plot, height =10, width =10, file =paste0("figures/", path, "F_NMDS_Betadiv_FigureA-D.png"))


# Plot only the three disease assays
ggarrange(g2, g3, g4, common.legend=T, nrow =1, ncol =3, labels  =c ("G. ultimum -cress", "G. ultimum -cucumber", "R. solani - cucumber"))
#ggsave(file ="figures/F_NMDS_Betadiv_FigureA-C_light.png", height = 4, width =12)

# PERMANOVA

# ATTENTION: Checkif correct df_BLW2.beta is loaded!

#sink(file =paste0("output/", path, "PERMANOVA_beta_div_without20.txt"))
paste("Composting site")
data <- ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K30", "K48", "K49")),
                          !(colnames(ISS.bray.avg.comp) %in% c("K30", "K48", "K49"))]
adonis2(data ~ site, data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)

paste("composting site + disease suppression")
adonis2(data ~ site + rel.cp  , data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)
adonis2(data ~ site + rel.cu.gu, data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)
adonis2(data ~ site +rel.cu.rs, data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)

paste("batch")
data <- ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K30")),
                          !(colnames(ISS.bray.avg.comp) %in% c("K30"))]
adonis2(data ~ batch, data = df_BLW2.beta, permutations = 999)

paste("Compost age")
data = ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K43", "K44", "K30")),
                         !(colnames(ISS.bray.avg.comp) %in% c("K43", "K44", "K30"))]
adonis2(data  ~ age, data = df_BLW2.beta %>% filter(!treatment %in% c("K43", "K44", "K30")), permutations = 999)

paste("Disease suppression")

data <- ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K30")),
                          !(colnames(ISS.bray.avg.comp) %in% c("K30"))]
adonis2(data ~ rel.cp  , data = df_BLW2.beta, permutations = 999)
adonis2(data ~ rel.cu.gu, data = df_BLW2.beta, permutations = 999)
adonis2(data ~ rel.cu.rs, data = df_BLW2.beta , permutations = 999)

#sink()

sink(file =paste0("output/", path, "PERMANOVA_beta_div.txt"))
paste("Composting site")
data <- ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K48", "K49")),
                          !(colnames(ISS.bray.avg.comp) %in% c("K48", "K49"))]
adonis2(data ~ site, data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)

paste("composting site + disease suppression")
adonis2(data ~ site + rel.cp  , data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)
adonis2(data ~ site + rel.cu.gu, data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)
adonis2(data ~ site +rel.cu.rs, data = df_BLW2.beta[!df_BLW2.beta$site %in% c("Frick", "Spreitenbach"),], permutations = 999)

paste("batch")
data <- ISS.bray.avg.comp
adonis2(data ~ batch, data = df_BLW2.beta, permutations = 999)

paste("Compost age")
data = ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K43", "K44" )),
                         !(colnames(ISS.bray.avg.comp) %in% c("K43", "K44"))]
adonis2(data  ~ age, data = df_BLW2.beta %>% filter(!treatment %in% c("K43", "K44")), permutations = 999)

paste("Disease suppression")

data <- ISS.bray.avg.comp
adonis2(data ~ rel.cp  , data = df_BLW2.beta, permutations = 999)
adonis2(data ~ rel.cu.gu, data = df_BLW2.beta, permutations = 999)
adonis2(data ~ rel.cu.rs, data = df_BLW2.beta, permutations = 999)

sink()

rm(data, final_plot, g1, g2, g3, g4, hulls, legend1, legend2, legends, nmds.BLW2.points, 
   nmds.BLW2, plots_without_legends, color)
## S13 Sup_F_tax-------

# Aggregate Phylum level
ISS.p <- aggregate(t(ISS.prop) ~ Phyla, data = tax.class.red[colnames(ISS.prop),], FUN = sum)  # aggregate by phylum
ISS.p <-ISS.p %>% relocate(Phyla, K7, K8, K9)
ISS.p <- as.matrix(data.frame(ISS.p[, -1], row.names = ISS.p[, 1], check.names = FALSE))  # move first col to row.names
ISS.p <- ISS.p[order(rowSums(ISS.p), decreasing = TRUE), ]  # order by rowSums
ISS.p.nc <- as.data.frame(ISS.p[rownames(ISS.p) != "unclassified", ])  #remove unclassified group
ISS.p.nc["unclassified", ] <- ISS.p[rownames(ISS.p) == "unclassified", drop = F]  # add unclassified at the end
ISS.p.plot <- as.matrix(ISS.p.nc)

# Aggregate by family level
ISS.f <- aggregate(t(ISS.prop) ~ Family, data = tax.class.red[colnames(ISS.prop),], FUN = sum)  # aggregate by family
ISS.f <-ISS.f %>% relocate(Family, K7, K8, K9)
ISS.f <- as.matrix(data.frame(ISS.f[, -1], row.names = ISS.f[, 1], check.names = FALSE))  # move first col to row.names
ISS.f <- ISS.f[order(rowSums(ISS.f), decreasing = TRUE), ]  # order by rowSums
ISS.f.nc <- as.data.frame(ISS.f[rownames(ISS.f) != "unclassified", ])  # remove unclassified group
ISS.f.nc["others", ] <- colSums(ISS.f.nc[min(nrow(ISS.f.nc), 21):nrow(ISS.f.nc), ])  # merge rare groups to others
ISS.f.nc["unclassified", ] <- ISS.f[rownames(ISS.f) == "unclassified", drop = F]  # add unclassified at the end
ISS.f.plot <- ISS.f.nc[-(min(nrow(ISS.f.nc), 21):(nrow(ISS.f.nc) - 2)), ]  # remove rare groups
ISS.f.plot <- as.matrix(ISS.f.plot)


# Aggreate by genera level

ISS.g <- aggregate(t(ISS.prop) ~ Genus, data = tax.class.red[colnames(ISS.prop),], FUN = sum)  # aggregate by family
ISS.g <-ISS.g %>% relocate(Genus, K7, K8, K9)
ISS.g <- as.matrix(data.frame(ISS.g[, -1], row.names = ISS.g[, 1], check.names = FALSE))  # move first col to row.names
ISS.g <- ISS.g[order(rowSums(ISS.g), decreasing = TRUE), ]  # order by rowSums
ISS.g.nc <- as.data.frame(ISS.g[rownames(ISS.g) != "unclassified", ])  # remove unclassified group
ISS.g.nc["others", ] <- colSums(ISS.g.nc[min(nrow(ISS.g.nc), 21):nrow(ISS.g.nc), ])  # merge rare groups to others
ISS.g.nc["unclassified", ] <- ISS.g[rownames(ISS.g) == "unclassified", drop = F]  # add unclassified at the end
ISS.g.plot <- ISS.g.nc[-(min(nrow(ISS.g.nc), 21):(nrow(ISS.g.nc) - 2)), ]  # remove rare groups
ISS.g.plot <- as.matrix(ISS.g.plot)


custom_x_labels <- compost_ID %>% dplyr::filter(!treatment %in% c("K6", "K10", "K17", "K27", "K37", "Std2",
                                                                  "Std3", "Std4", "Std5" )) %>% 
  dplyr::select(compost_ID) 

# Plot 1: Phyla

# Calculate mean percentages for each phylum
mean_percentages <- rowMeans(ISS.p.plot)
# Create labels with mean percentages
legend_labels <- paste(rownames(ISS.p.plot), sprintf("(%.3f%%)", mean_percentages))

#pdf(file =paste0(output_path, "Sup_F_tax_A.pdf"), width=10, height=5)
par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2))
color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
           "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#000000")
barplot(ISS.p.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
        width = 1, cex.axis = 0.8, cex.names = 0.8, axisnames = TRUE, names.arg = custom_x_labels$compost_ID)
par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.2, 0.2), new = TRUE)
plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topright", legend = legend_labels, pch = 22, cex = 0.7, bty = "n",
       ncol = 1, pt.bg = color, col = NA)
#dev.off()


# Plot 2: Family level

# Calculate mean percentages for each phylum
mean_percentages <- rowMeans(ISS.f.plot)
# Create labels with mean percentages
legend_labels <- paste(rownames(ISS.f.plot), sprintf("(%.2f%%)", mean_percentages))

#pdf(file =paste0(output_path, "Sup_F_tax_B.pdf"), width=10, height=5)
par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2))
color <-  c(replicate(20, generate_random_color()), "#999999", "#000000")
barplot(ISS.f.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
        width = 1, cex.axis = 0.8, cex.names = 0.8, axisnames = TRUE, names.arg = custom_x_labels$compost_ID)
par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.2, 0.2), new = TRUE)
plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topright", legend = legend_labels, pch = 22, cex = 0.7, bty = "n",
       ncol = 1, pt.bg = color, col = NA)

#dev.off()

# Plot 3: Genera level

# Calculate mean percentages for each phylum
mean_percentages <- rowMeans(ISS.g.plot)
# Create labels with mean percentages
legend_labels <- paste(rownames(ISS.g.plot), sprintf("(%.2f%%)", mean_percentages))
color <- c("#8998C8" ,"#2A2772","#6CF01B","#68322C","#6AFC75","#ADDAD8","#B887C3","#9F847C", "#C318C3",
           "#791B8C", "#0E48F8", "#E8991C", "#9C7186", "#E003D3", "#FAA1C3" ,"#33E8C0" ,"#0596C2" ,
           "#16B973", "#8BB0D2", "#014BAD", "#999999" ,"#000000")
#pdf(file =paste0(output_path, "Sup_F_tax_C.pdf"), width=10, height=5)
par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2))
#color <-  c(replicate(20, generate_random_color()), "#999999", "#000000")
barplot(ISS.g.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
        width = 1, cex.axis = 0.8, cex.names = 0.8, axisnames = TRUE, names.arg = custom_x_labels$compost_ID)
par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.2, 0.2), new = TRUE)
plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topright", legend = legend_labels, pch = 22, cex = 0.7, bty = "n",
       ncol = 1, pt.bg = color, col = NA)

#dev.off()

rm(mean_percentages, color, legend_labels, ISS.g, ISS.p, ISS.f,
   ISS.g.nc, ISS.p.nc, ISS.f.nc, ISS.g.plot, ISS.f.plot, ISS.p.plot, custom_x_labels)

## Sup. 4 Beta diversity all samples----

nmds.BLW2 <- metaMDS(ISS.bray)
nmds.BLW2.points = nmds.BLW2$points %>% as.data.frame()
colnames(nmds.BLW2.points) = c("nmds1", "nmds2")
df_BLW2.beta <- merge(design.BLW2, nmds.BLW2.points , by =0)
rownames(df_BLW2.beta) <- df_BLW2.beta$Row.names; df_BLW2.beta$Row.names <-NULL
nmds.BLW2$stress # 0.075 10000, 0.070 20000

# Plot

df_BLW2.beta$compost_ID <-droplevels(df_BLW2.beta$compost_ID)
levels <- compost_ID%>% filter(!treatment %in% compost.exclude) %>% pull(compost_ID) %>% as.character()
df_BLW2.beta$compost_ID<- factor(df_BLW2.beta$compost_ID, levels=levels)

levels <- compost_ID %>% pull(compost_ID) %>% as.character()
compost_ID$compost_ID<- factor(compost_ID$compost_ID, levels=levels)
rm(levels)

find_hull <- function(df_BLW2.beta) df_BLW2.beta[chull(df_BLW2.beta$nmds1, df_BLW2.beta$nmds2), ]
hulls <- plyr::ddply(df_BLW2.beta, "compost_ID", find_hull)

ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2))+
  geom_point(aes(color = compost_ID))+
  theme_classic(base_size = 14) + theme(legend.position = "right") +
  geom_text(x=0.6,y=-0.1, label='stress = 0.08', size=4)+
  ylab("NMDS2")+xlab("NMDS1")+
  geom_polygon(data=hulls,alpha=0.2,aes(fill=compost_ID))+
  coord_fixed(ratio =1)

#ggsave(filename = paste0("figures/", path, "F_Sup4_Beta_diverstiy_all_samples.png"), height = 7, width =12)

adonis2(ISS.bray[!grepl("Std", rownames(ISS.bray)), 
                 !grepl("Std", colnames(ISS.bray))] ~ treatment, 
        data = df_BLW2.beta %>% filter(!treatment %in% c("Std2", "Std3", "Std4", "Std5")),
        permutations = 999)

## Indicator species analysis------

### Summarize on taxonomic level----

path1 =paste0("output/", path, "Indicator_species_analysis/")

# ATTENTION: Only load one of them!
# Analysis based on ASV level 3/4 replicates
data.ISS = ISS.rob.comp.avg  %>% as.data.frame()
data.asv = asv.rob.comp.avg  %>% t() %>% as.data.frame()

# Analysis based on ASV level 2/4 replicates
#data.ISS = ISS.rob.comp.avg2  %>% as.data.frame()
#data.asv = asv.rob.comp.avg2  %>% t() %>% as.data.frame()

# Analysis based on genus level
data.tax = tax.class.red[rownames(data.ISS), ]
name = colnames(tax.class.red[2:ncol(tax.class.red)]) %>% rev()
summary =matrix(nrow = length(name), ncol =10) %>% as.data.frame()
colnames(summary) = c("tax level", "nclassfiedASVs", "perclassfiedASVS", "topnASVs1", "topnASVs2", "topnASVs3", "nclassications",
                      "nreadtop1", "nreadstop2", "nreadstop3")
list = list()

# ATTENTION: Either select ISS or asv! ISS is the default

# The following files are needed: data.tax, ISS.rob.comp.avg| asv.rob.comp.avg
# Only robustly detected ASVs are considered in the analysis

for (i in 1:length(name)) {
  data.tax$level = data.tax[, name[i]]
  tax.rob.id.genus =data.tax %>% filter(!grepl("unclassified", level))
  summary[i, 1] = name[i]
  summary[i, 2] = nrow(tax.rob.id.genus) # number of classified ASVs
  summary[i, 3] = 1-(nrow(data.tax) - nrow(tax.rob.id.genus)) / nrow(data.tax) # percentage of ASVs that are classified on the level
  summary[i,4:6] =table(tax.rob.id.genus$level) %>% as.data.frame() %>% arrange(desc(Freq)) %>% as.data.frame() %>% head(3) %>% pull(Var1) %>% as.vector() # Classifcation with most ASVs
  summary[i, 7] = nrow(table(tax.rob.id.genus$level) %>% as.data.frame() %>% arrange(desc(Freq))) # Number of classifcations on the level
  #ISS.rob.av.comp.id.genus = data.ISS[rownames(tax.rob.id.genus),] # ISS
  ISS.rob.av.comp.id.genus = data.asv[rownames(tax.rob.id.genus),] # ASV
  ISS.rob.av.comp.id.genus = merge(tax.rob.id.genus, ISS.rob.av.comp.id.genus, by = 0)
  rownames(ISS.rob.av.comp.id.genus) <- ISS.rob.av.comp.id.genus$Row.names
  ISS.rob.av.comp.id.genus$Row.names <- NULL
  result <- aggregate(. ~ level, data = ISS.rob.av.comp.id.genus[,8:ncol(ISS.rob.av.comp.id.genus)], sum)
  rownames(result) <- result$level
  result$level <- NULL
  summary[i, 8:10] = rowSums(result) %>% sort() %>% as.data.frame() %>% tail(3) %>% rownames() # Classication with the highest number of reads
  list[[i]] = result
}

# Save summary
# write.csv(summary, file = paste0(path1, "F_BLW2_summary_taxonomy_ASV.csv" ))
# write.csv(summary, file = paste0(path1, "F_BLW2_summary_taxonomy_ASV_2.csv" ))


# ISS
data.ISS.species = list[[1]] 
data.ISS.genus = list[[2]] 
data.ISS.family = list[[3]]
data.ISS.phyla <- list[[6]]

# ASV
data.asv.species = list[[1]]
data.asv.genus = list[[2]]
data.asv.family = list[[3]]
data.asv.phyla <- list[[6]]

list = list(nASV = NULL, nBLW1 = NULL, min2 = NULL, min3 = NULL, all = NULL)

## Gu-cress-----

top = df_BLW2.r %>% filter (treatment %in% c("K8", "K36", "K34", "K15", "K28", "K13", "K12", "K30", "K35"))
flop = df_BLW2.r %>% filter(treatment %in% c("K23", "K26", "K21", "K20", "K18", "K19", "K9", "K48", "K7"))
#flop = df_BLW2.r %>% filter(treatment %in% c("K23", "K26", "K21", "K20", "K18", "K9", "K48", "K7")) # Without K19
path.plant = "cress"
robdetect<-"_2_4"
#robdetect <-"_3_4"

# ASV level
level <- "ASV"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cress = ISA_PALMA(top, flop, data.ISS, data.asv, path2, path.plant = "cress", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cress, level =level, cutoff = TRUE)
# Without cut
output_cress = ISA_PALMA(top, flop, data.ISS, data.asv, path2, path.plant = "cress", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cress, level =level, cutoff = FALSE)

# Genera level
level <- "genus"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path2, path.plant = "cress",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cress, level =level, cutoff = TRUE)
# Without cut
output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path2, path.plant = "cress",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cress, level =level, cutoff = FALSE)

# Familiy level
level <- "family"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.family, data.asv.family, path2, path.plant = "cress",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cress, level = level, cutoff = TRUE)
# Without cut
output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.family, data.asv.family, path2, path.plant = "cress",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cress, level = level, cutoff = FALSE)

# Phyla
level <- "phyla"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.phyla, data.asv.phyla, path2, path.plant = "cress",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cress, level = level, cutoff = TRUE)
# Without cut
output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.phyla, data.asv.phyla, path2, path.plant = "cress",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cress, level = level, cutoff = FALSE)


## Gu-cuc-----
top = df_BLW2.r %>% filter (treatment %in% c("K7", "K9", "K12", "K13", "K15", "K28", "K34", "K45", "K47"))
flop = df_BLW2.r %>% filter(treatment %in% c("K18", "K20", "K21", "K22", "K23", "K25", "K26", "K43", "K49"))
path.plant = "cugu"
robdetect<-"_2_4"
#robdetect <-"_3_4"

# ASV level
level <- "ASV"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cugu = ISA_PALMA(top, flop, data.ISS, data.asv, path2, path.plant = "cugu", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cugu, level =level, cutoff = TRUE)
# Without cut
output_cugu = ISA_PALMA(top, flop, data.ISS, data.asv, path2, path.plant = "cugu", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cugu, level =level, cutoff = FALSE)

# Genera level
level <- "genus"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cugu = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path2, path.plant = "cugu",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cugu, level =level, cutoff = TRUE)
# Without cut
output_cugu = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path2, path.plant = "cugu",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cugu, level =level, cutoff = FALSE)

# Familiy level
level <- "family"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cugu = ISA_PALMA(top, flop, data.ISS = data.ISS.family, data.asv.family, path2, path.plant = "cugu",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cugu, level = level, cutoff = TRUE)
# Without cut
output_cugu = ISA_PALMA(top, flop, data.ISS = data.ISS.family, data.asv.family, path2, path.plant = "cugu",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cugu, level = level, cutoff = FALSE)

# Phyla
level <- "phyla"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_cugu = ISA_PALMA(top, flop, data.ISS = data.ISS.phyla, data.asv.phyla, path2, path.plant = "cugu",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_cugu, level = level, cutoff = TRUE)
# Without cut
output_cugu = ISA_PALMA(top, flop, data.ISS = data.ISS.phyla, data.asv.phyla, path2, path.plant = "cugu",
                         fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_cugu, level = level, cutoff = FALSE)

## Rs-cuc-------
top = df_BLW2.r %>% filter (treatment %in% c("K34", "K29", "K23", "K26", "K36", "K44", "K46", "K47", "K48"))
#flop = df_BLW2.r %>% filter(treatment %in% c("K8", "K12", "K20", "K21", "K22", "K30", "K33", "K45"))
flop = df_BLW2.r %>% filter(treatment %in% c("K8", "K12", "K19", "K20", "K21", "K22", "K30", "K33", "K45"))
path.plant = "curs"
#robdetect<-"_2_4"
robdetect <-"_3_4"

# ASV level
level <- "ASV"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_curs = ISA_PALMA(top, flop, data.ISS, data.asv, path2, path.plant = "curs", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_curs, level =level, cutoff = TRUE)
# Without cut
output_curs = ISA_PALMA(top, flop, data.ISS, data.asv, path2, path.plant = "curs", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_curs, level =level, cutoff = FALSE)

# Genera level
level <- "genus"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_curs = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path2, path.plant = "curs",
                        fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_curs, level =level, cutoff = TRUE)
# Without cut
output_curs = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path2, path.plant = "curs",
                        fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_curs, level =level, cutoff = FALSE)

# Familiy level
level <- "family"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_curs = ISA_PALMA(top, flop, data.ISS = data.ISS.family, data.asv.family, path2, path.plant = "curs",
                        fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_curs, level = level, cutoff = TRUE)
# Without cut
output_curs = ISA_PALMA(top, flop, data.ISS = data.ISS.family, data.asv.family, path2, path.plant = "curs",
                        fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_curs, level = level, cutoff = FALSE)

# Phyla
level <- "phyla"
path2 <- paste0(path1, level, robdetect, "/")
# With cut
output_curs = ISA_PALMA(top, flop, data.ISS = data.ISS.phyla, data.asv.phyla, path2, path.plant = "curs",
                        fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) 
comparison_Indicator_analysis(output_curs, level = level, cutoff = TRUE)
# Without cut
output_curs = ISA_PALMA(top, flop, data.ISS = data.ISS.phyla, data.asv.phyla, path2, path.plant = "curs",
                        fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) 
comparison_Indicator_analysis(output_curs, level = level, cutoff = FALSE)

## Comparison among pathogen-plant systems----

# Attention phylogenetical tree only works for ASVs!

comparison_pathogen_plant(level ="ASV", robdetect = "_3_4", cutoff = TRUE, nmethods =3)
comparison_pathogen_plant(level ="ASV", robdetect = "_3_4", cutoff = TRUE, nmethods =2)

comparison_pathogen_plant(level ="ASV", robdetect = "_3_4", cutoff = FALSE, nmethod =3)
comparison_pathogen_plant(level ="ASV", robdetect = "_3_4", cutoff = FALSE, nmethods =2)

comparison_pathogen_plant(level ="genus", robdetect = "_3_4", cutoff = TRUE)
comparison_pathogen_plant(level ="genus", robdetect = "_3_4", cutoff = FALSE)

comparison_pathogen_plant(level ="family", robdetect = "_3_4", cutoff = TRUE)
comparison_pathogen_plant(level ="family", robdetect = "_3_4", cutoff = FALSE)

comparison_pathogen_plant(level ="phyla", robdetect = "_3_4", cutoff = TRUE)
comparison_pathogen_plant(level ="phyla", robdetect = "_3_4", cutoff = FALSE)

## Phylogenetic tree all-------

# Not necessary since we do not have ASVs that were indicative in all three systems!!!
#ASV_present_all <-all %>% filter(assay_combination == "cress_cugu_curs") %>% dplyr::select(ASV) 
#ASV_all = all %>% filter(ASV %in% ASV_present_all$ASV & assay_combination == "cress_cugu_curs")
#all =all %>% filter(!ASV %in% ASV_present_all$ASV)
#all = rbind(all, ASV_all)

# Add the classification identification again
level <- "ASV"
cutoff =TRUE
cutoff =FALSE
all <- read.csv(file =paste0(path1, level,robdetect, "/comparisons",
                             ifelse(cutoff == TRUE, "_cut", ""), "_path_plants_systems.csv")) %>%
  filter(comparison ==1)



# Define prefixes
prefixes <- c(Class = "c_", Order = "o_", Family = "f_", Genus = "g_", Species = "s_")

# Apply prefixes conditionally: only if the entry does not contain "unclassified"
all <- all %>%
  mutate(across(all_of(names(prefixes)), 
                ~ ifelse(grepl("unclassified", ., ignore.case = TRUE), ., paste0(prefixes[cur_column()], .))))

colors <- data.frame(color = c("brown", "chartreuse", "darkgoldenrod", "darkseagreen", "darkslategrey", "darkorchid", "chocolate", "aquamarine2",
                               "bisque4", "cadetblue", "coral","darkolivegreen", "deeppink"),
                    phyla = tax.class.red$Phyla %>% unique())

clade_colors = colors %>% filter(phyla %in% unique(all$Phyla)) %>% dplyr::select(color)
phylo_tree <- phylo_tree_ASV_flat(all, sequences, clade_colors, all =TRUE, size =4, offset = 0.2, xlim.factor = 2)

# Mit GU-cress/RS-cuc§
#ggsave("Figures/All_phylogenetical_tree.png", height = 8, width =10, unit ="in")

# Without GU-cress/RS-cuc
all_red <- all %>% filter(assay_combination %in% c("cress_cugu", "cugu_curs", "cress_cugu_curs"))
clade_colors = colors %>% filter(phyla %in% unique(all_red$Phyla)) %>% dplyr::select(color)
phylo_tree_all <- phylo_tree_ASV_flat(all_red, sequences, clade_colors, all=TRUE, size =3.5, xlim.factor = 1.8)

p1 <- phylo_tree_all[[2]]
p2 <-ggarrange(p1, phylo_tree[[3]], ncol =2, nrow =1, widths = c(0.85, 0.15))
#ggsave(p2, filename= paste0("figures/", path, "Indicator_species_analysis/", level, robdetect, "/All_phylogenetical_tree_",
#                            ifelse(cutoff ==TRUE, "_cut", ""), ".png"), height =5, width =6, unit ="in")

## Phylogenetic tree seperatly------
# Add the classification identification again
level <- "ASV"
cutoff =TRUE
robdetect <- "_3_4"
#cutoff =FALSE
cress <- read.csv(file =paste0(path1, level,robdetect, "/cress",
                             ifelse(cutoff == TRUE, "_cut", ""), "_method_comparison.csv")) %>%
  filter(Freq > 1)
rownames(cress) <- cress$X; cress$X <-NULL
cress$ASV <- rownames(cress)

cugu <- read.csv(file =paste0(path1, level,robdetect, "/cugu",
                               ifelse(cutoff == TRUE, "_cut", ""), "_method_comparison.csv")) %>%
  filter(Freq > 1)
rownames(cugu) <- cugu$X; cugu$X <-NULL

curs <- read.csv(file =paste0(path1, level,robdetect, "/curs",
                               ifelse(cutoff == TRUE, "_cut", ""), "_method_comparison.csv"), sep =";") %>%
  filter(Freq > 1)
rownames(curs) <- curs$X; curs$X <-NULL


# Define prefixes
prefixes <- c(Class = "c_", Order = "o_", Family = "f_", Genus = "g_", Species = "s_")

# Apply prefixes conditionally: only if the entry does not contain "unclassified"
cress <- cress %>%
  mutate(across(all_of(names(prefixes)), 
                ~ ifelse(grepl("unclassified", ., ignore.case = TRUE), ., paste0(prefixes[cur_column()], .))))

colors <- data.frame(color = c("brown", "chartreuse", "darkgoldenrod", "darkseagreen", "darkslategrey", "darkorchid", "chocolate", "aquamarine2",
                               "bisque4", "cadetblue", "coral","darkolivegreen", "deeppink"),
                     phyla = tax.class.red$Phyla %>% unique())

clade_colors = colors %>% filter(phyla %in% unique(cress$Phyla)) %>% dplyr::select(color)
phylo_tree <- phylo_tree_ASV_flat(cress, sequences, clade_colors, all =FALSE, size =4, offset = 0.2, xlim.factor = 2)

for (i in 1:length(cress$ASV)) {
  data <- ISS.rob.comp.avg[cress$ASV[i],] %>% t() %>% as.data.frame()
  data$treatment <-rownames(data)
  colnames(data) <- c("ASV", "treatment")
  data <-merge(df_BLW2.r, data, by ="treatment")
  cor <-cor.test(data$rel.cp, data$ASV, method = "spearman")
  print(cor)
}




## Heat map and phylogenetical tree-------
asv_names <- all$ASV
ds <- c("rel.cp", "rel.cu.gu", "rel.cu.rs")
tip_order <-phylo_tree_all[[2]]$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

topflop <- list()
topflop[["T9_cress"]] <-df_BLW2.r %>% filter (treatment %in% c("K8", "K36", "K34", "K15", "K28", "K13", "K12", "K30", "K35"))
topflop[["F9_cress"]] <-df_BLW2.r %>% filter(treatment %in% c("K23", "K26", "K21", "K20", "K18", "K19", "K9", "K48", "K7"))
topflop[["T9_cugu"]] <-df_BLW2.r %>% filter (treatment %in% c("K7", "K9", "K12", "K13", "K15", "K28", "K34", "K45", "K47"))
topflop[["F9_cugu"]] <-df_BLW2.r %>% filter(treatment %in% c("K18", "K20", "K21", "K22", "K23", "K25", "K26", "K43", "K49"))
topflop[["T9_curs"]] <-df_BLW2.r %>% filter (treatment %in% c("K14", "K15", "K23", "K26", "K36", "K44", "K46", "K47", "K48"))
topflop[["F9_curs"]] <-df_BLW2.r %>% filter(treatment %in% c("K8", "K12", "K19", "K20", "K21", "K22", "K30", "K33", "K45"))
topflop[["T3"]] <-df_BLW2.r %>% filter(treatment %in% c("K34", "K46", "K47")) 
topflop[["F4"]] <-df_BLW2.r %>% filter(treatment %in% c("K18", "K19", "K20", "K21")) 

empty_plot <- ggplot() + theme_void()
heatmap <-heatmap_all(topflop, asv_names, ds, tip_order, presentation =TRUE, size=10, adjustcor = FALSE)
phylo_tree_all_l <-ggdraw(phylo_tree_all[[2]]+ draw_plot(phylo_tree_all[[3]],2, 35, 0, 0))


combined_plot1 <- cowplot::plot_grid(phylo_tree_all_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.07)) # Paper
combined_plot1 <- cowplot::plot_grid(phylo_tree_all_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.015)) # Presentation

combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap[[1]] , nrow =2, ncol =1, rel_heights = c(0.002, 1))
combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                     rel_widths = c(1.1,1)) # paper
combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                     rel_widths = c(1.3,1)) # Presentation
print(combined_plot3)

#ggsave(filename = "figures/F_All_not_cutphylogenetical_tree_heat_map_presentation.png", height =6, width =8) # Presentation

#Heatmap physcio-chemical properties
combined_plot1 <- cowplot::plot_grid(phylo_tree_all_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.1))
combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap_all_pc_prop, nrow =2, ncol =1, rel_heights = c(0.002, 1))
combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                     rel_widths = c(1.1,0.9))

#ggsave(filename = "Figures/All_phylogenetical_tree_heat_map_prop.png", height =8, width =11)

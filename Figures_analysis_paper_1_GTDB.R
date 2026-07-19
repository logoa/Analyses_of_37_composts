# Graphics and analysis selection paper 1
# Author: Anja Logo
# Last analysis changes: 16.01.24
# Last code organiszation changes: 19.07.26

# Loading environemnet
source(file = "setup.R")
output_path <- "../../../04_Dissemination/01_Publikationen/figures_pdf/"

# Load data------
  ## 0.1 General------
  
  set.seed(100)  # set random seed

  # compost to exclude and factors
  compost.exclude = c("K6","K10","K17","K27","K37")
  compost.error = c("K27", "K37")
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
  df_BLW2 <- read.csv(file ="data/meta_data_compost_BLW_without_errors.csv", sep =";")
  factors <- c("batch", "site", "company", "comp.system", "site_ID") # Change to factor
  for (i in factors) {
    df_BLW2[, i] <- as.factor(df_BLW2[,i])
  }
  df_BLW2.r <- df_BLW2 %>% filter(!treatment %in% compost.exclude) # Reduce normal distribution to 37 composts
  rm(factors)
  
  # Are the continuous factors normally distributed?
  # Function to perform Shapiro-Wilk tests and create the results data frame
  perform_shapiro_test <- function(df, factor.list) {
    results <- data.frame(
      factor = factor.list,
      p.value = NA,
      normal = NA
    )
    for (i in 1:length(factor.list)) {
      x <- factor.list[i]
      sp.test <- shapiro.test(df[, x])
      results[i, "p.value"] <- round(sp.test$p.value, 3)
      results[i, "normal"] <- ifelse(sp.test$p.value > 0.05, 1, 0)
    }
    return(results)
  }
  results <- list()
  
  # Perform Shapiro-Wilk tests for different factor lists
  results[["factor.continous"]] <- perform_shapiro_test(df_BLW2.r, factor.continous)
  results[["factor.bioassay"]] <- perform_shapiro_test(df_BLW2.r, factor.bioassay)
  results[["factor.site"]] <- perform_shapiro_test(df_BLW2.r, factor.site)
  # Bind to one data.frame
  normal.dist.factors = do.call(rbind, results) 
  rownames(normal.dist.factors) <- NULL
  rm(results)
  
  # Peat meta data
  # Std.meta <-read.csv("Std_meta_ds_gp.csv", header = TRUE, sep = ";")
  
  # Compost_ID for changing the labeling of the data frame
  compost_ID <- read.csv(file ="data/Key_treatment_compostID.csv", sep= ";")
  # Change substrate naming
  compost_ID$compost_ID <- gsub("S", "NC", compost_ID$compost_ID)
  compost_ID$treatment = as.factor(compost_ID$treatment)
  compost_ID$compost_ID = as.factor(compost_ID$compost_ID)
  # Change substrate naming
  compost_ID$site_ID <- gsub("Peat", "no compost", compost_ID$site_ID)
  compost_ID$site_ID = as.factor(compost_ID$site_ID)

  compost_ID <- compost_ID %>% filter(!treatment %in% compost.error) # Exclude two error composts
  compost_ID$treatment <-droplevels(compost_ID$treatment)
  
  # Key for sequencing samples
  design.BLW2 <- read.csv(file ="data/design_BLW2.csv")
  design.BLW2$site <- as.factor(design.BLW2$site)
  
  # Add the new numbers to the design.BLW2
  design.BLW2 = merge(design.BLW2, compost_ID, by ="treatment", ) # Add new labeling
  rownames(design.BLW2) <- design.BLW2$ID
  design.BLW2$ID <- NULL
  design.BLW2 = design.BLW2[order(rownames(design.BLW2)),]
  design.BLW2 = design.BLW2 %>% filter(!treatment %in% compost.exclude) # 37 compost
  
  # For plots
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
  
  site.colors =  c("#81d4fa", "#01579b", "#ffecb3", "#e65100","#b39ddb","#ff9800", "#0288d1", "lightgrey")
  site.colors.dark =  c("#81d4fa","#01579b","#9D8604FF","#e65100","#b39ddb", "#ff9800", "#0288d1")  
  
  
  bg_theme <- theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.3, vjust = 0, size = 14, color ="black"),
      axis.text.y = element_text(size = 14, color ="black"),
      text = element_text(size = 14, color ="black"),
      axis.title = element_text(size=14),
      legend.text = element_text(size=14)
    )
  
  # Color for bioassays
  color.cress.dark <- rgb(143, 170, 220,maxColorValue = 255)
  color.cress.light <- rgb(189, 215, 238, maxColorValue =255)
  color.cugu.dark <- rgb(169, 209, 142, maxColorValue =255)
  color.cugu.light <- rgb(197, 224, 180, maxColorValue =255)
  color.curs.dark <- rgb(244, 177, 131, maxColorValue =255)
  color.curs.light <- rgb(248, 203, 173, maxColorValue =255)
  color.dark <- c(color.cress.dark, color.cugu.dark, color.curs.dark)
  color.light <-c(color.cress.light, color.cugu.light, color.curs.light)
  
  ## 0.2 Bioassay data------
  # Gu-cress
  df_cp <-read.csv("data/df_cp.csv", header =TRUE, sep=",")
  factors <- colnames(df_cp[,1:9])
  df_cp[factors] <-lapply(df_cp[factors], factor)
  df_cp$treatment <- as.character(df_cp$treatment)
  df_cp$treatment <- factor(df_cp$treatment, levels=unique(df_cp$treatment))
  df_cp$batch <- as.integer(df_cp$batch)
  
  # Gu-cucumber
  df_cu_gu <- read.csv("data/df_cu_gu.csv", header = TRUE, sep = ",")
  factors <- colnames(df_cu_gu[,c(2:7)])
  df_cu_gu[factors] <-lapply(df_cu_gu[factors], factor)
  df_cu_gu$treatment <- as.character(df_cu_gu$treatment)
  df_cu_gu$treatment <- factor(df_cu_gu$treatment, levels=unique(df_cu_gu$treatment)) 
  
  levels_to_relevel <- c("Std5", "Std4", "Std3", "Std2")
  for (level in levels_to_relevel) {
    df_cu_gu$treatment <- relevel(df_cu_gu$treatment, level)
    }
  df_cu_gu$batch <- as.integer(df_cu_gu$batch)

  # Rs-cucumber
  df_cu_rs <-read.csv("data/df_cu_rs.csv", header = TRUE, sep = ",")
  df_cu_rs$treatment <- as.factor(df_cu_rs$treatment)
  #df_cu_rs$treatment <- factor(df_cu_rs$treatment, levels=unique(df_cu_rs$treatment))
  #levels_to_relevel <- c("Std5", "Std4", "Std3", "Std2")
  #for (level in levels_to_relevel) {
    #df_cu_rs$treatment <- relevel(df_cu_rs$treatment, level)}
  
  #levels(df_cu_rs$treatment) =c("Std2","Std3","Std4", "Std5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15", 
    #                            "K17","K18","K19","K20","K21","K22","K23","K24","K25","K26","K28","K29","K30" ,  
    #                           "K31","K32","K33","K34","K35","K36","K39","K40","K41","K42","K43","K44","K45",
    #                           "K46","K47","K48","K49")  # Change order again
  
  rm(levels_to_relevel, factors)
  
  # Normal distribution for bioassays
  normal.dist <-read.csv(file = "data/Bioassay_normal_distribution_filled.csv", sep =",")
  normal.dist$X <-NULL
  # Change transformation to actual formula
  normal.dist <- normal.dist %>% mutate(
    formula = case_when(
      transformation == "none" ~ "y ~ x",
      transformation =="log+1" ~ "log(y+1) ~x",
      transformation == "sqrt" ~ "sqrt(y) ~x"))
  
  ## 0.3 Sequencing data ----

  # asv file without rarefying
  asv <- read.table("data/asv.B.BLW2.txt")
  ISS <- read.table("data/ISS.B.BLW2.txt")
  tax <- read.table("data/tax.BLW2.txt")
  asv.rob.comp.avg <- read.table("data/asv.rob.comp.avg.txt")
  
  # rarefied
  ISS.rob <- read.table("data/ISS.rob.txt")
  ISS.rob.avg <- read.table("data/ISS.rob.avg.txt")
  ISS.rob.comp <- read.table("data/ISS.rob.comp.txt")
  ISS.rob.comp.avg <- read.table("data/ISS.rob.comp.avg.txt")
  ISS.comp.avg <- read.table("data/ISS.comp.avg.txt") #Rarefied but not only the robust ones
  
  # Clean-up classification
  tax.class <- read.table(file ="data/tax.red.txt") # Only classified names
  tax.class.rob <- tax.class[colnames(ISS.rob.comp.avg),]
  tax.class.rownames <- rownames(tax.class.rob)
  tax.class.rob <- lapply(tax.class.rob, function(x) gsub(".*?__", "", x, perl = TRUE)) %>% as.data.frame()
  rownames(tax.class.rob) <- tax.class.rownames
  tax.class.rob <- tax.class.rob[colnames(ISS.rob.comp.avg),] # Only robust ones with unclassified marked
  rm(tax.class.rownames)
  
  # Alpha diversity for 37 compost + peat substrate
  
  # Meta file with averages
  df_BLW2.div <- read.csv("data/df_BLW2.div.csv")
  df_BLW2.div$X <-NULL
  df_BLW2.div$batch <- as.factor(df_BLW2.div$batch)
  
  # Replicates separately with peat
  ISS.alpha <- read.csv("data/ISS.alpha.B.BLW2.csv")
  rownames(ISS.alpha) <- ISS.alpha$X; ISS.alpha$X <- NULL
  
  # Beta diversity 37 composts (using 100 rarefied abundance table)
  ISS.bray <- read.table("data/ISS.bray.txt")
  ISS.bray.avg <- read.table("data/ISS.bray.avg.txt") %>% as.matrix() #BC-dis mean of four replicates 
  ISS.bray.avg.comp <-ISS.bray.avg[!grepl("Std", rownames(ISS.bray.avg)),!grepl("Std", colnames(ISS.bray.avg))]
  
  ISS.jac <- read.table("data/ISS.jac.txt")
  ISS.jac.avg <- read.table("data/ISS.jac.avg.txt")
  ISS.jac.avg.comp <- ISS.jac.avg[!grepl("Std", rownames(ISS.jac.avg)),!grepl("Std", colnames(ISS.jac.avg))]
  
  # ind ASVs BLW1
  
  ASVs.names.blast.BLW1.all <-read.table(file ="data/ASVS_names_blast_BLW1_all.txt")
  ASVs.names.blast.BLW1.best <- read.table(file ="data/ASVS_names_blast_BLW1_best.txt")
  ASVs.names.blast.BLW1.all.398 <- read.table(file ="data/ASVS_names_blast_BLW1_all_389.txt")
  ASVs.names.blast.BLW1.best.398 <- read.table(file ="data/ASVS_names_blast_BLW1_best_389.txt")
  
  # Design table for all samples
  
  design.BLW2 <- read.csv(file ="data/design_BLW2.csv")
  design.BLW2 <- design.BLW2 %>% filter(!treatment %in% compost.exclude)
  design.BLW2 <- design.BLW2 %>% filter(!treatment == "Std1")
  rownames(design.BLW2) <- design.BLW2$ID
  design.BLW2.ID <-left_join(design.BLW2, compost_ID, by = "treatment" )
  rownames(design.BLW2.ID) <- design.BLW2.ID$ID; design.BLW2.ID$ID <-NULL
  
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
  plot_box_ast <- function(data_sum_all, data.test, color.code, variable, y.text) {
    
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
  
  # Cor.test and scatterplot
  # Scatter plot

  plot_scatter <-function(data, variable1, variable2, method = "spearman", xaxis = "Relative biomass\n[%]", yaxis = "Rating\n[5-25]", method.text = "rho = "){
    data$x <-data[,variable1];data$y <- data[, variable2]
    cor_res  <- cor.test(data$x, data$y, method = method)
    rho <- round(cor_res$estimate, 2)
    p_val <- ifelse(cor_res$p.value < 0.001, "< 0.001", paste("=", round(cor_res$p.value,3)))
    plot <- ggplot(data = data, aes(x = x, y = y)) +
      geom_point() +
      geom_text(x = min(data$x), y = max(data$y), 
                label = paste0(method.text, rho, ", p ", p_val),
                hjust = 0, vjust = 1, size = 4) +
      # Label the axes with the variable names
      labs(x = xaxis, y = yaxis) +
      theme_bw()
    return(plot)
  }
  
  # Function: select batch concentration
  batch.select.fun <- function(data, batch.select){
    data %>%
      filter(paste(batch, conc) %in% paste(batch.select$batch, batch.select$conc))
  }
  
  sig.differences <- function(x){
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
  
  generate_distinct_colors <- function(n) {
    colors <- vector("character", n)  # Initialize an empty vector for colors
    
    # Generate n distinct colors
    for (i in 1:n) {
      # Generate a random color
      new_color <- rgb(runif(1), runif(1), runif(1) * 0.8)  # RGB with some constraints
      # Check for uniqueness
      while (any(colors == new_color)) {
        new_color <- rgb(runif(1), runif(1), runif(1) * 0.8)
      }
      colors[i] <- new_color
    }
    
    return(colors)
  }
  
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
      write.csv(PBC, file = paste0(path ,path.plant, ifelse(cutoff == TRUE, "_cut",""), "_PBC.csv"))
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

  
  
  # Remove all unnessary factors etc.
  rm(level, i, size, background_pr, mypalette)

# Sequencing overview----
  # rarecurve(t(asv), step = 100, col = "blue", cex = 0.6) # ASVs

  # Phyla
  tax.comp.avg <- tax.class[ISS.comp.avg %>% colnames(),]
  tax.comp.avg$Phyla %>% table() %>% sort()
  
  # Families
  tax.comp.avg$Family %>% table() %>% sort() %>% nrow()
  1/nrow(tax.comp.avg)* (nrow(tax.comp.avg)-4731)*100
  
  # Genus
  tax.comp.avg$Genus %>% table() %>% sort() %>% View()
  1/nrow(tax.comp.avg)* (nrow(tax.comp.avg)-9902)*100
  
  # Present in all 37 composts
  ISS.comp.avg.all <- ISS.comp.avg[, apply(ISS.comp.avg, 2, function(col) !any(col == 0))]
  ncol(ISS.comp.avg.all)/ ncol(ISS.comp.avg) *100
  tax.class[colnames(ISS.comp.avg.all),"Family"] %>% table() %>% sort() %>% nrow()
  tax.class[colnames(ISS.comp.avg.all),"Genus"] %>% table() %>% sort() %>% nrow()
  
  # Present in all except 10, 20 and 23
  ISS.comp.avg.all <- ISS.comp.avg[!rownames(ISS.comp.avg) %in% c("K19", "K30", "K33"),
                                   apply(ISS.comp.avg[!rownames(ISS.comp.avg) %in% c("K19", "K30", "K33"),], 2, function(col) !any(col == 0))]
  
  
  ncol(ISS.comp.avg.all)/ ncol(ISS.comp.avg) *100
  
  # Present in 9/10 composts
  
  threshold <- 33
  
  # Removing columns that don't have at least 90% non-zero values
  ISS.comp.avg.all <- ISS.comp.avg[, apply(ISS.comp.avg, 2, function(col) sum(col != 0) >= threshold)]
  ncol(ISS.comp.avg.all)/ ncol(ISS.comp.avg) *100
  tax.class[colnames(ISS.comp.avg.all),"Family"] %>% table() %>% sort()
  tax.class[colnames(ISS.comp.avg.all),"Genus"] %>% table() %>% sort()
  
  # Venndiagram preperation
  
  ISS.comp.avg$treatment = rownames(ISS.comp.avg)
  merged_data =merge(df_BLW2[,c("treatment", "batch", "site")],ISS.comp.avg , by = "treatment")
  rownames(merged_data) <- merged_data$treatment
  merged_data$treatment <- NULL
  ISS.comp.avg$treatment <-NULL

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
  rm(custom_colors, asv_tables_by_site, asv_tables_by_site_filtered, sets, venn.plot,
     merged_data, composting_sites)

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
  rm(custom_colors, asv_tables_by_batch, asv_tables_by_batch_filtered, sets, venn.plot,
     merged_data, batch)
  
  png("Figures/Sup4_Venn_site_batch.png", width = 4000, height = 2000 , res=300)
  
  # Arrange and save the plots
  grid.arrange(
    grobs = list(venn1, venn2),
    ncol = 2  # Number of columns
  )
  
  # Close the PNG device
  dev.off()

# 1 Bioassay results-----
  ## S2_Sup_rel_rating-----
  
  data1 <-df_cp %>% filter(case_when(
                batch == 1 ~ conc == 0.45,
                batch == 2 ~ conc == 0.45,
                batch == 3 ~ conc == 0.45, 
                batch == 4 ~ conc == 0.45
  ))
  p1 <-plot_scatter(data1, "rel", "rating", yaxis = "Rating\n[(-1)-3]")
  data2 <-df_cu_gu %>% filter(case_when(
    batch == 1 ~ conc == 0.45,
    batch == 2 ~ conc == 0.45,
    batch == 3 ~ conc == 1.35, 
    batch == 4 ~ conc == 1.35
  ))
  p2 <-plot_scatter(data2, "rel", "rating")
  data3 <-df_cu_rs %>% filter(case_when(
    batch == 1 ~ conc == 1.4,
    batch == 2 ~ conc == 0.8,
    batch == 3 ~ conc == 0.8, 
    batch == 4 ~ conc == 0.8
  ))
  p3<-plot_scatter(data3, "rel", "rating")
  
  ggarrange(p1, p2,p3, nrow= 1, ncol =3, labels = c("A", "B", "C"))
  #ggsave(filename = "Figures/Sup2_Bioassays_cor_rating_rel.png", height = 3.5, width =10)
  #ggsave(filename = paste0(output_path, "Sup_rel_rating.pdf"), height = 3.5, width =10)
  rm(data1, data2, data3, p1,p2,p3)
  
  ## Figure1_Disease_suppression-----
  batch.select <- data.frame(batch = c(1, 2, 3, 4), conc = c(rep(0.45,4)))
  data.nm = batch.select.fun(normal.dist, batch.select) %>% filter(plant.pathogen == "globi-cress",
                                                                   variable == "rel")
  # loop for statistics
  list = list()
  for (i in 1:4) {
    data = df_cp %>% filter(batch == i, conc == data.nm$conc[i], !treatment %in% compost.exclude) 
    TEST_summary(data, data$treatment, data$rel, i)
    data_sum$batch <- i
    list[[i]] <- data_sum%>% arrange(mean)
  }
  data_sum_all =rbind(list[[1]], list[[2]], list[[3]], list[[4]]) # This has to be done manually
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x)) # change factor levels
  
  # For 2D- plot later & significant differences indicator
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, 1, 0))
  sig.rel.cp = data_sum_all[data_sum_all$sig ==1,]$x

  data_sum_all$sig =ifelse(data_sum_all$sig ==1, "*", "")
  data_sum_all = data_sum_all %>% filter(! x %in% compost.exclude)
  data_sum_all = data_sum_all[order(data_sum_all$mean),]
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x))
  
  # Data selection for plot & ordered by mean
  data.test = batch.select.fun(df_cp %>% filter(!treatment %in% compost.exclude), batch.select)
  data.test$variable = data.test$rel
  data.test = data.test %>% filter(!treatment %in% compost.exclude)
  data.test$treatment <-droplevels(data.test$treatment)
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  data.test = inner_join(data.test, compost_ID[, c("treatment", "site_ID")], by ="treatment")
  data.test$treatment = as.factor(data.test$treatment)
  
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  
  # Adaptions for new labeling
  order_treatment = compost_ID[match(levels(data.test$treatment), compost_ID$treatment), "compost_ID"]
  
  y.text <- "Relative biomass \n [%]"
  
  # Option with boxplot
  data.test1 = data.test
  value_max1 = value_max
  data_sum_all1 = data_sum_all
  order_treatment1 = order_treatment
  
  g1 =data.test1 %>%
    ggplot(aes(x = treatment, y = variable)) +
    geom_boxplot(aes(fill = site_ID), position = position_dodge(width = 0.9), width= 0.8, outlier.size = 0.8,
                 outlier.shape = NA, coef =0) +
    labs ( x = "", y = y.text)+
    graphic.style +
    theme(strip.placement = "outside", strip.text.x = element_text(size=16))+
    scale_fill_manual(name="Composting Site", values = site.colors)+
    geom_jitter(width=0.05,alpha = 0.4, shape =16, size=1)+
    scale_x_discrete(labels= order_treatment1)+
    geom_text(data = value_max1, aes(x=treatment, y = 10 + quantile, label = data_sum_all1$sig), vjust=0, size= 7, colour ="darkred")+
    scale_y_continuous(limits = c(0, 140), breaks = seq(0, 140, by = 20))
  
  #ggsave(g1,filename="Figures_all_conc/Publication1/Gu_cress_rel_boxplot_points_peat_colored.png", height =6, width = 14)
  
  ### 1B
  batch.select <- data.frame(batch = c(1, 2, 3, 4), conc = c(0.45, 0.45, 1.35, 1.35))
  data.nm = batch.select.fun(normal.dist, batch.select) %>% filter(plant.pathogen == "globi-cuc",
                                                                   variable == "rel")
  # loop for statistics
  list = list()
  for (i in 1:4) {
    data = df_cu_gu %>% filter(batch == i, conc == data.nm$conc[i], !treatment %in% compost.exclude) 
    TEST_summary(data, data$treatment, data$rel, i)
    data_sum$batch <- i
    list[[i]] <- data_sum%>% arrange(mean)
  }
  
  data =list[[2]]
  list[[2]] = data[c(4,1:3,5:nrow(data)),]
  data =list[[4]]
  list[[4]] = data[c(2,1,3:nrow(data)),]
  data_sum_all =rbind(list[[1]], list[[2]], list[[3]], list[[4]]) # This has to be done manually
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x)) # change factor levels

  # For 2D- plot later & significant differences indicator
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, 1, 0))
  sig.rel.cu.gu = data_sum_all[data_sum_all$sig ==1,]$x
  
  data_sum_all$sig =ifelse(data_sum_all$sig ==1, "*", "")
  data_sum_all = data_sum_all %>% filter(! x %in% compost.exclude)
  data_sum_all = data_sum_all[order(data_sum_all$mean),]
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x))
  
  # Data selection for plot & ordered by mean
  data.test = batch.select.fun(df_cu_gu %>% filter(!treatment %in% compost.exclude), batch.select)
  data.test$variable = data.test$rel
  data.test = data.test %>% filter(!treatment %in% compost.exclude)
  data.test$treatment <-droplevels(data.test$treatment)
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  data.test = inner_join(data.test, compost_ID[, c("treatment", "site_ID")], by ="treatment")
  
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  
  # Adaptions for new labeling
  order_treatment = compost_ID[match(levels(data.test$treatment), compost_ID$treatment), "compost_ID"]
  
  y.text <- "Relative biomass \n [%]"
  data.test2 = data.test
  value_max2 = value_max
  data_sum_all2 = data_sum_all
  order_treatment2 = order_treatment
  
  # Option with boxplot
  g2 =data.test2 %>%
    ggplot(aes(x = treatment, y = variable)) +
    geom_boxplot(aes(fill = site_ID), position = position_dodge(width = 0.9), width= 0.8, outlier.size = 0.8,
                 outlier.shape = NA, coef =0) +
    labs ( x = "", y = y.text)+
    graphic.style +
    theme(strip.placement = "outside", strip.text.x = element_text(size=16))+
    scale_fill_manual(name="Composting\nsite", values = site.colors)+
    geom_jitter(width=0.05,alpha = 0.4, shape =16, size=1)+
    scale_x_discrete(labels= order_treatment2)+
    geom_text(data = value_max2, aes(x=treatment, y = 10 + quantile, label = data_sum_all2$sig), vjust=0, size= 7, colour ="darkred")+
    scale_y_continuous(limits = c(-1, 140), breaks = seq(0, 140, by = 20))
  #ggsave(g2, filename="Figures_all_conc/Publication1/Gu_cuc_rel_boxplot_points_peat_colored.png", height =6, width = 14)
  
  ### 1C
  batch.select <- data.frame(batch = c(1, 2, 3, 4), conc = c(1.4, 0.8, 0.8, 0.8))
  data.nm = batch.select.fun(normal.dist, batch.select) %>% filter(plant.pathogen == "rsolani-cuc",
                                                                   variable == "rel")
  # loop for statistics
  list = list()
  for (i in 1:4) {
    data = df_cu_rs %>% filter(batch == i, conc == data.nm$conc[i], !treatment %in% compost.exclude) 
    TEST_summary(data, data$treatment, data$rel, i)
    data_sum$batch <- i
    list[[i]] <- data_sum%>% arrange(mean)
  }
  
  data =list[[1]]
  list[[1]] = data[c(5,1:4,6:nrow(data)),]
  data =list[[2]]
  list[[2]] = data[c(2,1,3:nrow(data)),]
  data_sum_all =rbind(list[[1]], list[[2]], list[[3]], list[[4]]) # This has to be done manually
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x)) # change factor levels
  
  # For 2D- plot later & significant differences indicator
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, 1, 0))
  sig.rel.cu.rs = data_sum_all[data_sum_all$sig ==1,]$x
  
  data_sum_all$sig =ifelse(data_sum_all$sig ==1, "*", "")
  data_sum_all = data_sum_all %>% filter(! x %in% compost.exclude)
  data_sum_all = data_sum_all[order(data_sum_all$mean),]
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x))
  
  # Data selection for plot & ordered by mean
  data.test = batch.select.fun(df_cu_rs %>% filter(!treatment %in% compost.exclude), batch.select)
  data.test$variable = data.test$rel
  data.test = data.test %>% filter(!treatment %in% compost.exclude)
  data.test$treatment <-droplevels(data.test$treatment)
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  data.test = inner_join(data.test, compost_ID[, c("treatment", "site_ID")], by ="treatment")
  
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  
  # Adaptions for new labeling
  order_treatment = compost_ID[match(levels(data.test$treatment), compost_ID$treatment), "compost_ID"]
  
  y.text <- "Relative biomass \n [%]"
  data.test3 = data.test
  value_max3 = value_max
  data_sum_all3 = data_sum_all
  order_treatment3 = order_treatment
  
  # Option with boxplot
  g3 =data.test3 %>%
    ggplot(aes(x = treatment, y = variable)) +
    geom_boxplot(aes(fill = site_ID), position = position_dodge(width = 0.9), width= 0.8, outlier.size = 0.8,
                 outlier.shape = NA, coef =0) +
    labs ( x = "", y = y.text)+
    graphic.style +
    theme(strip.placement = "outside", strip.text.x = element_text(size=16))+
    scale_fill_manual(name="Composting\nsite", values = site.colors)+
    geom_jitter(width=0.05,alpha = 0.4, shape =16, size=1)+
    scale_x_discrete(labels= order_treatment3)+
    geom_text(data = value_max3, aes(x=treatment, y = 10 + quantile, label = data_sum_all3$sig), vjust=0, size= 7, colour ="darkred")+
    scale_y_continuous(limits = c(-1, 140), breaks = seq(0, 140, by = 20))
  
  #ggsave(g3, filename="Figures_all_conc/Publication1/Rs_cuc_rel_boxplot_points_peat_colored.png", height =6, width = 14)
  
  ggarrange(g1, g2, g3, common.legend = TRUE, legend = "bottom", labels =c("A", "B", "C"), nrow=3)
  #ggsave(filename="Figures/Figure_1.png", height =14, width = 15)
  #ggsave(filename= paste0(output_path,"/Figure1_Disease_suppression.pfd"), height =14, width = 15)
  
  # Venndiagramm 1D
  #x <- list(cress = sig.rel.cp, cugu = sig.rel.cu.gu, curs =sig.rel.cu.rs)
  #ggVennDiagram(x)
  
  rm(g1, g2, g3, list, data.test1, data.test2, data.test3, batch.select, data_sum_all,
     data_sum_all1, data_sum_all2, data_sum_all3, y.text, data_sum, data, value_max,
     value_max1, value_max2, value_max3, i,
     order_treatment, order_treatment1, order_treatment2, order_treatment3, factors, x)
  
  ## S4_Sup_cor_pathplant ----

  # Correlation among pathogen plants systems
  # Data are roughly normal distributed
  #hist(df_BLW2.r$rel.cp)
  #hist(df_BLW2.r$rel.cu.gu)
  #hist(df_BLW2.r$rel.cu.rs)
  
  g1 <-plot_scatter(df_BLW2.r, "rel.cp", "rel.cu.gu",
               xaxis= "Cress-Gu \n Relative biomass [%]",
               yaxis ="Cucumber-Gu \n Relative biomass [%]",
               method ="pearson", method.text = "r = ")+
    geom_smooth(method ="lm", color="darkred", se=FALSE)
  
  g2 <-plot_scatter(df_BLW2.r, "rel.cp", "rel.cu.rs",
                    xaxis= "Cress-Gu \n Relative biomass [%]",
                    yaxis ="Cucumber-Rs \n Relative biomass [%]",
                    method ="pearson", method.text = "r = ")
  
  g3 <-plot_scatter(df_BLW2.r, "rel.cu.gu", "rel.cu.rs",
                    xaxis= "Cucumber-Gu \n Relative biomass [%]",
                    yaxis ="Cucumber-Rs \n Relative biomass [%]",
                    method ="pearson", method.text = "r = ")
  
  ggarrange(g1, g2, g3, labels =c("A", "B", "C"), nrow=1, ncol=3)
  
  #ggsave(filename = "Figures/Cor_plot_pearson_relative_biomass.png", height =4, width =12)
  #ggsave(filename = paste0(output_path,"Sup_cor_pathplant.pdf"), height =4, width =12)
  rm(g1, g2, g3)
  ## S5_Sup_dsbatchsite------
  df_BLW2_filtered <- df_BLW2.r %>% filter(!site %in% c("Spreitenbach", "Frick"))
  df_BLW2_filtered$site_ID <- droplevels(df_BLW2_filtered$site_ID)
  margin.plot = theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  ylim = ylim(0,110)
  color.batch = c("grey80", "grey60", "grey40", "grey25")
  
  anova_plot <- function(data, x, y, color) {
    
    data$x <- x
    data$x <-as.factor(data$x)
    data$y <- y
  
    # Statistics
    ANOVA= aov(y ~ x, data= data)
    print(summary(ANOVA))
    TUKEY = TukeyHSD(ANOVA, conf.level = 0.95)
    letters =multcompLetters4(ANOVA, TUKEY)
    summarized = data %>% group_by(x) %>% summarize(y= as.numeric(quantile( y, prob=c(.75))))
    par(mfrow =c(2,2))
    print(plot(ANOVA))
    
    order = levels(data$x)
    
    # Plot
    g = ggplot(data, aes(x= x, y= y, fill =x))+
      geom_text(data = summarized, aes(x = x, y= y, label =letters$x$Letters[order]), vjust =-1, hjust = -1.5, color ="darkred")+
      geom_boxplot()+ theme(legend.position = "none")+ theme_classic()+
      scale_fill_manual(values =color)+ theme(legend.position = "none")+
      ylab("Relative biomass\n[%]") + xlab("")+ theme(axis.text.x = element_text(angle = 0, hjust=0.3, vjust=0))+
      theme(text=element_text(size=14))+ margin.plot+ylim
    
    return(g)
    
  }
  
  g1 =anova_plot(df_BLW2.r, df_BLW2.r$batch_ID, df_BLW2.r$rel.cp, color.batch)
  g2 =anova_plot(df_BLW2_filtered, df_BLW2_filtered$site_ID, df_BLW2_filtered$rel.cp, site.colors)
  g3 =anova_plot(df_BLW2.r, df_BLW2.r$batch_ID, df_BLW2.r$rel.cu.gu, color.batch)
  g4 =anova_plot(df_BLW2_filtered, df_BLW2_filtered$site_ID, df_BLW2_filtered$rel.cu.gu, site.colors)
  g5 =anova_plot(df_BLW2.r, df_BLW2.r$batch_ID, df_BLW2.r$rel.cu.rs, color.batch)
  g6= anova_plot(df_BLW2_filtered, df_BLW2_filtered$site_ID, df_BLW2_filtered$rel.cu.rs, site.colors)
  
  ggarrange(g1, g2, g3, g4, g5, g6, nrow=3, ncol=2,
            labels = c("A", "B", "C", "D", "E", "F"),
            font.label = list(size =14, face= "plain"))
  
  #ggsave(filename = "Figures/Sup_disease_suppression_batch_compostorigin.png", height = 10, width =8)
  #ggsave(filename = paste0(output_path,"Sup_dsbatchsite.pdf"), height = 10, width =8)
  
  ## S1_Sup_growth_promotion----
  ### A
  batch.select <- data.frame(batch = c(1, 2, 3, 4), conc = c(rep(0,4)))
  data.nm = batch.select.fun(normal.dist, batch.select) %>% filter(plant.pathogen == "globi-cress",
                                                                   variable == "gp")
  # loop for statistics
  list = list()
  for (i in 1:4) {
    data = df_cp %>% filter(batch == i, conc == data.nm$conc[i], !treatment %in% compost.error) 
    TEST_summary(data, data$treatment, data$gp, i)
    data_sum$batch <- i
    list[[i]] <- data_sum%>% arrange(mean)
  }
  data_sum_all =rbind(list[[1]], list[[2]], list[[3]], list[[4]]) # This has to be done manually
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x)) # change factor levels
  
  # For 2D- plot later & significant differences indicator
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, 1, 0))
  sig.gp.cp = data_sum_all[data_sum_all$sig ==1,]$x
  
  data_sum_all$sig =ifelse(data_sum_all$sig ==1, "*", "")
  data_sum_all = data_sum_all %>% filter(! x %in% compost.error)
  data_sum_all = data_sum_all[order(data_sum_all$mean),]
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x))
  
  # Data selection for plot & ordered by mean
  data.test = batch.select.fun(df_cp %>% filter(!treatment %in% compost.error), batch.select)
  data.test$variable = data.test$gp
  data.test$treatment <-droplevels(data.test$treatment)
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  
  # Here is the problem because K6, K10, K17 are excluded, maybe call for example Ex1, Ex2, Ex3??
  data.test = inner_join(data.test, compost_ID[, c("treatment", "site_ID")], by ="treatment")
  
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  
  # Adaptions for new labeling
  order_treatment = compost_ID[match(levels(data.test$treatment), compost_ID$treatment), "compost_ID"]
  
  y.text <- "Growth promotion \n [%]"
  
  # Option with boxplot
  STD = c("Std2", "Std3", "Std4", "Std5")
  data.test1 = data.test %>% filter(!treatment %in% STD)
  value_max1 = value_max %>% filter(!treatment %in% STD)
  data_sum_all1 = data_sum_all %>% filter(!x %in% STD)
  order_treatment1 =  order_treatment[!grepl("^NC", order_treatment)]
  
  g1 =data.test1 %>%
    ggplot(aes(x = treatment, y = variable)) +
    geom_boxplot(aes(fill = site_ID), position = position_dodge(width = 0.9), width= 0.8, outlier.size = 0.8,
                 outlier.shape = NA, coef =0) +
    labs ( x = "", y = y.text)+
    graphic.style +
    theme(strip.placement = "outside", strip.text.x = element_text(size=16))+
    scale_fill_manual(name="Composting \n site", values = site.colors)+
    geom_jitter(width=0.05,alpha = 0.4, shape =16, size=1)+
    scale_x_discrete(labels= order_treatment1)+
    geom_text(data = value_max1, aes(x=treatment, y = 10 + quantile, label = data_sum_all1$sig), vjust=0, size= 7, colour ="darkred")+
    scale_y_continuous(limits = c(-80, 180), breaks = seq(-80, 180, by = 40))
  
  #ggsave(g1,filename="Figures_all_conc/Publication1/Gu_cress_rel_boxplot_points_peat_colored.png", height =6, width = 14)
  
  ### 1B
  batch.select <- data.frame(batch = c(1, 2, 3, 4), conc = c(rep(0,4)))
  data.nm = batch.select.fun(normal.dist, batch.select) %>% filter(plant.pathogen == "globi-cuc",
                                                                   variable == "gp")
  # loop for statistics
  list = list()
  for (i in 1:4) {
    data = df_cu_gu %>% filter(batch == i, conc == data.nm$conc[i], !treatment %in% compost.error) 
    TEST_summary(data, data$treatment, data$gp, i)
    data_sum$batch <- i
    list[[i]] <- data_sum%>% arrange(mean)
  }
  
  data_sum_all =rbind(list[[1]], list[[2]], list[[3]], list[[4]]) # This has to be done manually
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x)) # change factor levels
  
  # For 2D- plot later & significant differences indicator
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, 1, 0))
  sig.gp.cu.gu = data_sum_all[data_sum_all$sig ==1,]$x
  
  data_sum_all$sig =ifelse(data_sum_all$sig ==1, "*", "")
  data_sum_all = data_sum_all %>% filter(! x %in% compost.error)
  data_sum_all = data_sum_all[order(data_sum_all$mean),]
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x))
  
  # Data selection for plot & ordered by mean
  data.test = batch.select.fun(df_cu_gu %>% filter(!treatment %in% compost.error), batch.select)
  data.test$variable = data.test$gp
  data.test = data.test %>% filter(!treatment %in% compost.error)
  data.test$treatment <-droplevels(data.test$treatment)
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  data.test = inner_join(data.test, compost_ID[, c("treatment", "site_ID")], by ="treatment")
  
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  
  # Adaptions for new labeling
  order_treatment = compost_ID[match(levels(data.test$treatment), compost_ID$treatment), "compost_ID"]
  
  y.text <- "Growth promotion \n [%]"
  STD = c("Std2", "Std3", "Std4", "Std5")
  data.test2 = data.test %>% filter(!treatment %in% STD)
  value_max2 = value_max %>% filter(!treatment %in% STD)
  data_sum_all2 = data_sum_all %>% filter(!x %in% STD)
  order_treatment2 =  order_treatment[!grepl("^NC", order_treatment)]
  
  
  # Option with boxplot
  g2 =data.test2 %>%
    ggplot(aes(x = treatment, y = variable)) +
    geom_boxplot(aes(fill = site_ID), position = position_dodge(width = 0.9), width= 0.8, outlier.size = 0.8,
                 outlier.shape = NA, coef =0) +
    labs ( x = "", y = y.text)+
    graphic.style +
    theme(strip.placement = "outside", strip.text.x = element_text(size=16))+
    scale_fill_manual(name="Composting\n site", values = site.colors)+
    geom_jitter(width=0.05,alpha = 0.4, shape =16, size=1)+
    scale_x_discrete(labels= order_treatment2)+
    geom_text(data = value_max2, aes(x=treatment, y = 10 + quantile, label = data_sum_all2$sig), vjust=0, size= 7, colour ="darkred")+
    scale_y_continuous(limits = c(-80, 180), breaks = seq(-80, 180, by = 40))
  #ggsave(g2, filename="Figures_all_conc/Publication1/Gu_cuc_gu_boxplot_points_peat_colored.png", height =6, width = 14)
  
  ### 1C
  batch.select <- data.frame(batch = c(1, 2, 3, 4), conc = rep(0,4))
  data.nm = batch.select.fun(normal.dist, batch.select) %>% filter(plant.pathogen == "rsolani-cuc",
                                                                   variable == "gp")
  # loop for statistics
  list = list()
  for (i in 1:4) {
    data = df_cu_rs %>% filter(batch == i, conc == data.nm$conc[i], !treatment %in% compost.error) 
    TEST_summary(data, data$treatment, data$gp, i)
    data_sum$batch <- i
    list[[i]] <- data_sum%>% arrange(mean)
  }
  data_sum_all =rbind(list[[1]], list[[2]], list[[3]], list[[4]]) # This has to be done manually
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x)) # change factor levels
  
  # For 2D- plot later & significant differences indicator
  STD = data_sum_all[grepl("S", data_sum_all$x),]$letters %>% sapply( transform_characters)
  data_sum_all <- data_sum_all %>%
    mutate(sig = case_when(
      batch == "1" ~ !grepl(STD[1], data_sum_all$letters),
      batch == "2" ~ !grepl(STD[2], data_sum_all$letters),
      batch == "3" ~ !grepl(STD[3], data_sum_all$letters),
      batch == "4" ~ !grepl(STD[4], data_sum_all$letters)
    ))  %>% mutate(sig = ifelse(sig, 1, 0))
  sig.gp.cu.rs = data_sum_all[data_sum_all$sig ==1,]$x
  
  data_sum_all$sig =ifelse(data_sum_all$sig ==1, "*", "")
  data_sum_all = data_sum_all %>% filter(! x %in% compost.error)
  data_sum_all = data_sum_all[order(data_sum_all$mean),]
  data_sum_all$x <- as.character(data_sum_all$x) %>% factor(levels=unique(data_sum_all$x))
  
  # Data selection for plot & ordered by mean
  data.test = batch.select.fun(df_cu_rs %>% filter(!treatment %in% compost.error), batch.select)
  data.test$variable = data.test$gp
  data.test = data.test %>% filter(!treatment %in% compost.error)
  data.test$treatment <-droplevels(data.test$treatment)
  data.test$treatment <- factor(data.test$treatment, levels = levels(data_sum_all$x))
  data.test = inner_join(data.test, compost_ID[, c("treatment", "site_ID")], by ="treatment")
  
  value_max = data.test %>% group_by(treatment) %>% dplyr::summarize(max_value = max(variable),
                                                                     quantile = quantile(variable, probs = 0.75),
                                                                     batch = mean(batch))
  value_max$batch <- as.factor(value_max$batch)
  
  # Adaptions for new labeling
  order_treatment = compost_ID[match(levels(data.test$treatment), compost_ID$treatment), "compost_ID"]
  
  y.text <- "Growth promotion \n [%]"
  STD = c("Std2", "Std3", "Std4", "Std5")
  data.test3 = data.test %>% filter(!treatment %in% STD)
  value_max3 = value_max %>% filter(!treatment %in% STD)
  data_sum_all3 = data_sum_all %>% filter(!x %in% STD)
  order_treatment3 =  order_treatment[!grepl("^NC", order_treatment)]
  
  # Option with boxplot
  g3 =data.test3 %>%
    ggplot(aes(x = treatment, y = variable)) +
    geom_boxplot(aes(fill = site_ID), position = position_dodge(width = 0.9), width= 0.8, outlier.size = 0.8,
                 outlier.shape = NA, coef =0) +
    labs ( x = "", y = y.text)+
    graphic.style +
    theme(strip.placement = "outside", strip.text.x = element_text(size=16))+
    scale_fill_manual(name="Composting\nsite", values = site.colors)+
    geom_jitter(width=0.05,alpha = 0.4, shape =16, size=1)+
    scale_x_discrete(labels= order_treatment3)+
    geom_text(data = value_max3, aes(x=treatment, y = 10 + quantile, label = data_sum_all3$sig), vjust=0, size= 7, colour ="darkred")+
    scale_y_continuous(limits = c(-80, 180), breaks = seq(-80, 180, by = 40))
  
  #ggsave(g3, filename="Figures_all_conc/Publication1/Rs_cuc_gp_boxplot_points_peat_colored.png", height =6, width = 14)
  
  ggarrange(g1, g2, g3, common.legend = TRUE, legend = "bottom", labels =c("A", "B", "C"), nrow=3)
  #ggsave(filename="Figures_all_conc/Publication1/Sup_Growthpromotion.png", height =18, width = 14)
  #ggsave(filename=paste0(output_path, "Sup_growth_promotion.pdf"), height =12, width = 13)
  
  rm(g1, g2, g3, list, data.test1, data.test2, data.test3, batch.select, data_sum_all,
     data_sum_all1, data_sum_all2, data_sum_all3, y.text, data_sum, data, value_max,
     value_max1, value_max2, value_max3, i, level,
     order_treatment, order_treatment1, order_treatment2, order_treatment3, factors,data.test)
  
  ## S3_Sup_pathcon-----
  
  # select peat substrate from all three data frames and create a new
  factors = c("treatment", "path.plant", "conc", "rel", "biomass", "batch")
  
  # Create new data.frame to faciliate the plotting
  df_cu_gu$path.plant = paste0(df_cu_gu$pathogen, "-", df_cu_gu$plant)
  df_cu_rs$path.plant = paste0(df_cu_rs$pathogen, "-", df_cu_rs$plant)
  cress = df_cp[grepl("S", df_cp$treatment), factors]
  cress$conc <- as.character(cress$conc)
  gucu = df_cu_gu[grepl("S", df_cu_gu$treatment),factors]
  gucu$conc <- as.character(gucu$conc)
  rscu = df_cu_rs[grepl("S", df_cu_rs$treatment), factors]
  rscu$conc <- as.character(rscu$conc)
  df_std = rbind(cress, gucu, rscu)
  df_std = df_std %>% filter(!conc %in% c("0", "0.15"))
  df_std[df_std$conc %in% c("0.6","0.45"),"conc"] <- "low"
  df_std[df_std$conc %in% c("1.4","1.35"),"conc"] <- "high"
  df_std[df_std$conc == 0.8 & df_std$batch %in% c(1,2,3), "conc"] <- "low"
  df_std[df_std$conc == 0.8,"conc"] <- "high"
  df_std$path.plant <- as.factor(df_std$path.plant)
  df_std$conc <- factor(df_std$conc, levels = c("low", "high"))
  
  df_std_select =df_std %>% 
    mutate(select_conc = case_when(
      path.plant == "globi-cress" & conc =="low" ~ "Selection",
      path.plant == "globi-cucumber" & batch %in% c(1,2) & conc == "low" ~ "Selection",
      path.plant == "globi-cucumber" & batch %in% c(3,4) & conc =="high" ~ "Selection",
      path.plant == "rsolani-cucumber" & batch %in% c(1,4)& conc == "high" ~ "Selection",
      path.plant == "rsolani-cucumber" & batch %in% c(2,3)& conc =="low"~ "Selection",
      TRUE ~ "n" # Default case
    )) %>% filter(select_conc == "Selection")
  
  # Statistics
  data = df_std %>% filter(path.plant =="globi-cress" & conc =="low") 
  data = df_std %>% filter(path.plant =="globi-cress" & conc =="high") 
  data = df_std %>% filter(path.plant =="globi-cucumber" & conc =="low") 
  data = df_std %>% filter(path.plant =="globi-cucumber" & conc =="high") 
  data = df_std %>% filter(path.plant =="rsolani-cucumber" & conc =="low") 
  data = df_std %>% filter(path.plant =="rsolani-cucumber" & conc =="high") 
  
  data = df_std_select %>% filter(path.plant =="globi-cucumber")
  data = df_std_select %>% filter(path.plant =="rsolani-cucumber")
  
  ANOVA= aov(rel ~ treatment, data= data)
  summary(ANOVA)
  TUKEY = TukeyHSD(ANOVA, conf.level = 0.95)
  multcompLetters4(ANOVA, TUKEY)
  par(mfrow =c(2,2))
  plot(ANOVA)
  par(mfrow=c(1,1))
  
  KRUSKAL =kruskal.test(rel ~treatment, data =data)
  pairwise_results <- pairwise.wilcox.test(data$rel, data$treatment, p.adjust.method = "BH")
  KRUSKAL =tri.to.squ(pairwise_results$p.value)
  multcompLetters(KRUSKAL)
  
  df_std <- df_std %>% mutate(path.plant = case_when(
    path.plant == "globi-cress" ~ "cress-Gu",
    path.plant == "globi-cucumber" ~ "cucumber-Gu",
    path.plant =="rsolani-cucumber" ~"cucumber-Rs"
  ))
  
  
  df_labels <- df_std %>%
    group_by(treatment, conc, path.plant) %>%
    summarize(
      max_y = max(rel, na.rm = TRUE),
      label = NA  # Example label: number of observations
    ) %>%
    ungroup()
  
  df_labels <- df_labels[order(df_labels$path.plant,df_labels$conc), ]
  df_labels$label <- c("a", "a", "a", "a",
                       "a", "b", "ab", "a",
                       "b", "b", "a", "ab",
                       "bc", "c", "a", "ab",
                       "a", "b", "b", "b",
                       "a", "b", "b", "ab"
                       )
  
  g1 =ggplot(data = df_std, aes(x = treatment, y = rel, fill = conc))+
    geom_boxplot()+ facet_grid(conc ~ path.plant)+ xlab("")+
    ylab("Relative biomass\n [%]")+
    graphic.style+
    scale_fill_manual(values = c("lightgrey", "gray42"))+
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.text = element_text(size = 14, color = "white"),
          strip.background = element_rect(fill = "black" ),
          axis.text.x = element_blank())+
    ylim(0,120) +
    geom_text(
      data = df_labels,
      aes(x = treatment, y = max_y + 5, label = label), # Adjust y for spacing above boxplot
      inherit.aes = FALSE, # Prevent overwriting aesthetics
      size = 4, vjust = 0
    )

  df_labels2 <- df_std_select %>%
    group_by(treatment, conc, path.plant) %>%
    summarize(
      max_y = max(rel, na.rm = TRUE),
      label = NA  # Example label: number of observations
    ) %>%
    ungroup()
  df_labels2 <- df_labels2[order(df_labels2$path.plant,df_labels2$conc), ]
  df_labels2$label <- c(rep("a",12))
  
  g2 =ggplot(data = df_std_select, aes(x = treatment, y = rel, fill = conc))+
    geom_boxplot()+ facet_grid(select_conc~ path.plant)+
    xlab("")+
    ylab("Relative biomass\n [%]")+
    graphic.style+
    scale_fill_manual(values = c("lightgrey", "gray42"))+
    scale_x_discrete(labels= c("NC-I", "NC-II","NC-III", "NC-IV"))+
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.text = element_text(size = 14, color = "white"),
          strip.background = element_rect(fill = "black"),
          strip.text.x = element_blank())+
    ylim(0,120)+
    geom_text(
      data = df_labels2,
      aes(x = treatment, y = max_y + 5, label = label), # Adjust y for spacing above boxplot
      inherit.aes = FALSE, # Prevent overwriting aesthetics
      size = 4, vjust = 0
    )
  
  ggarrange(g1, g2, nrow =2, heights = c(1.9,1.1))
  #ggsave("Bacteria_BLW2_23/BLW2_Compost_only_paper1/Figures/Sup_Peat_concentrations.png", height = 12, width =12)
  #ggsave(paste0(output_path, "Sup_pathcon.pdf"), height = 12, width =12)
  
  
  rm(factors, cress, gucu, rscu, df_std, g1, g2, df_std_select, KRUSKAL, TUKEY,
     ANOVA, data, data.test, pairwise_results, STD, df_labels, df_labels2)


# 2 Compost properties ----
  ## S6_Sup_PCA_pc-----
  factor.select <- factor.continous[!factor.continous %in% c("NO3.Nmin", "Nmin", "Corg.N", "Corg", "max_WHC", "basal", "Ntot")]
  comp.meta.pca <- prcomp(df_BLW2.r[,factor.select], center = TRUE, scale. = TRUE) 
  
  loadings <- comp.meta.pca$rotation %>% as.data.frame() # Loadings, indicator for importance of a factors
  scores <- as.data.frame(comp.meta.pca$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  data <-cbind(df_BLW2.r, scores)
  
  # Extract percentages
  var <- comp.meta.pca$sdev^2 / sum(comp.meta.pca$sdev^2)
  perc <- var * 100
  #color = c("orange", "firebrick","orchid4","blue3","cornflowerblue")
  #color = c( "#F1F1F1","#FFC59E" ,"#E1BB4E", "yellow", "#9BB306", "#26A63A", "darkgreen")
  color = c("white", "black")
  list = list()
  for (i in c("rel.cp", "rel.cu.gu", "rel.cu.rs", "age")) {
    data$x <- data[,i]
    list[[i]] =ggplot(data, aes(x = PC1, y = PC2, color =x)) +
      xlim(-6.2, 2.9) + ylim(-4, 4.3)+
      geom_point(aes(fill = x), size=4, pch=21, colour ="darkgrey") + # Add points
      bg_theme+
      labs(x = paste0("PC1 (",round(perc[1],2), "%)"), y = paste0("PC2 (",round(perc[2],2), "%)")) +
      coord_fixed(1)+
      theme(legend.position = "right", legend.title = element_text(size =12))+
      scale_fill_gradientn(colours = color,  name= "Relative\nbiomass [%]")
     #scale_fill_gradientn(colours = c("yellow", "darkgreen"), name= "Compost age [%]")
    
  }
  g2 = list[[1]]
  g3 = list[[2]]
  g4 = list[[3]]
  
  #ggsave(ggarrange(g2, g3, g4, nrow=3, labels =c("A", "B", "C")),
  #       file= paste0(output_path,"Sup_PCA_pc.pdf"), height =9, width =6)
  
  ## Figure2_PCOA_compost_properties-----
  
  scaling_factor <- min(
    (max(scores$PC1) - min(scores$PC1) / (max(loadings$PC1) - min(loadings$PC1))),
    (max(scores$PC2) - min(scores$PC2) / (max(loadings$PC2) - min(loadings$PC2)))
  )
  loadings <- loadings*scaling_factor
  loadings$var <- rownames(loadings)
  color <- c( "#E5F6DA", "#264612")
  
  g5 <-ggplot(data, aes(x = PC1, y = PC2, color =age)) +
    geom_point(data = data[!is.na(data$age),], aes(fill = x), size=4, pch=21, colour ="black") + 
    scale_fill_gradientn(colours = color, name= "Compost age\n[days]")+
    geom_point(data = data[is.na(data$age), ],
               size = 4, pch=21,  fill = "white", colour ="gray73") +
    bg_theme+
    xlim(-6.2, 2.9) + ylim(-4, 4.3)+
    labs(x = paste0("PC1 (",round(perc[1],2), "%)"), y = paste0("PC2 (",round(perc[2],2), "%)")) +
    coord_fixed(1)+
    theme(legend.position = "bottom")+
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "darkred", size=0.8) +
    geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = var), 
                    color = "black", size = 4)
  
  #ggsave(g5, file ="Figures/PCOA_compost_proterties_age.png", height = 7, width =7)
  
  # PERMANOVA
  dist <-vegdist(scale(data[,factor.select]), method = "euclidean") 
  adonis_result <- adonis(dist ~ rel.cp, data)
  adonis_result <- adonis(dist ~ rel.cu.gu, data)
  adonis_result <- adonis(dist ~ rel.cu.rs, data)
  adonis_result <- adonis(dist ~ batch, data)
  adonis_result <- adonis(dist ~ site, data)
 
  # Compost age without K34 and K44
  data1 = data %>% filter(!treatment %in% c("K43", "K44"))
  dist <-vegdist(scale(data1[,factor.select]), method = "euclidean") 
  adonis_result <- adonis(dist ~ age, data1)

  adonis_result$aov.tab 
  
  # PCOA & plot
  factor.select <- factor.continous[!factor.continous %in% c("NO3.Nmin", "Nmin", "Corg.N", "Ctot", "max_WHC", "basal", "Ntot")]
  comp.meta.pca <- prcomp(df_BLW2.r[,factor.select], center = TRUE, scale. = TRUE) 
  loadings <- comp.meta.pca$rotation # Loadings, indicator for importance of a factors
  scores <- as.data.frame(comp.meta.pca$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  data =cbind(df_BLW2.r, scores)
  
  # Extract percentages
  var <- comp.meta.pca$sdev^2 / sum(comp.meta.pca$sdev^2)
  perc <- var * 100
  
  # Draw polygon
  find_hull <- function(data) data[chull(data$PC1, data$PC2), ]
  hulls <- plyr::ddply(data, "site_ID", find_hull)
  
  g1 =ggplot(data, aes(x = PC1, y = PC2)) +
    xlim(-6.2, 2.9) + ylim(-4, 4.3)+
    bg_theme+
    geom_point(aes(color = site_ID), size=2) + # Add points
    labs(x = paste0("PC1 (",round(perc[1],2), "%)"), y = paste0("PC2 (",round(perc[2],2), "%)")) +
    scale_fill_manual(values = site.colors.dark, name = "Composting\nsite")+
    scale_color_manual(values = site.colors.dark, name = "Composting\nsite")+
    theme(legend.position = "bottom")+
    geom_text_repel(aes(label = Number, color = site_ID), size = 4, show.legend = FALSE)+
    coord_fixed(ratio =1)+
    geom_polygon(data=hulls,alpha=0.2,aes(fill=site_ID))

  
  # PERMANOVA
  data1 = data %>% filter(!site_ID %in% c("F", "G"))
  dist <-vegdist(scale(data1[,factor.select]), method = "euclidean") 
  adonis_result <- adonis(dist ~ site_ID, data1)
  adonis_result$aov.tab 
  
  #ggsave(ggarrange(g1, g5, ncol=2, labels = c("A", "B")), file=paste0(output_path,"Figure2_PCOA_compost_properties.pdf"), height = 5, width = 9)
  
  rm(g1, g2, g3, g4, legend1, legend2, final_plot, data1, dist, adonis_result, plots_without_legends,
     comp.meta.pca, data, hulls, list, legends, loadings, scores, color, dist, i, perc, var, factor.select,find_hull)
  
  ## S7_Sup_coramongcp------
  data <-df_BLW2.r
  factor.list = c("age", "DS","max_WHC","pH","sal","OD550","NO2","NO3","NH4","Nmin","NO3.Nmin","PO4","Ntot","Corg","Corg.N", "basal", "FDA")
  correlation_estimates <- data.frame(Variable1 = character(0), Variable2 = character(0),
                                      Estimate = numeric(0),
                                      P_Value = numeric(0),
                                      P_adj = numeric(0),
                                      n = numeric(0))
  correlation_matrix <- matrix(NA, nrow = length(factor.list), ncol = length(factor.list))
  
  for (i in 1:length(factor.list)){
    for (j in 1:length(factor.list)) {
      var1 <- factor.list[i]
      var2 <- factor.list[j]
      
      cor <-cor.test(data[[var1]], data[[var2]], method = "spearman", use= "p")
      n <- data[c(var1, var2)] %>% na.omit() %>% nrow()
      estimate <- cor$estimate
      p_value <- cor$p.value
      correlation_estimates <- rbind(correlation_estimates, data.frame(Variable1 = var1, Variable2 = var2,
                                                                       Estimate = estimate,
                                                                       P_value = p_value, n = n))
      }
  }
  
  # adjust for multiple testing
  correlation_estimates$padj <-correlation_estimates$P_value %>% p.adjust(method ="BH")
  rownames(correlation_estimates) = NULL
  correlation_estimates[,3] = round(correlation_estimates[,3], 2)
  correlation_estimates[,4] = round(correlation_estimates[,4], 3)
  correlation_estimates[,6] = round(correlation_estimates[,6], 3)
  correlation_estimates <-correlation_estimates %>% arrange(Variable1)
  
  #write.csv(correlation_estimates, file ="Output/cor_comp_cp_pvalue_corrected.csv", row.names = FALSE)
  
  #Heat map (done by ChatGPT)
  
  corr_res <- Hmisc::rcorr(as.matrix(data %>% select(factor.list)), type = "spearman")
  cor_matrix <- corr_res$r
  p_matrix <- corr_res$P
  
  p_adj_matrix <- matrix(p.adjust(as.vector(p_matrix), method = "BH"), nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  rownames(p_adj_matrix) <- rownames(p_matrix)
  colnames(p_adj_matrix) <- colnames(p_matrix)
  
  cor_df <- reshape2::melt(cor_matrix, varnames = c("Var1", "Var2"), value.name = "Correlation")
  p_df <- reshape2::melt(p_adj_matrix, varnames = c("Var1", "Var2"), value.name = "p.adj")
  cor_df <- merge(cor_df, p_df, by = c("Var1", "Var2"))
  
  mask_insignificant <- function(cor_matrix, p_adj_matrix, alpha = 0.05) {
    cor_matrix[p_adj_matrix > alpha] <- NA
    return(cor_matrix)
  }
  cor_matrix_masked <- mask_insignificant(cor_matrix, p_adj_matrix)
  #pdf(file = paste0(output_path, "Sup_coramongcp.pdf"), width = 9, height = 9)
  corrplot(cor_matrix_masked, 
           method = "color", 
           type = "upper", 
           tl.col = "black", 
           tl.srt = 45, 
           na.label = "NA", 
           addCoef.col = "black", 
           p.mat = p_adj_matrix, 
           sig.level = 0.05, 
           insig = "blank",
           diag =FALSE)
  
  #dev.off()

## Tbl2_Disease_suppression overview_tbl-----

  # Composting properties Pearson and Spearman mixed, without multiple testing correction
  data <-df_BLW2.r
  
  factor.list = c(factor.continous, "age")
  correlation_estimates <- data.frame(Variable1 = character(0), Variable2 = character(0),
                                      Estimate = numeric(0), Method = character(0),
                                      P_Value = numeric(0), n = numeric(0))
  correlation_matrix <- matrix(NA, nrow = length(factor.list), ncol = length(factor.bioassay))
  
  
  for (i in 1:length(factor.list)){
    for (j in 1:length(factor.bioassay)) {
      var1 <- factor.list[i]
      var2 <- factor.bioassay[j]
      is_var1_normal <- normal.dist.factors[normal.dist.factors$factor == factor.list[i], 3]
      is_var2_normal <- normal.dist.factors[normal.dist.factors$factor == factor.bioassay[j], 3]
      
      # Determine the correlation method based on normality
      if (is_var1_normal == 1 && is_var2_normal == 1) {
        corr_method <- "pearson"
      } else {
        corr_method <- "spearman"
      }
      cor <-cor.test(data[[var1]], data[[var2]], method = corr_method, use= "p")
      n <- data[c(var1, var2)] %>% na.omit() %>% nrow()
      estimate <- cor$estimate
      p_value <- cor$p.value
      correlation_estimates <- rbind(correlation_estimates, data.frame(Variable1 = var1, Variable2 = var2,
                                                                       Method = corr_method, Estimate = estimate,
                                                                       P_value = p_value, n = n))
      
      if (p_value < 0.05) {
        # Calculate the correlation and store it in the appropriate position
        correlation <- cor(data[[var1]], data[[var2]], method = corr_method, use= "p")
        correlation_matrix[i, j] <- correlation
      }
    }
  }
  rownames(correlation_matrix) = factor.list; colnames(correlation_matrix) = factor.bioassay
  correlation_matrix[is.na(correlation_matrix)] <-0
  corrplot(correlation_matrix, method = "number", diag= TRUE)
  
  rownames(correlation_estimates) = NULL
  correlation_estimates[,4:5] = round(correlation_estimates[,4:5], 3)
  
  #write.csv(correlation_estimates, file= "Output/bioassay_compost_prop_cor_selection2", row.names =FALSE)

  # Composting properties only Spearman, with multiple testing correction
  data <-df_BLW2.div
  #data <- df_BLW2.r %>% filter(batch != 2)
  factor.list = c("age", factor.continous, "mean.sobs", "mean.even", "mean.shannon", "mean.ivsimp")
  factor.bioassay.red <- c("rel.cp", "rel.cu.gu", "rel.cu.rs")
  correlation_estimates <- data.frame(Variable1 = character(0), Variable2 = character(0),
                                      Estimate = numeric(0),
                                      P_Value = numeric(0),
                                      P_adj = numeric(0),
                                      n = numeric(0))
  correlation_matrix <- matrix(NA, nrow = length(factor.list), ncol = length(factor.bioassay.red))
  
  for (i in 1:length(factor.list)){
    for (j in 1:length(factor.bioassay.red)) {
      var1 <- factor.list[i]
      var2 <- factor.bioassay.red[j]
      
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
  rownames(correlation_matrix) = factor.list; colnames(correlation_matrix) = factor.bioassay.red
  correlation_matrix[is.na(correlation_matrix)] <-0
  corrplot(correlation_matrix, method = "number", diag= TRUE)
  
  rownames(correlation_estimates) = NULL
  correlation_estimates[,3] = round(correlation_estimates[,3], 2)
  correlation_estimates[,4] = round(correlation_estimates[,4], 3)
  correlation_estimates <-correlation_estimates %>% arrange(Variable2)
  
  # write.csv(correlation_estimates, file= "Output/bioassay_compost_prop_cor_selection2_spearman.csv", row.names =FALSE)
  #write.csv(correlation_estimates, file= "Output/bioassay_compost_prop_cor_selection2_spearman_without_batchII.csv", row.names =FALSE)
  
  # Adjusting all the p.values together
  p.values <-read.csv(file ="Output/Indicator_disease_suppression_multiple_testing.csv", sep =";")
  
  list <-list()
  for (i in c("pc", "alpha", "beta")) {
  list[[i]] <-p.values %>% filter(topic == i)  %>%
      mutate(rel.cp.adj = p.adjust(rel.cp.con, "BH") %>% round(3),
             rel.cu.gu.adj = p.adjust(rel.cu.gu.con, "BH") %>% round(3),
             rel.cu.rs.adj = p.adjust(rel.cu.rs.con, "BH") %>% round(3))  
  }
  p.values.adj <-do.call(rbind, list)
  rownames(p.values.adj) <- NULL
 # write.csv(p.values.adj, file ="Output/Indicator_disease_suppression_multiple_testing_filled.csv")
  
  
## Tbl Compost set & site and compost properties
  data <- df_BLW2.r %>% select(batch, site, "age", all_of(factor.continous)); j <-1
  data <- data %>% filter(! site %in% c("Spreitenbach", "Frick")); j <- 2
  
  colums <- c("variable", "Chi", "p.value", "p.adj")
  results <- matrix(ncol = length(colums), nrow =(ncol(data)-4)) %>% as.data.frame()
  colnames(results) <- colums; rm(colums)
  for(i in 3:20){
    data$x <- data[,j]
    data$y <- data[,i]
    results[i-2, 1] <- colnames(data)[i]
    kruskal_results = kruskal.test(y~x, data =data[!is.na(data$x),])
    results[i -2, 2] <-kruskal_results$statistic %>% round(2)
    results[i-2, 3] <- kruskal_results$p.value %>% round(3)}
  
  results$p.adj <-p.adjust(results$p.value, method = "BH") %>% round(3)
  
  #write.csv(results, file= "Output/Batch_effect_compost_prop.csv")
 # write.csv(results, file= "Output/Compostingsite_effect_compost_prop.csv")
  
## S8-S10_Sup_corcpcress /cugu /curs-----
  
  # Assuming `data` is your dataset with 17 parameters and 3 response variables
  # Replace "data" with your actual data frame name
  # params: Columns 1 to 17
  # responses: Columns 18 to 20
  params <- c("age", "DS", "max_WHC", "pH",  "sal", "OD550","NO2","NO3","NH4","Nmin",
              "NO3.Nmin", "PO4", "Ntot","Corg", "Corg.N", "basal", "FDA") 
  responses <- c("rel.cp", "rel.cu.gu", "rel.cu.rs")
  
  # Function to create scatter plots for a response variable
  create_plots <- function(data, response_var, params) {
    plots <- list()  # To store individual plots
    
    for (param in params) {
      # Remove rows where the current parameter or response variable is missing
      clean_data <- data %>% filter(!is.na(.data[[param]]), !is.na(.data[[response_var]]))
      
      # Create scatter plot using ggplot2
      p <- ggplot(clean_data, aes_string(x = param, y = response_var)) +
        geom_point() +
        ylab("Relative biomass [%]")+
        bg_theme +
        theme_minimal()
      
      # Append plot to list
      plots[[param]] <- p
    }
    
    return(plots)
  }
  
  # Generate plots for each response variable
  for (response in responses) {
    # Create plots for this specific response variable
    plots <- create_plots(df_BLW2.r, response, params)
    
    # Arrange all 17 plots in a grid (4x5 grid in this case)
    grid_arrangement <- do.call(grid.arrange, c(plots, ncol = 4, nrow = 5))
    
    # Save the plot grid to a file or display it
    ggsave(paste0(output_path, "plot", response, ".pdf"), grid_arrangement, width = 10, height = 13)
  }
  
# Figure4A (only bacteria) Alpha diversity-----
  
  # Good's coverage

  Goods <- rep(NA, ncol(asv))
  for (i in 1:ncol(asv)) {
    Goods[i] = 1 - (sum(asv[,i] == 1) / sum(asv[,i]))
  }
  
  cbind(colnames(asv), Goods) %>% View()
  
  
  summary(Goods)
  sd(Goods)*100

  data = design.BLW2.ID %>% filter(!treatment %in% compost.exclude)
  data1 = ISS.alpha[rownames(ISS.alpha) %in% rownames(data),]
  rownames(data) == rownames(data1)
  
  data1[data1 %>% rownames() %>% startsWith("Std"),] %>% pull(sobs) %>% mean()
  
  data1[data1 %>% rownames() %>% startsWith("K"),] %>% select(sobs, evenness, shannon, invsimpson) %>% summary()

  select <- data1[data1 %>% rownames() %>% startsWith("K"),] %>% select(sobs, evenness, shannon, invsimpson) 
  sd(select$sobs); sd(select$evenness); sd(select$shannon); sd(select$invsimpson)
  
# Figure 3
  g1 =ggplot(data1 , aes(x = reorder(data$compost_ID, sobs, na.rm =T), y = sobs, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
    labs(x = "", y ="")+
    scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
    scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
    scale_y_continuous(limits = c(800, 5300), breaks = seq(1000, 5000, by = 500))+
    bg_theme +
    theme(legend.position = "right", legend.key = element_blank(), legend.title = element_text(size =14),
          legend.text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(color = "lightgrey"),
          axis.title = element_text(size=14))

  
  #ggsave(file ="Figures/Alpha_diversity_colored_by_site.png", height =6, width = 14)
  ggsave(g1, file ="Figures/Alpha_diversity_colored_by_site.png", height =5, width = 14)
  
  ## S14_alpha_diversity (only bacteria)------
  
  s1 <-ggplot(data1 , aes(x = reorder(data$compost_ID, evenness, na.rm =T), y = evenness, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
    labs(x = "", y ="Evenness")+
    scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
    scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
    scale_y_continuous(limits = c(0.6, 0.9), breaks = seq(0.6, 0.9, by = 0.05))+
    bg_theme +
    theme(legend.position = "none", legend.key = element_blank(), legend.title = element_text(size =14),
          legend.text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(color = "lightgrey"))
  
  s2 <- ggplot(data1 , aes(x = reorder(data$compost_ID, shannon, na.rm =T), y = shannon, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
    labs(x = "", y ="Shannon diversity")+
    scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
    scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
    scale_y_continuous(limits = c(4.5, 7.5), breaks = seq(4.5, 7.5, by = 0.5))+
    bg_theme +
    theme(legend.position = "none", legend.key = element_blank(), legend.title = element_text(size =14),
          legend.text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(color = "lightgrey"))
  
  s3 <- ggplot(data1 , aes(x = reorder(data$compost_ID, invsimpson, na.rm =T), y = invsimpson, color =data$site_ID, fill = data$site_ID))+ geom_boxplot()+
    labs(x = "", y ="Inversed Simpson diversity")+
    scale_color_manual(values= c(site.colors.dark, "grey22"), name= "Composting\nsite")+
    scale_fill_manual(values= adjustcolor(c(site.colors.dark,"grey22"), alpha.f = 0.2),  name= "Composting\nsite")+
    scale_y_continuous(limits = c(10, 480), breaks = seq(10, 480, by = 50))+
    bg_theme +
    theme(legend.position = "none", legend.key = element_blank(), legend.title = element_text(size =14),
          legend.text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(color = "lightgrey"))
  
  legend <- get_legend(s1+ theme(legend.position = "bottom"))
  plot_grid(s1, s2, s3, legend, nrow =4, axix ="none", rel_heights = c(1,1,1,0.15),
            ncol =1,labels =c("A", "B", "C")) 
  
  plot_grid(s1, s2, legend, nrow =3, axix ="none", rel_heights = c(1,1,0.15),
            ncol =1,labels =c("A", "B")) 
  
  #ggsave(filename ="Figures/Sup4_Alpha_diversity.png", height =12, width =13)
  ggsave(filename ="Figures/Sup4_Alpha_diversity.png", height =10, width =14)
  
  #Statistics for Supplementary file
  
  #sink(paste0("Output/Composting_site_batch_effec_alpha_div.txt"))
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
  sink()
  
  #paste("Correlation with disease suppression")
  
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
    #ggsave(filename = paste0("Figures/", "F_BLW2_cor_ds_",alpha.div[j], ".png"), height =4, width =12, units= "in")
  }
  
  g5 <- plot_scatter(df_BLW2.div, variable1 = "rel.cp", variable2 = "mean.sobs", xaxis = "Disease suppression GU-cress [%] ", yaxis = "Bacterial richness" )+
    geom_text_repel(aes(label = Number, show.legend =F,  max.overlaps =20))+
    geom_smooth(se =F, method ="lm", colour = "#64769A" )+ bg_theme
  

  color.dark.dark <- c("#64769A", "#769263", "#AA7B5B")
  
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
  

  #annotate(geom = "text", label = paste("rho =",round(cor1$estimate,2),"p =", round(cor1$p.value,3)), x =4500, y =95 ) +
    #annotate(geom = "text", label = paste("rho =",round(cor2$estimate,2),"p =", round(cor2$p.value,3)), x =4500, y =60 ) +
    #annotate(geom = "text", label = paste("rho =",round(cor3$estimate,2),"p =", round(cor3$p.value,3)), x =4500, y =35 ) +
   
  #ggsave(g2, filename = "Figures/Alpha_diversity_correlation_disease_suppression.png", height =5, width = 6)
  #ggsave(g2, filename = "Figures/Alpha_diversity_correlation_disease_suppression_light.png", height =5, width = 6)
  

  # Composting site
  
  # Statistics
  kruskal_result <- kruskal.test(mean.sobs ~ site_ID, data = df_BLW2.div %>% filter(!site_ID %in% c("F", "G")))
  
  # Plot
  g3 = ggplot(df_BLW2.div %>% filter(!site_ID %in% c("F", "G")), aes(x= site_ID, y= mean.sobs, fill = site_ID, color =site_ID))+
    geom_boxplot()+ theme(legend.position = "none")+
    scale_color_manual(values= c(site.colors.dark, "grey22"))+
    scale_fill_manual(values= adjustcolor(c(site.colors.dark[1:5],"grey22"), alpha.f = 0.2))+
    ylab("") + xlab("")+
    bg_theme +theme(axis.text.x = element_text(angle = 0, hjust=0.3, vjust=0))
  
  # Figure 3D
  # Statistics
  kruskal_result <- kruskal.test(mean.sobs ~ batch, data = df_BLW2.div)
  pairwise_results = pairwise.wilcox.test(df_BLW2.div$mean.sob, df_BLW2.div$batch, p.adj ="BH")
  # comparison_letters <- multcompLetters(as.vector(pairwise_results$p.value))
  KRUSKAL =tri.to.squ(pairwise_results$p.value)
  compact_letters_df =multcompLetters(KRUSKAL)
  
  
  value_max = df_BLW2.div %>% group_by(batch) %>% dplyr::summarize(max_value = max(mean.sobs),
                                                                     quantile = quantile(mean.sobs, probs = 0.75))
  value_max$batch <- as.factor(value_max$batch)
  value_max$letter <- c("a", "b", "ab", "a")
  
  # Plot wiht letters
  g4 = ggplot(df_BLW2.div , aes(x= batch, y= mean.sobs, fill = batch))+
    geom_boxplot()+
    scale_fill_manual(values =c (rep("white",4)))+
    scale_x_discrete(labels= c("I", "II", "III", "IV"))+
    ylab("Bacterial richness") + xlab("Compost set")+
    bg_theme + theme(axis.text.x = element_text(angle = 0, hjust=0.3, vjust=0))+
    geom_text(data =value_max, aes(x=batch, y = 100 + quantile, label = letter), hjust=-0.2, vjust =-0.2,  size= 4, colour ="black")
    
  g = plot_grid(g2, g3, g4, nrow =1, axix ="none", ncol =3, rel_widths = c(2,2,2),labels =c("B", "C", "D")) 
  plot_grid(g1, g, nrow=2, ncol =1, rel_heights = c(1.1,1), labels =c("A", ""))
  
  #ggsave(filename ="Figures/Alpha_diversity_plot_A_D.png", height =9, width =14)
  
  g = plot_grid(g5, g4, nrow =1, axis ="non", ncol =2, rel_widths = c(1,0.4), labels = c("C", "D"))
  ggsave(filename ="Figures/Alpha_diversity_plot_C_D.png", height =4.5, width =7)

  # Correlations batch separatly
  for (i in 1:4) {
    data <- df_BLW2.div %>% filter(batch ==i)
    cor <-cor.test(data$mean.sobs, data$rel.cp, method = "spearman")
    print(cor)
  }

  rm(g1, g2,g3, g4, data, data1, g, cor1, cor2, cor3, kruskal_result, pairwise_results, KRUSKAL, compact_letters,
     compact_letters_df, color1, color2, color3)
## Figure 3 Alpha diversity with fungi
  
  
# Tbl_S7_CP_BC------
  # Correlation analysis compost properties + compost age and alpha diversity metrics
  
  data <-df_BLW2.div
  data <-df_BLW2.div %>% filter(batch !=2)
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
  #write.csv(correlation_estimates, file ="Output/Alpha_diversity_compost_prop_multiple_testing_filled.csv", row.names = FALSE)
  #write.csv(correlation_estimates, file ="Output/Alpha_diversity_compost_prop_multiple_testing_filled_without_batchII.csv", row.names = FALSE)
  
  
  # Beta diversity
  
  data <- df_BLW2.r
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
  
  #write.csv(correlation_estimates, file ="Output/Beta_diversity_compost_prop_multiple_testing_filled.csv", row.names = FALSE)
  
  

# Figure3 Beta diversity ASV level (only bacteria)-----

  nmds.BLW2 <- metaMDS(ISS.bray.avg.comp)
  nmds.BLW2.points = nmds.BLW2$points %>% as.data.frame()
  colnames(nmds.BLW2.points) = c("nmds1", "nmds2")
  nmds.BLW2.points$treatment = rownames(nmds.BLW2.points)
  df_BLW2.beta <- merge(df_BLW2.r, nmds.BLW2.points , by ="treatment")
  nmds.BLW2$stress # 0.10
 
  # composting site
  find_hull <- function(df_BLW2.beta) df_BLW2.beta[chull(df_BLW2.beta$nmds1, df_BLW2.beta$nmds2), ]
  hulls <- plyr::ddply(df_BLW2.beta, "site_ID", find_hull)
  
  g1 =ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2))+
    geom_point(aes(color = site_ID))+
    theme_classic(base_size = 14) + theme(legend.position = "none") +
    scale_fill_manual(values = site.colors.dark, name = "Composting\nsite")+
    scale_color_manual(values = site.colors.dark, name = "Composting\nsite")+
    geom_text(x=0.6,y=-0.2,label='stress = 0.1', size=4)+
    geom_polygon(data=hulls,alpha=0.2,aes(fill=site_ID))+
    ylab("NMDS2")+xlab("NMDS1")+
    coord_fixed(ratio =1)+
    geom_text_repel(data = df_BLW2.beta, aes(x = nmds1 , y = nmds2 , label = Number, color = site_ID), show.legend = FALSE, max.overlaps =20)
  #ggsave(filename = "Figures/N;Ds_Betadiv_composting_site.png", height = 8, width =12)
  
  
  ## Disease suppression
  color = c("orange", "firebrick","orchid4","blue3","cornflowerblue")
  color = c( "#F1F1F1","#FFC59E" ,"#E1BB4E", "yellow", "#9BB306", "#26A63A", "darkgreen")
  
    g2 = ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cp))+
      theme_classic(base_size = 14) +
      geom_point(aes(fill = rel.cp), size =4, pch=21, colour ="darkgrey")+
      coord_fixed(ratio =1) +
      theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
    scale_fill_gradientn(colours = color,
                           limits = c(0,110), name= "Disease\n suppression [%]")
    
    g3= ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cu.gu))+
      theme_classic(base_size = 14) +
      geom_point(aes(fill = rel.cu.gu), size =4, pch=21, colour="darkgrey")+
      coord_fixed(ratio =1) +
      theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
      scale_fill_gradientn(colours = color,
                             limits =c(0,110), name= "Disease\nsuppression [%]")

    g4= ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cu.rs))+
      theme_classic(base_size = 14) +
      geom_point(aes(fill = rel.cu.rs), size =4, pch =21, colour="darkgrey")+
      coord_fixed(ratio =1) + 
      theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
      scale_fill_gradientn(colours = color,
                             limits = c(0,110),name= "Disease\nsuppression [%]")
    
  legend1 = ggpubr::get_legend(g2+ theme(legend.position = "bottom"))
  legend2 = ggpubr::get_legend(g1+ theme(legend.position = "bottom"))

  plots_without_legends <- plot_grid(g2, g3, g4, g1, ncol = 2, align = 'hv', labels = c("A", "B", "C", "D"))
  legends <- plot_grid(legend1, legend2, ncol = 2, align ='v')
  final_plot <- plot_grid(plots_without_legends, legends, ncol = 1, align = 'v', rel_heights = c(1, 0.1))
  # ggsave(final_plot, height =10, width =14, file ="Figures/NMDS_Betadiv_FigureA-D.png")
  # ggsave(final_plot, height =10, width =14, file ="Figures/NMDS_Betadiv_FigureA-D_light.png")
  
  # Plot only the three disease assays
  ggarrange(g2, g3, g4, common.legend=T, nrow =1, ncol =3, labels  =c ("G. ultimum -cress", "G. ultimum -cucumber", "R. solani - cucumber"))
  #ggsave(file ="Bacteria_BLW2_23/BLW2_Compost_only_paper1/Figures/NMDS_Betadiv_FigureA-C.png", height = 4, width =14)
  #ggsave(file ="Figures/NMDS_Betadiv_FigureA-C_light.png", height = 4, width =14)
  
  # PERMANOVA
  
  #sink(file ="Output/B_ANOVA_composting_site_batch.txt")
  paste("Composting site")
  data <- ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K48", "K49")),
                                 !(colnames(ISS.bray.avg.comp) %in% c("K48", "K49"))]
  
  adonis2(data ~ site, data = df_BLW2.r[!df_BLW2.r$site %in% c("Frick", "Spreitenbach"),], permutations = 999)
  
  paste("batch")
  adonis2(ISS.bray.avg.comp ~ batch, data = df_BLW2.r, permutations = 999)
  
  paste("Compost age")
  data = ISS.bray.avg.comp[!(rownames(ISS.bray.avg.comp) %in% c("K43", "K44")),
                           !(colnames(ISS.bray.avg.comp) %in% c("K43", "K44"))]
  adonis2(data  ~ age, data = df_BLW2.r %>% filter(!treatment %in% c("K43", "K44")), permutations = 999)
  
  paste("Disease suppression")
  adonis2(ISS.bray.avg.comp ~ rel.cp, data = df_BLW2.r, permutations = 999)
  adonis2(ISS.bray.avg.comp ~ rel.cu.gu, data = df_BLW2.r, permutations = 999)
  adonis2(ISS.bray.avg.comp ~ rel.cu.rs, data = df_BLW2.r, permutations = 999)
  #sink()
  
  rm(data, final_plot, g1, g2, g3, g4, hulls, legend1, legend2, legends, nmds.BLW2.points, 
     nmds.BLW2, plots_without_legends, color)

  ## S12_Sup_B_tax/S13_sup_F_tax-------
  
  # Preperation
  tax.class.red <- lapply(tax.class, function(x) sub(".*__", "", x)) %>% as.data.frame()
  rownames(tax.class.red) <- rownames(tax.class)
  ISS.comp.avg = as.matrix(ISS.comp.avg)
  ISS.prop = prop.table(ISS.comp.avg, margin = 1) * 100  # get proportions
  ISS.prop = as.data.frame(ISS.prop)
  ISS.comp.avg = as.data.frame(ISS.comp.avg)
  
  # Aggregate Phylum level
  ISS.p <- aggregate(t(ISS.prop) ~ Phyla, data = tax.class.red[colnames(ISS.prop),], FUN = sum)  # aggregate by phylum
  ISS.p <-ISS.p %>% relocate(Phyla, K7, K8, K9)
  ISS.p <- as.matrix(data.frame(ISS.p[, -1], row.names = ISS.p[, 1], check.names = FALSE))  # move first col to row.names
  ISS.p <- ISS.p[order(rowSums(ISS.p), decreasing = TRUE), ]  # order by rowSums
  ISS.p.nc <- as.data.frame(ISS.p[rownames(ISS.p) != "unclassified", ])  #remove unclassified group
  ISS.p.nc["others", ] <- colSums(ISS.p.nc[min(nrow(ISS.p.nc), 13):nrow(ISS.p.nc),
  ])  # merge rare groups to others
  ISS.p.nc["unclassified", ] <- ISS.p[rownames(ISS.p) == "unclassified", drop = F]  # add unclassified at the end
  ISS.p.plot <- ISS.p.nc[-(min(nrow(ISS.p.nc), 13):(nrow(ISS.p.nc) - 2)), ]  # remove rare groups
  ISS.p.plot <- as.matrix(ISS.p.plot)
  
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
  
  custom_x_labels <- compost_ID %>% dplyr::filter(!treatment %in% c("K6", "K10", "K17", "K27", "K37", "Std2",
                                                                  "Std3", "Std4", "Std5" )) %>% 
    dplyr::select(compost_ID) 
  
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
  mean_percentages <- rowMeans(ISS.g.plot)
  
  
  # Plot 1: Phyla
  
  # Calculate mean percentages for each phylum
  mean_percentages <- rowMeans(ISS.p.plot)
  # Create labels with mean percentages
  legend_labels <- paste(rownames(ISS.p.plot), sprintf("(%.1f%%)", mean_percentages))
  

  #pdf(file =paste0(output_path, "Sup_B_tax_A.pdf"), width=10, height=5)
  par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2))
  color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
             "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#999999", "#000000")
  barplot(ISS.p.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.8, axisnames = TRUE, names.arg = custom_x_labels$compost_ID)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.1, 0.2), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = legend_labels, pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  #dev.off()


  # Plot 2: Family level
  
  # Calculate mean percentages for each phylum
  mean_percentages <- rowMeans(ISS.f.plot)
  # Create labels with mean percentages
  legend_labels <- paste(rownames(ISS.f.plot), sprintf("(%.1f%%)", mean_percentages))
  
  #pdf(file =paste0(output_path, "Sup_B_tax_B.pdf"), width=10, height=5)
  par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2))
  color <-  c(generate_distinct_colors(20), "#999999", "#000000")
  barplot(ISS.f.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.8, axisnames = TRUE, names.arg = custom_x_labels$compost_ID)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.1, 0.2), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = legend_labels, pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  
  #dev.off()
  
  # Plot 3: Genera level
  
  # Calculate mean percentages for each phylum
  mean_percentages <- rowMeans(ISS.g.plot)
  # Create labels with mean percentages
  legend_labels <- paste(rownames(ISS.g.plot), sprintf("(%.1f%%)", mean_percentages))
  
  #pdf(file =paste0(output_path, "Sup_B_tax_C.pdf"), width=10, height=5)
  par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2))
  color <-  c(generate_distinct_colors(20), "#999999", "#000000")
  barplot(ISS.g.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.8, axisnames = TRUE, names.arg = custom_x_labels$compost_ID)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.1, 0.2), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = legend_labels, pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  
  #dev.off()
  

  ## S11_Beta (only bacteria) diversity all samples----
  
  nmds.BLW2 <- metaMDS(ISS.bray)
  nmds.BLW2.points = nmds.BLW2$points %>% as.data.frame()
  colnames(nmds.BLW2.points) = c("nmds1", "nmds2")
  df_BLW2.beta <- merge(design.BLW2.ID, nmds.BLW2.points , by =0)
  rownames(df_BLW2.beta) <- df_BLW2.beta$Row.names; df_BLW2.beta$Row.names <-NULL
  nmds.BLW2$stress # 0.06
  
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
    geom_text(x=0.75,y=-0.15, label='stress = 0.06', size=4)+
    ylab("NMDS2")+xlab("NMDS1")+
    geom_polygon(data=hulls,alpha=0.2,aes(fill=compost_ID))+
    coord_fixed(ratio =1)
  
  #ggsave(filename = "Figures/Sup4_Beta_diverstiy_all_samples.png", height = 7, width =12)
  
  adonis2(ISS.bray[!grepl("Std", rownames(ISS.bray)), 
                   !grepl("Std", colnames(ISS.bray))] ~ treatment, 
          data = df_BLW2.beta %>% filter(!treatment %in% c("Std2", "Std3", "Std4", "Std5")),
          permutations = 999)
  
# Figure5A&B Indicator Species analysis------
  
  ## Summarize on taxonomic level----
  # Match with previous project asv_names_BLW1_99, more stringent only best match
  # Sort by abundance in the data set (which abundance)
  # match among the different systems (venn diagram with four groups)
  
  # BLW1 indicative ASVs
  #asv_names_BLW1_99 %>% length()
  #colnames(ISS.rob.comp.avg) %in% asv_names_BLW1_99 %>% sum() # The majority of the ASVs was present present in the rarefied dataset!
  
  path ="Output/Indicator_species_analysis/"
  
  # Analysis based on ASV level
  data.ISS = ISS.rob.comp.avg %>% t() %>% as.data.frame()
  data.asv = asv.rob.comp.avg %>% t() %>% as.data.frame()
  
  # Analysis based on genus level
  data.tax =tax.class.rob # The ASVs are the same for ISS and asv: maybe wrong!! Attention!
  name = colnames(tax.class.rob[2:ncol(tax.class.rob)]) %>% rev()
  summary =matrix(nrow = length(name), ncol =10) %>% as.data.frame()
  colnames(summary) = c("tax level", "nclassfiedASVs", "perclassfiedASVS", "topnASVs1", "topnASVs2", "topnASVs3", "nclassications",
                        "nreadtop1", "nreadstop2", "nreadstop3")
  list = list()
  
  # ATTENTION: Either select ISS or asv! ISS is the default
  
  # The following files are needed: data.tax, ISS.rob.comp.avg| asv.rob.comp.avg
  for (i in 1:length(name)) {
    data.tax$level = data.tax[, name[i]]
    tax.rob.id.genus =data.tax %>% filter(!grepl("unclassified", level))
    summary[i, 1] = name[i]
    summary[i, 2] = nrow(tax.rob.id.genus) # number of classified ASVs
    summary[i, 3] = 1-(nrow(data.tax) - nrow(tax.rob.id.genus)) / nrow(data.tax) # percentage of ASVs that are classified on the level
    summary[i,4:6] =table(tax.rob.id.genus$level) %>% as.data.frame() %>% arrange(desc(Freq)) %>% as.data.frame() %>% head(3) %>% pull(Var1) %>% as.vector() # Classifcation with most ASVs
    summary[i, 7] = nrow(table(tax.rob.id.genus$level) %>% as.data.frame() %>% arrange(desc(Freq))) # Number of classifcations on the level
    #ISS.rob.av.comp.id.genus = ISS.rob.comp.avg[,rownames(tax.rob.id.genus)] # ISS
    ISS.rob.av.comp.id.genus = asv.rob.comp.avg[,rownames(tax.rob.id.genus)] # ASV
    ISS.rob.av.comp.id.genus = merge(tax.rob.id.genus, t(ISS.rob.av.comp.id.genus), by = 0)
    rownames(ISS.rob.av.comp.id.genus) <- ISS.rob.av.comp.id.genus$Row.names
    ISS.rob.av.comp.id.genus$Row.names <- NULL
    result <- aggregate(. ~ level, data = ISS.rob.av.comp.id.genus[,8:ncol(ISS.rob.av.comp.id.genus)], sum)
    rownames(result) <- result$level
    result$level <- NULL
    summary[i, 8:10] = rowSums(result) %>% sort() %>% as.data.frame() %>% tail(3) %>% rownames() # Classication with the highest number of reads
    list[[i]] = result
  }
  # ISS
  data.ISS.species = list[[1]] 
  data.ISS.genus = list[[2]] 
  data.ISS.family = list[[3]]
  
  # ASV
  data.asv.species = list[[1]]
  data.asv.genus = list[[2]]
  data.asv.family = list[[3]]
  
  list = list(nASV = NULL, nBLW1 = NULL, min2 = NULL, min3 = NULL, all = NULL)
  ## Beta diversity on genus and family level------
  
  # On genera level
  data.ISS.genus %>% dim()
  dist <-vegdist(data.ISS.genus %>% t(), method = "bray")
  
  # On family level
  data.ISS.family %>% dim()
  dist <-vegdist(data.ISS.family %>% t(), method = "bray")
  
  # Only load one or the other
  nmds <- metaMDS(dist)
  nmdspoints = nmds$points %>% as.data.frame()
  colnames(nmdspoints) = c("nmds1", "nmds2")
  nmdspoints$treatment = rownames(nmdspoints)
  df_BLW2.beta <- merge(df_BLW2.r, nmdspoints , by ="treatment")
  stress <- round(nmds$stress,3)
  
  # composting site
  find_hull <- function(df_BLW2.beta) df_BLW2.beta[chull(df_BLW2.beta$nmds1, df_BLW2.beta$nmds2), ]
  hulls <- plyr::ddply(df_BLW2.beta, "site_ID", find_hull)
  
  g1 =ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2))+
    geom_point(aes(color = site_ID))+
    theme_classic(base_size = 14) + theme(legend.position = "none") +
    scale_fill_manual(values = site.colors.dark, name = "Composting\nsite")+
    scale_color_manual(values = site.colors.dark, name = "Composting\nsite")+
    #geom_text(x=1.7,y=-0.7,label=paste("stress =", stress), size=4)+# genus
    geom_text(x=1.3,y=-0.55, label=paste("stress =", stress), size=4)+# family
    geom_polygon(data=hulls,alpha=0.2,aes(fill=site_ID))+
    ylab("NMDS2")+xlab("NMDS1")+
    coord_fixed(ratio =1)+
    geom_text_repel(data = df_BLW2.beta, aes(x = nmds1 , y = nmds2 , label = Number, color = site_ID), show.legend = FALSE, max.overlaps =20)
  
  color = c( "#F1F1F1","#FFC59E" ,"#E1BB4E", "yellow", "#9BB306", "#26A63A", "darkgreen")
  g2 <-ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cp))+
    theme_classic(base_size = 14) +
    geom_point(aes(fill = rel.cp), size =4, pch=21, colour ="darkgrey")+
    coord_fixed(ratio =1) +
    theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
    scale_fill_gradientn(colours = color,
                         limits = c(0,110), name= "Disease\n suppression [%]")
  g3 <-ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cu.gu))+
    theme_classic(base_size = 14) +
    geom_point(aes(fill = rel.cu.gu), size =4, pch=21, colour ="darkgrey")+
    coord_fixed(ratio =1) +
    theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
    scale_fill_gradientn(colours = color,
                         limits = c(0,110), name= "Disease\n suppression [%]")
  g4 <-ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2, color= rel.cu.rs))+
    theme_classic(base_size = 14) +
    geom_point(aes(fill = rel.cu.rs), size =4, pch=21, colour ="darkgrey")+
    coord_fixed(ratio =1) +
    theme(legend.position = "none", legend.title = element_text(size =12))+  ylab("NMDS2")+xlab("NMDS1")+
    scale_fill_gradientn(colours = color,
                         limits = c(0,110), name= "Disease\n suppression [%]")
  
  legend1 = ggpubr::get_legend(g2+ theme(legend.position = "bottom"))
  legend2 = ggpubr::get_legend(g1+ theme(legend.position = "bottom"))
  
  plots_without_legends <- plot_grid(g2, g3, g4, g1, ncol = 2, align = 'hv', labels = c("A", "B", "C", "D"))
  legends <- plot_grid(legend1, legend2, ncol = 2, align ='v')
  final_plot <- plot_grid(plots_without_legends, legends, ncol = 1, align = 'v', rel_heights = c(1, 0.1))
   #ggsave(final_plot, height =10, width =14, file ="Figures/NMDS_Betadiv_FigureA-D_genus.png")
  #ggsave(final_plot, height =10, width =14, file ="Figures/NMDS_Betadiv_FigureA-D_family.png")
  
  # Plot only the three disease assays
  ggarrange(g2, g3, g4, common.legend=T, nrow =1, ncol =3, labels  =c ("G. ultimum -cress", "G. ultimum -cucumber", "R. solani - cucumber"))
  #ggsave(file ="Bacteria_BLW2_23/BLW2_Compost_only_paper1/Figures/NMDS_Betadiv_FigureA-C.png", height = 4, width =14)
  #ggsave(file ="Figures/NMDS_Betadiv_FigureA-C_light.png", height = 4, width =14)
  
  
  # PERMANOVA
  
  # Composting site
  dist <- as.matrix(dist)
  data <- dist[!(rownames(dist) %in% c("K48", "K49")),
                            !(colnames(dist) %in% c("K48", "K49"))]
  
  adonis2(data ~ site, data = df_BLW2.r[!df_BLW2.r$site %in% c("Frick", "Spreitenbach"),], permutations = 999)
  
  # batch
  adonis2(dist ~ batch, data = df_BLW2.r, permutations = 999)
  
  # Compost age
  data = dist[!(rownames(dist) %in% c("K43", "K44")),
                           !(colnames(dist) %in% c("K43", "K44"))]
  adonis2(data  ~ age, data = df_BLW2.r %>% filter(!treatment %in% c("K43", "K44")), permutations = 999)
  
  # Disease suppression
  adonis2(dist ~ rel.cp, data = df_BLW2.r, permutations = 999)
  adonis2(dist ~ rel.cu.gu, data = df_BLW2.r, permutations = 999)
  adonis2(dist ~ rel.cu.rs, data = df_BLW2.r, permutations = 999)
  
  rm(data, final_plot, g1, g2, g3, g4, hulls, legend1, legend2, legends, nmds.BLW2.points, 
     nmds.BLW2, plots_without_legends, color)
  
  
  ## Gu-cress, S16A-----
  top = df_BLW2.r %>% filter (treatment %in% c("K8", "K36", "K34", "K15", "K28", "K13", "K12", "K30", "K35"))
  flop = df_BLW2.r %>% filter(treatment %in% c("K23", "K26", "K21", "K20", "K18", "K19", "K9", "K48", "K7"))
  
  #top = df_BLW2.r %>% filter (Number %in% c("2", "26", "24", "8", "18", "6", "5", "20", "25"))
  # flop = df_BLW2.r %>% filter(Number %in% c("14", "17", "12", "11", "9", "10", "3", "36", "1"))
  
  path.plant = "cress"
  output_cress = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "cress", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) # with cut
  #output_cress = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "cress", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) # without cut
  #output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path, path.plant = "cress", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) # Genus level # with cut
  #output_cress = ISA_PALMA(top, flop, data.ISS = data.ISS.genus, data.asv.genus, path, path.plant = "cress", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) # Genus level # without cut
  
  all = c(output_cress$all, output_cress$min2, output_cress$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  tax_select =tax.class.rob[unique(all),] # for ASVS
  tax_select = data.asv.genus[unique(all),] # for genera
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  # Order by abundance in 9 top composts ASV level
  order =ISS.rob.comp.avg[top$treatment, rownames(tax_select_count)] %>% colSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  # Order by abundance in 9 top compost genus level
  order =data.ISS.genus[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  tax_select_count_order =tax_select_count[rownames(order),]
  # tax_select_count_order = tax_select_count_order[order(-tax_select_count_order$Freq),]
  tax_select_count_order$BLW1_best = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best$x
  tax_select_count_order$BLW1_all = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all$x
  tax_select_count_order$BLW1_best_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best.398$x
  tax_select_count_order$BLW1_all_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all.398$x
  
  #write.csv(tax_select_count_order, file = paste0(path, path.plant, "_method_comparison.csv"))
    
  ## Gu-cuc, S16B-----
  top = df_BLW2.r %>% filter (treatment %in% c("K7", "K9", "K12", "K13", "K15", "K28", "K34", "K45", "K47"))
  flop = df_BLW2.r %>% filter(treatment %in% c("K18", "K20", "K21", "K22", "K23", "K25", "K26", "K43", "K49"))

  #top = df_BLW2.r %>% filter (Number %in% c("6", "24", "8", "18", "33", "1", "35", "5", "3"))
  # flop = df_BLW2.r %>% filter(Number %in% c("17", "31", "16", "13", "11", "12", "37", "9", "14"))
  path.plant = "cugu"
  output_cugu = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant ="cugu", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE)
  #output_cugu = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant ="cugu", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE)
  #output_cugu = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path, path.plant ="cugu", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE)
  #output_cugu = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path, path.plant ="cugu", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE)
  
  all = c(output_cugu$all, output_cugu$min2, output_cugu$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  tax_select =tax.class.rob[unique(all),] # for ASVS
  tax_select = data.asv.genus[unique(all),] # for genera
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  # Order by abundance in 9 top composts for ASVS
  order =ISS.rob.comp.avg[top$treatment, rownames(tax_select_count)] %>% colSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  # Order by abundance in 9 top compost for genera
  order =data.ISS.genus[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  tax_select_count_order =tax_select_count[rownames(order),]
  #tax_select_count_order = tax_select_count_order[order(-tax_select_count_order$Freq),] # This sorts based on frequncy that was detected
  
  tax_select_count_order$BLW1_best = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best$x
  tax_select_count_order$BLW1_all = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all$x
  tax_select_count_order$BLW1_best_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best.398$x
  tax_select_count_order$BLW1_all_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all.398$x
  #write.csv(tax_select_count_order, file = paste0(path, path.plant, "_method_comparison.csv"))
  
  ## Rs-cuc, S16C-------
  #top = df_BLW2.r %>% filter (treatment %in% c("K14", "K15", "K23", "K26", "K36", "K44", "K46", "K47", "K48"))
  #flop = df_BLW2.r %>% filter(treatment %in% c("K8", "K12", "K19", "K20", "K21", "K22", "K30", "K33", "K45"))
  
  top = df_BLW2.r %>% filter (Number %in% c("32", "17", "35", "14", "34", "26", "36", "24", "19"))
  flop = df_BLW2.r %>% filter(Number %in% c("2", "11", "5", "13", "33", "20", "12", "10", "23"))
  
  
  path.plant = "curs"
  output_curs = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant ="curs", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE)
  #output_curs = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant ="curs", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE)
  #output_curs = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path, path.plant ="curs", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE)
  #output_curs = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path, path.plant ="curs", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE)
  
  all = c(output_curs$all, output_curs$min2, output_curs$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  tax_select =tax.class.rob[unique(all),] # for ASVS
  #tax_select = data.asv.genus[unique(all),] # for genera
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  # Order by abundance in 9 top composts for ASVS
  order =ISS.rob.comp.avg[top$treatment, rownames(tax_select_count)] %>% colSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  # Order by abundance in 9 top compost for genera
  order =data.ISS.genus[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  tax_select_count_order =tax_select_count[rownames(order),]
  # tax_select_count_order = tax_select_count_order[order(-tax_select_count_order$Freq),]
  tax_select_count_order$BLW1_best = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best$x
  tax_select_count_order$BLW1_all = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all$x
  tax_select_count_order$BLW1_best_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best.398$x
  tax_select_count_order$BLW1_all_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all.398$x
  
  write.csv(tax_select_count_order, file = paste0(path, path.plant, "_method_comparison.csv"))
  
  ## Not included: 4 vs 4 all systems ----
  
  top = df_BLW2 %>% filter (treatment %in% c("K34", "K46", "K47", "K15"))
  flop = df_BLW2 %>% filter(treatment %in% c("K20", "K21", "K22", "K19"))
  path.plant = "4v4"
  #output_4v4 = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "4v4", fct = "all", ncompres = 3, nperm =9999, cutoff = TRUE) # with cut
  #output_4v4 = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "4v4", fct = "all", ncompres = 3, nperm =9999, cutoff = FALSE) # without cut
  #output_4v4 = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path, path.plant = "4v4", fct = "all", ncompres = 3, nperm =9999, cutoff = TRUE) # with cut
  #output_4v4 = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path, path.plant = "4v4", fct = "all", ncompres = 3, nperm =9999, cutoff = FALSE) # without cut
  
  all = c(output_4v4$all, output_4v4$min2, output_4v4$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  #tax_select =tax[unique(all),] # for ASVS
  tax_select = data.asv.genus[unique(all),] # for genera
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  # Order by abundance in 9 top composts for ASVS
  order =ISS.rob.comp.avg[top$treatment, rownames(tax_select_count)] %>% colSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  # Order by abundance in 9 top compost for genera
  order =data.ISS.genus[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  tax_select_count_order =tax_select_count[rownames(order),]
  tax_select_count_order = tax_select_count_order[order(-tax_select_count_order$Freq),]
  tax_select_count_order$BLW1 = rownames(tax_select_count_order) %in% asv_names_BLW1_99
  # write.csv(tax_select_count_order, file = paste0(path, path.plant, "_method_comparison.csv"))
  
  ## Not included: 14 vs 6 Gu-systems -----
  
  top = df_BLW2 %>% filter (treatment %in% c("K8", "K11, K12", "K13", "K14","K15", "K28", "K34","K39", "K40", "K42", "K45", "K46", "K47"))
  flop = df_BLW2 %>% filter(treatment %in% c("K18", "K19", "K20", "K21", "K23", "K26"))
  path.plant = "14v6"
  #output_14v6 = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "14v6", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) # with cut
  #output_14v6 = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "14v6", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) # without cut
  output_14v6 = ISA_PALMA(top, flop, data.ISS.genus, data.asv.genus, path,  path.plant = "14v6", fct = "all", ncompres = 6, nperm =9999, cutoff = TRUE) # with cut
  #output_14v6 = ISA_PALMA(top, flop,data.ISS.genus, data.asv.genus, path, path.plant = "14v6", fct = "all", ncompres = 6, nperm =9999, cutoff = FALSE) # without cut
  
  all = c(output_14v6$all, output_14v6$min2, output_14v6$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  #tax_select =tax[unique(all),] # for ASVS
  tax_select = data.asv.genus[unique(all),] # for genera
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  # Order by abundance in 9 top composts for ASVS
  order =ISS.rob.comp.avg[top$treatment, rownames(tax_select_count)] %>% colSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  # Order by abundance in 9 top compost for genera
  order =data.ISS.genus[rownames(tax_select_count),top$treatment] %>% rowSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  tax_select_count_order =tax_select_count[rownames(order),]
  tax_select_count_order = tax_select_count_order[order(-tax_select_count_order$Freq),]
  tax_select_count_order$BLW1 = rownames(tax_select_count_order) %in% asv_names_BLW1_99
  write.csv(tax_select_count_order, file = paste0(path, path.plant, "_method_comparison.csv"))
  
  
  ## Comparison among pathogen-plant systems----
  
  # Initialize the list to store the results
  results_list <- list()
  
  # Define the comparison types
  comparison_types <- c("min2", "min3", "all")
  
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
        output_cress[[comparison]]
      } else if(assay == "cugu") {
        output_cugu[[comparison]]
      } else if(assay == "curs") {
        output_curs[[comparison]]
      }
    }), include_assays)
    
    return(comparison_list)
  }
  
  # Function to process overlap
  process_overlap <- function(x, comparison, assay_combination) {
    
    # For ASVS
    overlapp <- tax.class.rob[Reduce(intersect, x),]
    overlapp$ASV <- rownames(overlapp)
    
    # For genus
    #overlapp <- Reduce(intersect, x) %>% as.data.frame()
    
    overlapp$comparison <- comparison
    
    # Add the assay combination as a string (e.g., "cress_cugu")
    overlapp$assay_combination <- paste(assay_combination, collapse = "_")
    
    # Generate and save the Venn diagram
    venn_plot <- ggVennDiagram(x) # Modify this based on your actual Venn diagram function
    filename <- paste0(path, "Venndiagram_", overlapp$comparison, "_", overlapp$assay_combination, ".png")
    ggsave(filename = filename, plot = venn_plot, height = 5, width = 5)
    
    return(overlapp)
  }
  
  path = "Output/Indicator_species_analysis/ASV/"
  # Loop over each comparison type
  for (comparison in comparison_types) {
    # Perform comparisons for all assay combinations
    for (assay_pair in assay_combinations) {
      comparison_list <- create_comparison_list(assay_pair, comparison)
      results_list <- c(results_list, list(process_overlap(comparison_list, comparison, assay_pair)))
    }
  }
  
  # Combine all results into a single data frame
  rm(final_results)
  final_results <- do.call(rbind, results_list)
  # Only for ASVs
  final_results$BLW1 <- final_results$ASV %in% ASVs.names.blast.BLW1.all$x
  # Only for genus
  colnames(final_results) = c("Genus", "comparison", "assay_combination")
  final_results$BLW1 <- final_results$Genus %in% c("Algoriphagus", "Pirellula", "Pseudoxanthomonas", "Sphingopyxis","Ureibacillus",
                                                 "Cytophaga", "IMCC26134", "Thermogutta", "Caldalkalibacillus")
    
  #write.csv(final_results, file =paste0(path,"comparisons_path_plants_systems.csv"), row.names = FALSE)

 # Version from the fungi
  path2 <- paste0(path, "ASV/")
  level <- "ASV"
  cutoff =TRUE
  cress <- read.csv(file =paste0(path2, "cress", ifelse(cutoff ==TRUE, "_cut","" ), "_method_comparison.csv"), sep =";")
  rownames(cress) <-cress[,1]; cress[,1] <-NULL
  
  cugu <- read.csv(file =paste0(path2, "cugu", ifelse(cutoff ==TRUE, "_cut","" ),"_method_comparison.csv"), sep =";")
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
      overlapp <- tax.class.rob[Reduce(intersect, x),]
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
  
  # Only for ASVs
  final_results$BLW1 <- final_results$ASV %in% ASVs.names.blast.BLW1.all$x
  #write.csv(final_results, file =paste0(path2,"comparisons", ifelse(cutoff ==TRUE, "_cut","" ), "_path_plants_systems.csv"), row.names = FALSE)
  
  
  
## Not included: Comparison 14v6 and 4v4----
  # For all, min3, min2
  
  assay_combinations <- list(
    c("cress", "cugu", "14v6"), 
    c("cress", "curs", "cugu", "4v4"))
  
  # Function to create comparison list
  create_comparison_list <- function(include_assays, comparison) {
    comparison_list <- setNames(lapply(include_assays, function(assay) {
      # Fetch the correct dataset based on the assay name
      if(assay == "cress") {
        output_cress[[comparison]]
      } else if(assay == "cugu") {
        output_cugu[[comparison]]
      } else if(assay == "curs") {
        output_curs[[comparison]]
      } else if (assay == "14v6") {
        output_14v6[[comparison]]
      } else if (assay == "4v4") {
        output_4v4[[comparison]]
      }
    }), include_assays)
    
    return(comparison_list)
  }
  
  
  # Loop over each comparison type
  for (comparison in comparison_types) {
    # Perform comparisons for all assay combinations
    for (assay_pair in assay_combinations) {
      comparison_list <- create_comparison_list(assay_pair, comparison)
      results_list <- c(results_list, list(process_overlap(comparison_list, comparison, assay_pair)))
    }
  }
  rm(final_results)
  final_results <- do.call(rbind, results_list)
  # Only for ASVs
  #final_results$BLW1 <- final_results$ASV %in% asv_names_BLW1_99 
  # Only for genus
  colnames(final_results) = c("Genus", "comparison", "assay_combination")
  final_results$BLW1 <- final_results$Genus %in% c("Algoriphagus", "Pirellula", "Pseudoxanthomonas", "Sphingopyxis","Ureibacillus",
                                                   "Cytophaga", "IMCC26134", "Thermogutta", "Caldalkalibacillus") 
  
  final_results =final_results %>% filter(grepl( "4v4", assay_combination) | grepl("14v6", assay_combination))
  #write.csv(final_results, file =paste0(path,"comparisons_14v6_4v4.csv"), row.names = FALSE)
  
  
  ## Comparison ASV and genus level----
  # Ony consider at least present with 3 methods (feq =2) 
  comp = list()
  system = list()
  ASV_list = list()
  path.plant = c("cress", "cugu", "curs")
  # path.plant = c("cress", "cugu", "curs", "4v4", "14v6")
  for (i in 1:3 ) {
    #ASV <- read.csv(file=paste0(path, "CUT/", path.plant[i], "_cut_method_comparison.csv"), sep =";")
    ASV = read.csv(file = paste0("Output/Indicator_species_analysis/ASV/", path.plant[i], "_method_comparison.csv"), sep =";")
    ASV = ASV %>% filter(Freq >= 2)
    ASV_list[[path.plant[i]]] <-  ASV
    #GENUS = read.csv(file=paste0(path,"GENUS/CUT/",  path.plant[i], "_genus_cut_method_comparison.csv"), sep =";")
    GENUS = read.csv(file=paste0("Output/Indicator_species_analysis/GENUS/",  path.plant[i], "_method_comparison.csv"), sep =";")
    GENUS = GENUS %>% filter(Freq >= 2)
    system[[path.plant[i]]] <- unique(GENUS)
    # Comparison
    comp[[path.plant[i]]] <- intersect(unique(ASV$Genus), GENUS$Genus)
  }
  # Venn diagram of genus which were indicative on genus an ASV level
  x <- list(gucress = comp[[1]], gucuc = comp[[2]], rscuc = comp[[3]])
  ggVennDiagram(x)
  intersect(comp[[1]], comp[[2]])
  intersect(comp[[1]], comp[[3]])
  intersect(comp[[2]], comp[[3]])
  intersect(comp[[3]],  intersect(comp[[1]], comp[[2]]))

  
  # Three most abundant for each system. must be adapted!
  #cress =system[[1]]
  #rownames(cress) = cress$Genus
  #cress[,c("K14", "K15", "K23", "K26", "K36", "K44", "K46", "K47", "K48")] %>% rowSums %>%
  #  sort(decreasing =TRUE) %>% as.data.frame() %>% rownames()

  
# Figure5C (old version) Indicative ASVs for GU-cress and GU-cuc----
  
  # Create function that generates the plot based on ASV input and fasta file
  cress <- read.csv(file = paste0("Output/Indicator_species_analysis/ASV/cress_cut_method_comparison.csv"), sep =";")%>%
    filter(Freq > 1)
  cugu <- read.csv(file = paste0("Output/Indicator_species_analysis/ASV/cugu_cut_method_comparison.csv"), sep =";")%>%
    filter(Freq > 1)
  curs <- read.csv(file = paste0("Output/Indicator_species_analysis/ASV/curs_cut_method_comparison.csv"), sep =";")%>%
    filter(Freq > 1)
  
  colors <- data.frame(color = c(
    "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "lightblue",
    "#911EB4", "#42D4F4", "#F032E6", "#BFEF45", "#FABEBE",
    "coral4", "#E6BEFF","mistyrose" , "#FFFAC8", "#800000",
    "#AAFFC3", "#808000", "#A9A9A9", "#000075"),
    phlya = c(unique(cress$Phyla), unique(cugu$Phyla), unique(curs$Phyla)) %>% unique() %>% sort() )

  sequences <- read.fasta(file = "ind_ASV_blast_search/8.all.ASV_metaxa.fasta", seqtype = "DNA")
  
  phylo_tree_ASV <- function(ASV, sequences, clade_colors, all =FALSE){
    if(all == FALSE){ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species, BLW1_all)}
  else{ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species, BLW1_all, assay_combination)}
  
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
    lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"],
    BLW1 = ASV[match(names(subset_sequences), ASV$ASV), "BLW1_all"])}
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

    p1 =p +ggnewscale::new_scale_colour()  + geom_text2(aes(label = lowest, angle = angle2, hjust = hjust, color = BLW1), data = p_data, size = 3)+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      theme(plot.margin = margin(2, 2, 2, 2, "cm"),
            axis.line = element_blank(),  # Remove axis lines
            axis.text = element_blank(),  # Remove axis text
            axis.ticks = element_blank(),  # Remove axis ticks
            panel.grid = element_blank(),  # Remove grid lines
            legend.position = "none")

  return(list(phylo_tree, p1, legend1))
  }
  phylo_tree_ASV_flat <- function(ASV, sequences, clade_colors, all =FALSE, size =2, xlim.factor = 1.35){
    if(all == FALSE){ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species, BLW1_all)}
    else{ASV <- ASV %>% dplyr::select(ASV, Phyla, Class, Order, Family, Genus, Species, BLW1_all, assay_combination)}
    
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
        lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"],
        BLW1 = ASV[match(names(subset_sequences), ASV$ASV), "BLW1_all"])}
    else{
      metadata <- data.frame(
        ASV = names(subset_sequences),  # Ensure this matches the sequence IDs
        Phylum = ASV[match(names(subset_sequences), ASV$ASV), "Phyla"],
        lowest = ASV[match(names(subset_sequences), ASV$ASV), "lowest"],
        BLW1 = ASV[match(names(subset_sequences), ASV$ASV), "BLW1_all"],
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
    
    p1 =p +ggnewscale::new_scale_colour()  + geom_tiplab(aes(label = lowest, color = BLW1), data = p_data, size = size, offset =0.5, hjust =0, align =TRUE, linetype =NA)+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      theme(axis.line = element_blank(),  # Remove axis lines
            axis.text = element_blank(),  # Remove axis text
            axis.ticks = element_blank(),  # Remove axis ticks
            panel.grid = element_blank(),  # Remove grid lines
            legend.position = "none")+
      xlim(NA, max(p$data$x) * xlim.factor)
    
    return(list(phylo_tree, p1, legend1))
  }
  
  clade_colors = colors %>% filter(phlya %in% unique(cress$Phyla)) %>% dplyr::select(color)
  phylo_tree_ASV(cress, sequences, clade_colors, all=FALSE)
  phylo_tree_cress <-phylo_tree_ASV_flat(cress, sequences, clade_colors, all=FALSE)
  #ggsave("Figures/GU_cress_phylogenetical_tree.png", height = 10, width =12, unit ="in")
  
  clade_colors = colors %>% filter(phlya %in% unique(cugu$Phyla)) %>% dplyr::select(color)
  phylo_tree_ASV(cugu, sequences, clade_colors)
  phylo_tree_cugu <-phylo_tree_ASV_flat(cugu, sequences, clade_colors, all=FALSE)
  #ggsave("Figures/GU_cuc_phylogenetical_tree.png", height = 10, width =12, unit ="in")
  
  clade_colors = colors %>% filter(phlya %in% unique(curs$Phyla)) %>% dplyr::select(color)
  phylo_tree_ASV(curs, sequences, clade_colors)
  phylo_tree_curs <-phylo_tree_ASV_flat(curs, sequences, clade_colors, all=FALSE)
 #ggsave("Figures/Rs_cuc_phylogenetical_tree.png", height = 10, width =12, unit ="in")
 
  # Only ASVs that were at least indicative in 2/3 systems
  
  all <- read.csv(file = paste0("Output/Indicator_species_analysis/ASV/comparisons_cut_path_plants_systems.csv"), sep =";") %>%
    filter(comparison %in% c(2))
  
  ASV_present_all <-all %>% filter(assay_combination == "cress_cugu_curs") %>% dplyr::select(ASV) 
  ASV_all = all %>% filter(ASV %in% ASV_present_all$ASV & assay_combination == "cress_cugu_curs")
  all =all %>% filter(!ASV %in% ASV_present_all$ASV)
  all = rbind(all, ASV_all)
  all$BLW1_all = all$BLW1
  all_BLW1 <-NULL
  
  # Add the classification indicator again
  
  # Define prefixes
  prefixes <- c(Class = "c_", Order = "o_", Family = "f_", Genus = "g_", Species = "s_")
  
  # Apply prefixes conditionally: only if the entry does not contain "unclassified"
  all <- all %>%
    mutate(across(all_of(names(prefixes)), 
                  ~ ifelse(grepl("unclassified", ., ignore.case = TRUE), ., paste0(prefixes[cur_column()], .))))
  
  
  clade_colors = colors %>% filter(phlya %in% unique(all$Phyla)) %>% dplyr::select(color)
  phylo_tree_all <- phylo_tree_ASV_flat(all, sequences, clade_colors, all=TRUE, size =3.5, xlim.factor = 1.8)
  
  # Mit GU-cress/RS-cuc§
  #ggsave("Figures/All_phylogenetical_tree.png", height = 8, width =10, unit ="in")
  
  # Without GU-cress/RS-cuc
  all_red <- all %>% filter(assay_combination %in% c("cress_cugu", "cugu_curs", "cress_cugu_curs"))
  clade_colors = colors %>% filter(phlya %in% unique(all_red$Phyla)) %>% dplyr::select(color)
  phylo_tree <- phylo_tree_ASV(all_red, sequences, clade_colors, all =TRUE)
  phylo_tree_all <- phylo_tree_ASV_flat(all_red, sequences, clade_colors, all=TRUE, size =3.5, xlim.factor = 1.8)
    
  p1 <- phylo_tree[[2]]
  p2 <-ggarrange(p1, phylo_tree[[3]], ncol =2, nrow =1, widths = c(0.85, 0.15))
  # ggsave(p2, filename= "Figures/All_phylogenetical_tree_red.png", height =8, width =9, unit ="in")
  
  ## Venn-diagram ASV and genus among pathogen levels-----
  
  # Venn.diagram ASV level

  venn.plot <- venn.diagram(
    x = list(cress = cress$ASV, cugu = cugu$ASV, curs = curs$ASV),
    #category.names = c("cress" = paste("GU-cress \n",nrow(cress), "ASVs" ),
    #                   "cugu" = paste("GU-cucumber \n", nrow(cugu),"ASVs"),
    #                   "curs" = paste("RS-cucumber \n", nrow(curs), "ASVs")),
    category.names = c("cress" = " ",
                       "cugu" = " ",
                       "curs" = " "),
    filename = NULL,  # Display on screen
    col = color.dark,
    fill = color.light,
    cat.col = "black",
    cat.cex = 1.5,
    cex = 1.5,
    #cat.pos = c(360-45, 45, 180),  # Adjust these values to change label positions
    #cat.dist = c(0.15, 0.15, 0.15),  # Adjust these values to change label distances
    #margin = 0.1,
    resolution = 300,  # Increase resolution for better clarity
    cat.fontface = "bold"
  )
  venn_ASV <- grid.grabExpr(grid.draw(venn.plot))
  
  venn_ASV <- ggplotify::as.ggplot(venn_ASV) +
    theme(
      aspect.ratio = 1
    )
  
  # Venn.diagram on genus level
  cress_g <- read.csv(file = paste0("Output/Indicator_species_analysis/Genus//cress_method_comparison.csv"), sep =";")%>%
    filter(Freq > 1)
  cugu_g <-read.csv(file = paste0("Output/Indicator_species_analysis/Genus/cugu_method_comparison.csv"), sep =";")%>%
    filter(Freq > 1)
  curs_g <- read.csv(file = paste0("Output/Indicator_species_analysis/Genus/curs_method_comparison.csv"), sep =";")%>%
    filter(Freq > 1)

  
  venn.plot <- venn.diagram(
    x = list(cress = cress_g$Genus, cugu = cugu_g$Genus, curs = curs_g$Genus),
    #category.names = c("cress" = paste("GU-cress \n",nrow(cress_g), "genera" ),
    #                   "cugu" = paste("GU-cucumber \n", nrow(cugu_g),"genera"),
    #                   "curs" = paste("RS-cucumber \n", nrow(curs_g), "genera")),
    category.names = c("cress" = " ",
                       "cugu" = " ",
                       "curs" = " "),
    filename = NULL,  # Display on screen
    col = color.dark,
    fill = color.light,
    cat.col = "black",
    cat.cex = 1.5,
    cex = 1.5,
    #cat.pos = c(360-45, 45, 180),  # Adjust these values to change label positions
    # cat.dist = c(0.15, 0.15, 0.15),  # Adjust these values to change label distances
    #margin = 0.1,
    resolution = 300,  # Increase resolution for better clarity
    cat.fontface = "bold"
  )
  venn_genus <- grid.grabExpr(grid.draw(venn.plot))
  
  venn_genus <- ggplotify::as.ggplot(venn_genus) +
    theme(
      aspect.ratio = 1
    )


  legend_df <- data.frame(
    x = 1,
    y = 1:3,
    label = c("GU-cress", "GU-cucumber", "RS-cucumber"),
    color = color.dark
  )
  
  legend_plot <- ggplot(legend_df, aes(x, y)) +
    geom_point(aes(color = color), size = 6, shape =15) +
    scale_color_identity() +
    geom_text(aes(label = label), hjust = 0, vjust = 0.5, nudge_x = 0.1, size = 3.5) +
    theme_void() +
    theme(legend.position = "none") +
    xlim(0.8, 1.8)
  
  empty_plot <- ggplot() + theme_void()
  
  
  legend <- ggarrange(empty_plot,legend_plot, empty_plot,
    ncol = 1,
    nrow = 3,
    heights =  c(0.5, 1,0.5)
  )
  
  Venn_diagramms =ggarrange(venn_ASV, venn_genus, legend, ncol =3,
            widths = c(1,1,0.5), labels = c("B", "C", ""))
  

  ggarrange(p2, Venn_diagramms, ncol=1, nrow =2,
            heights = c(1, 0.4), labels =c("A", ""))

  #ggsave("Figure5_Indicative_taxa.png", height =10, width =9)
  
  
  ### Heat map compost properties--------

  all %>% filter(Genus == "Sphingopyxis") %>% dplyr::select(ASV)
  all %>% filter(Genus == "Algoriphagus") %>% dplyr::select(ASV)
  all %>% filter(Genus == "Cellvibrio") %>% dplyr::select(ASV)
  all %>% filter(Genus == "Luteimonas_D") %>% dplyr::select(ASV)
  
  ind_bacteria = c("ASV146", "ASV467", "ASV778", "ASV2767" )
  
  asv_ind_sel <-asv.rob.comp.avg[,ind_bacteria]
  asv_ind_sel$treatment <-rownames(asv_ind_sel)
  
  df_BLW2_ind_ASV <- merge(df_BLW2.r, asv_ind_sel, by = "treatment") %>%
    dplyr::select(all_of(c(factor.continous, ind_bacteria)))
  
  # Calculate the correlations
  
  cor_results <- df_BLW2_ind_ASV %>%
    gather(key = "bacteria", value = "bacteria_value", all_of(ind_bacteria)) %>%
    gather(key = "characteristic", value = "char_value", -bacteria, -bacteria_value) %>%
    group_by(bacteria, characteristic) %>%
    dplyr::summarize(cor = stats::cor(bacteria_value[!is.na(char_value)], char_value[!is.na(char_value)], method = "spearman"),
              p_value = cor.test(bacteria_value[!is.na(char_value)], char_value[!is.na(char_value)], method = "spearman")$p.value) %>%
    ungroup()
  
  # Filter significant correlations
  cor_results <- cor_results %>%
  mutate(significant = ifelse(p_value < 0.05, "significant", "not_significant"))

  # Prepare data for heatmap
  cor_heatmap <- cor_results %>%
    mutate(color = ifelse(significant == "significant", ifelse(cor > 0, "red", "blue"), "grey"))
  
  
  ggplot(cor_heatmap, aes(x = characteristic, y = bacteria, fill = color)) +
    geom_tile(color = "white") +
    scale_fill_identity() +
    geom_text(aes(label = ifelse(significant == "significant", round(cor, 2), "")), color = "black", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Spearman Rank Correlation Heatmap",
         x = "Compost Characteristics",
         y = "Bacteria",
         fill = "Correlation")
  
  
  # Heat map for all indicative ASVs
  
  #Discard GU-cress & RS-cuc comparison
  select_ASVs <-all %>% filter(assay_combination %in% c("cress_cugu", "cugu_curs", "cress_cugu_curs")) %>% pull(ASV)
  asv_ind_sel <-asv.rob.comp.avg[, select_ASVs]
  asv_ind_sel$treatment <-rownames(asv_ind_sel)
  
  df_BLW2_ind_ASV <- merge(df_BLW2.r, asv_ind_sel, by = "treatment") %>%
    dplyr::select(all_of(c(factor.continous, select_ASVs)))
  
  # Calculate the correlations
  
  cor_results <- df_BLW2_ind_ASV %>%
    gather(key = "bacteria", value = "bacteria_value", all_of(select_ASVs)) %>%
    gather(key = "characteristic", value = "char_value", -bacteria, -bacteria_value) %>%
    group_by(bacteria, characteristic) %>%
    dplyr::summarize(cor = stats::cor(bacteria_value[!is.na(char_value)], char_value[!is.na(char_value)], method = "spearman"),
                     p_value = cor.test(bacteria_value[!is.na(char_value)], char_value[!is.na(char_value)], method = "spearman")$p.value) %>%
    ungroup()
  
  # Here the order ist still correct
  
  # Filter significant correlations
  cor_results <- cor_results %>%
    mutate(significant = ifelse(p_value < 0.05, "significant", "not_significant"))
  
  # Prepare data for heatmap
  cor_heatmap <- cor_results %>%
    mutate(color = ifelse(significant == "significant", ifelse(cor > 0, "blue", "red"), "grey"))
  
  # Order by phylogenetical tree
  tip_order <-phylo_tree_all[[2]]$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)
  
  cor_heatmap <-cor_heatmap %>% mutate(bacteria = factor(bacteria, levels = tip_order))


  heatmap_all_pc_prop <-ggplot(cor_heatmap, aes(x = characteristic, y = bacteria, fill = color)) +
    geom_tile(color = "white") +
    scale_fill_identity() +
    geom_text(aes(label = ifelse(significant == "significant", round(cor, 2), "")), color = "black", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "",
         x = "",
         fill = "Correlation")
  #ggsave(filename = "Figures/Heat_map_all_indicative_ASVs.png", height =12, width =10)
  
  ## PERMANOVA ds + addional effects----
  ISS.red <-ISS.rob.comp.avg[!rownames(ISS.rob.comp.avg) %in% c("K48", "K49"),] %>% as.matrix()
  df_BLW_red <- df_BLW2.r %>% filter(!site %in% c("Spreitenbach", "Frick"))
  
  ISS.red <-ISS.rob.comp.avg[!rownames(ISS.rob.comp.avg) %in% c("K43", "K44"),] %>% as.matrix()
  df_BLW_red <- df_BLW2.r %>% filter(!treatment %in% c("K43", "K44"))
  
  
  # Run PERMANOVA

  # GU-cress
  
  adonis_result <- adonis(ISS.rob.comp.avg  ~ rel.cp, data = df_BLW2.r, permutations = 999)
  adonis_result <- adonis(ISS.red  ~ site + rel.cp, data = df_BLW_red, permutations = 999) # site
  adonis_result <- adonis(ISS.rob.comp.avg  ~ batch + rel.cp, data = df_BLW2.r, permutations = 999) # batch
  adonis_result <- adonis(ISS.red  ~ age + rel.cp, data = df_BLW_red, permutations = 999) # age
  
  print(adonis_result$aov.tab)
  
  # GU-cuc
  adonis_result <- adonis(ISS.rob.comp.avg  ~ rel.cu.gu, data = df_BLW2.r, permutations = 999)
  adonis_result <- adonis(ISS.red  ~ site + rel.cu.gu, data = df_BLW_red, permutations = 999)
  adonis_result <- adonis(ISS.rob.comp.avg  ~ batch + rel.cu.gu, data = df_BLW2.r, permutations = 999)
  adonis_result <- adonis(ISS.red  ~ age + rel.cu.gu, data = df_BLW_red, permutations = 999) # age
  
  print(adonis_result$aov.tab)
  
  # RS-cuc
  adonis_result <- adonis(ISS.rob.comp.avg  ~ rel.cu.rs, data = df_BLW2.r, permutations = 999)
  adonis_result <- adonis(ISS.red  ~ site + rel.cu.rs, data = df_BLW_red, permutations = 999)
  adonis_result <- adonis(ISS.rob.comp.avg  ~ batch + rel.cu.rs, data = df_BLW2.r, permutations = 999)
  adonis_result <- adonis(ISS.red  ~ age + rel.cu.rs, data = df_BLW_red, permutations = 999) # age
  
  print(adonis_result$aov.tab)
  
  # When you keep age stable disease suppression of rel.cu.rs is better explained (4.2%). The other two are less good explained!
  
  # CAP analysis
  
  dissimilarity_matrix <- vegdist(ISS.red, method = "bray")
  
  # Perform CAP analysis
  cap_result <- capscale(dissimilarity_matrix ~ rel.cp + site, data = df_BLW_red, add = TRUE)
  cap_result <- capscale(dissimilarity_matrix ~ rel.cu.gu + site, data = df_BLW_red, add = TRUE)
  cap_result <- capscale(dissimilarity_matrix ~ rel.cu.rs + site, data = df_BLW_red, add = TRUE)
  
  # Extract scores for plotting
  scores <- vegan::scores(cap_result, display = "sites")
  scores_df <- as.data.frame(scores)
  
  # Add metadata for plotting
  scores_df$disease_suppression <- df_BLW_red$rel.cp
  scores_df$disease_suppression <- df_BLW_red$rel.cu.gu
  scores_df$disease_suppression <- df_BLW_red$rel.cu.rs
  
  scores_df$compost_site <- df_BLW_red$site
  
  # Plot the CAP results using ggplot2
  color <- c( "#F1F1F1","#FFC59E" ,"#E1BB4E", "yellow", "#9BB306", "#26A63A", "darkgreen")
  ggplot(scores_df, aes(x = CAP1, y = CAP2, color = disease_suppression)) +
    geom_point(aes(fill = disease_suppression), size=4, pch=21, colour ="darkgrey") +
    theme_classic(base_size = 14) +
    labs(x = "CAP1",
         y = "CAP2") +
    theme_minimal() +
    coord_fixed(ratio =1)+
    scale_fill_gradientn(colours = color,
                         limits = c(0,110), name= "Disease\n suppression [%]")

  #ggsave(filename = "Figures/CAP_site_rel_cp.png", height =4, width =4)
  #ggsave(filename = "Figures/CAP_site_rel_cugu.png", height =4, width =4)
  #ggsave(filename = "Figures/CAP_site_rel_curs.png", height =4, width =4)


# Not included: Figure Complex heat map, abundance, site-specific of indicative ASVS/genera-----

  # The following data should be loaded
  # asv.rob.comp.avg # robustly detected asvs, unrarefied
  # tax.class.rob # taxonomic classification of robustly detected ASVs
  # df_BLW2.r # meta file
  empty_plot <- ggplot() + theme_void()
  heatmap_ab_sites <-function(top, flop, asv_names, ds, tip_order){
    data.asv <-asv.rob.comp.avg %>% as.matrix() # Based on unrarefied asv table
    data.asv.TSS <- prop.table(data.asv, margin = 1) %>% as.data.frame() # Calculate relative abundance
    ind.ASV.sum <- matrix(NA, nrow = length(asv_names), ncol = 12) %>% as.data.frame() # create empty data.frame
    colnames(ind.ASV.sum)<- c("ASV", "tax", "n", "T3", "F4", "ratio.3.4", "T9", "F9", "ratio.9.9", "sites", "rho", "p.value")
    ind.ASV.sum$ASV <- asv_names
    min.value <- min(data.asv.TSS[data.asv.TSS >0])/2 # Half of the minimum abundance (to avoid zeros in the data set)
    data.asv.TSS[data.asv.TSS == 0] <- min.value
    
    ind.ASV.sum$tax <-tax.class.rob[asv_names,] %>%
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
    
    # Abundance in top 3 composts
    for (i in 1:nrow(ind.ASV.sum)) {
      ind.ASV.sum[i, "top3"] <-data.asv.TSS[  c("K34", "K46", "K47"), asv_names[i]]  %>% mean() %>% log10() %>% round(2)
    }
    
    # Abundance in flop 4 composts
    for (i in 1:nrow(ind.ASV.sum)) {
      ind.ASV.sum[i, "flop4"] <-data.asv.TSS[  c("K18", "K19", "K20", "K21"), asv_names[i]]  %>% mean() %>% log10() %>% round(2)
    }
    
    # Abundance in top 9 composts
    for (i in 1:nrow(ind.ASV.sum)) {
      ind.ASV.sum[i, "top9"] <-data.asv.TSS[top$treatment, asv_names[i]]  %>% mean() %>% log10() %>% round(2)
    }
    
    # Abundance in flop 9 composts
    for (i in 1:nrow(ind.ASV.sum)) {
      ind.ASV.sum[i, "flop9"] <-data.asv.TSS[flop$treatment, asv_names[i]]  %>% mean() %>% log10() %>% round(2)
    }
    # Calculates ratios
    ind.ASV.sum$ratio.3.4 <-(10^ind.ASV.sum$top3)/(10^ind.ASV.sum$flop4)
    ind.ASV.sum$ratio.3.4 <-ind.ASV.sum$ratio.3.4 %>% round(0)
    ind.ASV.sum$ratio.9.9<-(10^ind.ASV.sum$top9)/(10^ind.ASV.sum$flop9)
    ind.ASV.sum$ratio.9.9 <-ind.ASV.sum$ratio.9.9 %>% round(0)
    
    # present in which sites?
    for (i in 1:nrow(ind.ASV.sum)) {
      composts <-data.asv.TSS[data.asv.TSS[, asv_names[i]] >min.value,] %>% rownames()
      sites <-df_BLW2.r %>% filter(treatment %in% composts) %>% select(site_ID) %>% unique()
      ind.ASV.sum[i, "sites"] <- sites$site_ID %>% sort() %>% paste0(collapse = " ")
    }
    # Correlation test
    for (i in 1:nrow(ind.ASV.sum)) {
      cor <-cor.test(data.asv.TSS[ , asv_names[i]], df_BLW2.r[,ds], method ="spearman")
      ind.ASV.sum[i, "rho"] <- cor$estimate %>% round(2)
      ind.ASV.sum[i, "p.value"] <- cor$p.value %>% round(3)
    }
    
    ind.ASV.sum$rho.sig <-ifelse(ind.ASV.sum$p.value < 0.05 , ind.ASV.sum$rho, NA)
    
    rownames(ind.ASV.sum) <- ind.ASV.sum$ASV; ind.ASV.sum$ASV <- NULL
    
    # Sort by the order in the phylogenetical Tree
    ind.ASV.sum <- ind.ASV.sum[rev(tip_order),]
    
    # Create heat maps
    column_width <- unit(0.5, "cm")
    heatmap_data <- ind.ASV.sum[, c("top3", "flop4", "top9", "flop9")]
    heatmap_data[heatmap_data==log10(min.value) %>% round(2)] <-NA
    
    color_ramp <- colorRampPalette(c("steelblue", "darkblue"))(100)
    
    H1 <-Heatmap(
      as.matrix(ind.ASV.sum %>% select(rho.sig)),
      na_col = "grey",
      name = "Sperman's\nRho",
      col = color_ramp,
      show_row_names = TRUE,
      show_column_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 5),
      row_names_side = "left",
      width = column_width, 
      heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"))
    )
    
    
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
      column_names_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 5, fontface = "bold"),
      width = 4*column_width, 
      heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"))
    )
    color_ramp <-colorRampPalette(c("lightgreen", "darkgreen"))(100)
    
    H3 <- Heatmap(
      as.matrix(ind.ASV.sum %>% select(n)),
      name = "Number of\ncomposts",
      col = color_ramp,
      show_row_names = FALSE,
      show_column_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 5, fontface = "bold"),
      width = column_width, 
      heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"))
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
      annotation_name_gp = gpar(fontsize = 10),
      annotation_name_rot = 0, 
      show_annotation_name = TRUE,
      show_legend = c("bar" = FALSE))
    
    heat <-H1+H3+H2 +left_annotation 
    
    heatmap_ggplot <- draw(heat, heatmap_legend_side = "right")
    heatmap_grob <- grid::grid.grabExpr(draw(heatmap_ggplot))
    return(heatmap_grob)
  }
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
    
    ind.ASV.sum$tax <-tax.class.rob[asv_names,] %>%
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
      sites <-df_BLW2.r %>% filter(treatment %in% composts) %>% select(site_ID) %>% unique()
      ind.ASV.sum[i, "sites"] <- sites$site_ID %>% sort() %>% paste0(collapse = " ")
    }
    
    # Correlation test
    rho.name <- c("rho_cress", "rho_cugu", "rho_curs")
    p.name <- c("p_cress", "p_cugu", "p_curs")
    p.adj <- c("padj_cress", "padj_cugu", "padj_curs")
    rho.sig <- c("rho.sig.cress", "rho.sig.cugu", "rho.sig.curs")
    
    for (j in 1:3) {
    for (i in 1:nrow(ind.ASV.sum)) {
      cor <-cor.test(data.asv.TSS[ , asv_names[i]], df_BLW2.r[,ds[j]], method ="spearman")
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

  ## GU-cress------
  top <- df_BLW2.r %>% filter (treatment %in% c("K8", "K36", "K34", "K15", "K28", "K13", "K12", "K30", "K35"))
  flop <- df_BLW2.r %>% filter(treatment %in% c("K23", "K26", "K21", "K20", "K18", "K19", "K9", "K48", "K7"))
  asv_names <-cress$ASV
  ds <- "rel.cp"
  tip_order <- phylo_tree_cress[[2]]$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)
  
  heatmap_cress <-heatmap_ab_sites(top, flop, asv_names, ds, tip_order)

  phylo_tree_cress_l <-ggdraw(phylo_tree_cress[[2]]+ draw_plot(phylo_tree_cress[[3]], 5, 110, 0, 0))
  combined_plot1 <- cowplot::plot_grid(phylo_tree_cress_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.03))
  combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap_cress, nrow =2, ncol =1, rel_heights = c(0.001, 1))
  combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                       rel_widths = c(1,1))
  print(combined_plot3)
  #ggsave(filename = "Figures/Gu_cress_phylogenetical_tree_heat_map.png", height =10, width =8.1)
  
  ## GU-cuc-----
  top = df_BLW2.r %>% filter (treatment %in% c("K7", "K9", "K12", "K13", "K15", "K28", "K34", "K45", "K47"))
  flop = df_BLW2.r %>% filter(treatment %in% c("K18", "K20", "K21", "K22", "K23", "K25", "K26", "K43", "K49"))
  asv_names <-cugu$ASV
  ds <- "rel.cu.gu"
  tip_order <- phylo_tree_cugu[[2]]$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)
  
  heatmap_cugu <-heatmap_ab_sites(top, flop, asv_names, ds, tip_order)
  
  phylo_tree_cugu_l <-ggdraw(phylo_tree_cugu[[2]]+ draw_plot(phylo_tree_cugu[[3]], 5, 110, 0, 0))
  combined_plot1 <- cowplot::plot_grid(phylo_tree_cugu_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.03))
  combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap_cugu, nrow =2, ncol =1, rel_heights = c(0.001, 1))
  combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                       rel_widths = c(1,1))
  print(combined_plot3)
  #ggsave(filename = "Figures/Gu_cuc_phylogenetical_tree_heat_map.png", height =10, width =8.1)
  
  ## RS-cuc------
  top = df_BLW2.r %>% filter (treatment %in% c("K14", "K15", "K23", "K26", "K36", "K44", "K46", "K47", "K48"))
  flop = df_BLW2.r %>% filter(treatment %in% c("K8", "K12", "K19", "K20", "K21", "K22", "K30", "K33", "K45"))
  asv_names <-curs$ASV
  ds <- "rel.cu.rs"
  tip_order <- phylo_tree_curs[[2]]$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)
  
  heatmap_curs <- heatmap_ab_sites(top, flop, asv_names, ds, tip_order)
  
  phylo_tree_curs_l <-ggdraw(phylo_tree_curs[[2]]+ draw_plot(phylo_tree_curs[[3]], 5, 110, 0, 0))
  combined_plot1 <- cowplot::plot_grid(phylo_tree_curs_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.02))
  combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap_curs, nrow =2, ncol =1, rel_heights = c(0.001, 1))
  combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                       rel_widths = c(1,1))
  print(combined_plot3)
  #ggsave(filename = "Figures/Rs_cuc_phylogenetical_tree_heat_map.png", height =15, width =8.1)

  ## All three----
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

  heatmap <-heatmap_all(topflop, asv_names, ds, tip_order, presentation =TRUE, size=10, adjustcor = TRUE)
  phylo_tree_all_l <-ggdraw(phylo_tree_all[[2]]+ draw_plot(phylo_tree_all[[3]],2, 35, 0, 0))
  
  
  combined_plot1 <- cowplot::plot_grid(phylo_tree_all_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.07)) # Paper
  combined_plot1 <- cowplot::plot_grid(phylo_tree_all_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.015)) # Presentation
  
  combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap[[1]], nrow =2, ncol =1, rel_heights = c(0.002, 1))
  combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                       rel_widths = c(1.1,1)) # paper
  combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                       rel_widths = c(1.3,1)) # Presentation
  print(combined_plot3)
  #ggsave(filename = "Figures/All_phylogenetical_tree_heat_map.png", height =8, width =11) # Paper
  #ggsave(filename = "Figures/All_phylogenetical_tree_heat_map_1024.png", height =10, width =8) # Paper Version 27.08.24
  #ggsave(filename = "Figures/All_phylogenetical_tree_heat_map_presentation.png", height =8, width =9) # Presentation
  
  #Heatmap physcio-chemical properties
  combined_plot1 <- cowplot::plot_grid(phylo_tree_all_l, empty_plot, nrow =2, ncol =1, rel_heights = c(1, 0.1))
  combined_plot2 <- cowplot::plot_grid(empty_plot, heatmap_all_pc_prop, nrow =2, ncol =1, rel_heights = c(0.002, 1))
  combined_plot3 <- cowplot::plot_grid(combined_plot1, combined_plot2, ncol = 2,
                                       rel_widths = c(1.1,0.9))
  
  #ggsave(filename = "Figures/All_phylogenetical_tree_heat_map_prop1024.png", height =10, width =10)
  
  
  ## Taxonomic summary association study-------
  
  cress$Phyla %>% table() %>% sort() %>% as.data.frame()
  cress$Family %>% table() %>% sort() %>% as.data.frame()
  cress$Genus %>% table() %>% sort() %>% as.data.frame()
  
  cugu$Phyla %>% table() %>% sort() %>% as.data.frame()
  cugu$Family %>% table() %>% sort() %>% as.data.frame()
  cugu$Genus %>% table() %>% sort() %>% as.data.frame()
  
  curs$Phyla %>% table() %>% sort() %>% as.data.frame()
  curs$Family %>% table() %>% sort() %>% as.data.frame()
  curs$Genus %>% table() %>% sort() %>% as.data.frame()
  
  all %>% filter(assay_combination == "cress_cugu") %>% pull(Phyla) %>% table() %>% sort() %>% as.data.frame()
  all %>% filter(assay_combination == "cress_cugu") %>% pull(Family) %>% table() %>% sort() %>% as.data.frame()
  all %>% filter(assay_combination == "cress_cugu") %>% pull(Genus) %>% table() %>% sort() %>% as.data.frame()
  
  all %>% filter(assay_combination == "cugu_curs") %>% pull(Phyla) %>% table() %>% sort() %>% as.data.frame()
  all %>% filter(assay_combination == "cugu_curs") %>% pull(Family) %>% table() %>% sort() %>% as.data.frame()
  all %>% filter(assay_combination == "cugu_curs") %>% pull(Genus) %>% table() %>% sort() %>% as.data.frame()
  
  all %>% filter(assay_combination == "cress_cugu_curs") %>% pull(Family) %>% table() %>% sort() %>% as.data.frame()
  
  
# Not included: High quality composts------
  
  df_BLW.quality <-df_BLW2.r %>% filter(pH < 8.2 & OD550 < 0.6 & sal <20 & DS >50 & Nmin >100 & Corg < 50 & NO3 > 80 & NH4 < 200 & NO3.Nmin > 0.4 & NO2 < 20)
  
  # Assuming df_BLW2.r is your original data frame
  original_count <- nrow(df_BLW2.r)
  
  # Define the conditions
  conditions <- list(
    "pH < 8.2" = df_BLW2.r$pH < 8.2,
    "OD550 < 0.6" = df_BLW2.r$OD550 < 0.6,
    "sal < 20" = df_BLW2.r$sal < 20,
    "DS > 50" = df_BLW2.r$DS > 50,
    "Nmin > 100" = df_BLW2.r$Nmin > 100,
    "Corg < 50" = df_BLW2.r$Corg < 50,
    "NO3 > 80" = df_BLW2.r$NO3 > 80,
    "NH4 < 200" = df_BLW2.r$NH4 < 200,
    "NO3.Nmin > 0.4" = df_BLW2.r$NO3.Nmin > 0.4,
    "NO2 < 20" = df_BLW2.r$NO2 < 20
  )
  
  # Initialize a vector to store the number of rows remaining after each filter
  remaining_counts <- sapply(conditions, function(condition) {
    sum(condition, na.rm = TRUE)
  })
  
  reduction_counts <- original_count - remaining_counts
  max_reduction_filter <- names(reduction_counts)[which.max(reduction_counts)]
  print(reduction_counts)

  df_BLW.noqu <- df_BLW2.r %>% filter(!treatment %in% df_BLW.quality$treatment)
  
  df_BLW2.r$quality <- ifelse(df_BLW2.r$treatment %in% df_BLW.quality$treatment, 1, 0) %>% as.factor()
  
  # Differences in disease suppression low and high quality composts
  g1 <- ggplot(data = df_BLW2.r, aes(x = quality, y = rel.cp)) + geom_boxplot()
  g2 <- ggplot(data = df_BLW2.r, aes(x = quality, y = rel.cu.gu)) + geom_boxplot()
  g3 <- ggplot(data = df_BLW2.r, aes(x = quality, y = rel.cu.rs)) + geom_boxplot()
  
  ANOVA <-aov(rel.cp ~ quality, data = df_BLW2.div)
  ANOVA <-aov(rel.cu.gu ~ quality, data = df_BLW2.div)
  ANOVA <-aov(rel.cu.rs ~ quality, data = df_BLW2.div)
  summary(ANOVA)
  
  
  ggarrange(g1, g2, g3, nrow=1, ncol =3, labels = c("GUcress", "GUcuc", "RScuc"))
  # ggsave(filename = "Quality composts/Differnces_ds_compost_quality.png", height =4, width =9)
  
  qcomp <-df_BLW2.r %>% filter(quality ==1) %>% pull(treatment)
  
  df_BLW2.div.q <-df_BLW2.div %>% filter(treatment %in% qcomp) 
  
  # Correlation with alpha diversity
  df_BLW2.div$quality <- ifelse(df_BLW2.div$treatment %in% df_BLW.quality$treatment, 1, 0) %>% as.factor()
  
  ANOVA <-aov(mean.sobs ~ quality, data = df_BLW2.div)
  summary(ANOVA)
  
  cor.test(df_BLW2.div.q$mean.sobs, df_BLW2.div.q$rel.cp, method ="spearman")
  cor.test(df_BLW2.div.q$mean.sobs, df_BLW2.div.q$rel.cu.gu, method ="spearman")
  cor.test(df_BLW2.div.q$mean.sobs, df_BLW2.div.q$rel.cu.rs , method = "spearman")
  
  # Correlation with beta diversity
  
  data <- ISS.bray.avg.comp[(rownames(ISS.bray.avg.comp) %in% df_BLW.quality$treatment),
                            (colnames(ISS.bray.avg.comp) %in% df_BLW.quality$treatment)]
  
  
  
  adonis2(data ~ rel.cp, data = df_BLW.quality, permutations = 999)
  adonis2(data ~ rel.cu.gu, data = df_BLW.quality, permutations = 999)
  adonis2(data ~ rel.cu.rs, data = df_BLW.quality, permutations = 999)
  
  adonis2(ISS.bray.avg.comp ~quality, data =df_BLW2.r, permutations =999)
  
  # NMDS
  
  nmds.BLW2 <- metaMDS(ISS.bray.avg.comp)
  nmds.BLW2.points = nmds.BLW2$points %>% as.data.frame()
  colnames(nmds.BLW2.points) = c("nmds1", "nmds2")
  nmds.BLW2.points$treatment = rownames(nmds.BLW2.points)
  df_BLW2.beta <- merge(df_BLW2.r, nmds.BLW2.points , by ="treatment")
  nmds.BLW2$stress # 0.10
  
  ggplot(data = df_BLW2.beta, aes(x = nmds1, y = nmds2))+
    geom_point(aes(color = quality))+
    theme_classic(base_size = 14) + theme(legend.position = "right") +
    geom_text(x=0.6,y=-0.2,label='stress = 0.1', size=4)+
    ylab("NMDS2")+xlab("NMDS1")+
    coord_fixed(ratio =1)+
    geom_text_repel(data = df_BLW2.beta, aes(x = nmds1 , y = nmds2 , label = Number), show.legend = FALSE, max.overlaps =20)
  #ggsave(filename = "Figures/NMDs_Betadiv_quality_composts.png", height = 8, width =12)
  
  
  # Also when only considering the quality compost we cannot explain differences in disease suppression with the bacterial community strucutre
  
  # Indicator species analysis
  
  # CUGU
  top = df_BLW2.r %>% filter (Number %in% c(1,6,7,8,22,33,35))
  flop = df_BLW2.r %>% filter(Number %in% c(12,13,15,16,19,21,25))
  output = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "cugu", fct = "all", ncompres = 4, nperm =9999, cutoff = TRUE)
  
  # GURS
  
  top = df_BLW2.r %>% filter (Number %in% c(19, 35, 8, 26))
  flop = df_BLW2.r %>% filter(Number %in% c(13, 12, 33, 1))
  output = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "curs", fct = "all", ncompres = 3, nperm =9999, cutoff = TRUE)
  
  # CRESS
  top = df_BLW2.r %>% filter (Number %in% c(6,8,22,25))
  flop = df_BLW2.r %>% filter(Number %in% c(1,16,23,33))
  output = ISA_PALMA(top, flop, data.ISS, data.asv, path, path.plant = "cress", fct = "all", ncompres = 3, nperm =9999, cutoff = TRUE)
  
  all = c(output$all, output$min2, output$min3)
  count = all %>% table() %>% as.data.frame()
  rownames(count) = count$.; count$. <- NULL
  tax_select =tax.class.rob[unique(all),] # for ASVS
  tax_select_count = merge(count, tax_select, by =0)
  rownames(tax_select_count) <- tax_select_count$Row.names;  tax_select_count$Row.names <- NULL
  
  # Order by abundance in 9 top composts ASV level
  order =ISS.rob.comp.avg[top$treatment, rownames(tax_select_count)] %>% colSums() %>%
    sort(decreasing = TRUE) %>% as.data.frame()
  
  tax_select_count_order =tax_select_count[rownames(order),]
  # tax_select_count_order = tax_select_count_order[order(-tax_select_count_order$Freq),]
  tax_select_count_order$BLW1_best = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best$x
  tax_select_count_order$BLW1_all = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all$x
  tax_select_count_order$BLW1_best_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.best.398$x
  tax_select_count_order$BLW1_all_389 = rownames(tax_select_count_order) %in% ASVs.names.blast.BLW1.all.398$x
  
  #write.csv(tax_select_count_order,  file ="Output/Indicator_species_analysis/ASV_quality/cugu_method_comparison.csv")
  #write.csv(tax_select_count_order,  file ="Output/Indicator_species_analysis/ASV_quality/curs_method_comparison.csv")
  #write.csv(tax_select_count_order,  file ="Output/Indicator_species_analysis/ASV_quality/cress_method_comparison.csv")
  
  all  <-tax_select_count_order
  
  # Summarizing indicator species analysis
  

  folder <- "ASV_quality"
  folder <- "ASV"
  
  pathogen.plant <- "cress"
  pathogen.plant <- "cugu"
  pathogen.plant <- "curs"

  #sink(paste0("Output/Indicator_species_analysis/ASV/Summary_ASV_", pathogen.plant,".txt"))
  
  # Simple for loop which creates a txt file for each of the pathogen-plant system within seconds!
  for (i in c("all", "PBC", "ISA", "ALDEX", "MAASLIN")) {
     method <- i
    if (method == "all") {
      data <-read.csv(file=paste0("Output/Indicator_species_analysis/",folder,"/", pathogen.plant, "_method_comparison.csv"), sep =";")
      rownames(data) <- data$ASV; data$ASV <-NULL
      data <-data %>% filter(Freq >1)
      data_tax <- data
    } else {
      data <-read.csv(file=paste0("Output/Indicator_species_analysis/", folder, "/", pathogen.plant, "_cut_", method,".csv"))
      
      if (method == "MAASLIN") {
        rownames(data) <- data$feature; data$feature <-NULL
      } else{
        rownames(data) <- data$X; data$X <-NULL
      }
      
      data_tax <-merge(data, tax, by = "row.names", all.x =FALSE)
      rownames(data_tax) <- data_tax$Row.names; data_tax$Row.names <-NULL
      
      if (method %in% c("PBC", "ISA")) {
        data_tax %>% arrange(desc(stat)) %>% select(stat, Phyla, Family, Genus) %>% head()
      }
      
      if (method == "ALDEX") {
        data_tax %>% arrange(desc(effect)) %>% select(effect, Phyla, Family, Genus) %>% head() 
      }
      
      if (method =="MAASLIN"){
        data_tax %>% arrange(desc(coef)) %>% select(coef, Phyla, Family, Genus) %>% head()
      }
    }
    
    # for all comparisons
    
    paste("Method used:", method) %>% print()
    paste("Number of ASVs that were indicative:", nrow(data_tax)) %>% print()
    if(method == "all") {
      paste("Number of ASVs shared with BLW1 99%:", data %>% filter(BLW1_all== TRUE) %>% nrow()) %>% print()
      data %>% filter(BLW1_all== TRUE) %>% select(Phyla, Family, Genus) %>% print()
    }
    paste("Summary for phyla") %>% print()
    data_tax %>% pull(Phyla) %>% table() %>% as.data.frame() %>% print()
    paste("Most abundant families") %>% print()
    data_tax %>% pull(Family) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% head() %>% print()
    paste("Most abundant genera") %>% print()
    data_tax %>% pull(Genus) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% head() %>% print()
  }
  
  #sink()
  

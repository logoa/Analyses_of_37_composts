# Preperation for data analysis fungi BLW2

# Author: Anja Logo
# Started on 09.09.24
# Last changes: 24.09.24

source(file = "setup_F24.R")

# Functions-----

# Summarize the read count respecitve the number of ASVs
smry_asv_seq <-function(data){
  print(paste("Number of sequences", sum(data)))
  print(paste(" Number of ASVs", nrow(data)))
  print(paste("Number of sequences per sample"))
  print(data %>% colSums() %>% sort() %>% as.data.frame() %>% summary()) # sequences per samples
  print(paste("Number of ASVs per sample"))
  print(ifelse(data >0, 1,0) %>% colSums() %>% sort() %>% as.data.frame() %>% summary())
  
}

# Function to calculate the vegdist for the average of 100 rarefied data table for each sample
avg.vegdist = function(x, method ="bray") {
  data = merge(design.BLW2, x, by= 0) # merge with design file
  rownames(data) = data[,1] # change rownames
  data[,1] = NULL
  data2 = as.data.frame(aggregate(data[,(ncol(design.BLW2)+1):ncol(data)], list(data$treatment), mean))
  rownames(data2) <-data2[,1]
  data2[,1] <- NULL
  vegdist(data2, method = "bray")
}

# ASV & Tax file------

compost.exclude = c("K6", "K10", "K17", "K27", "K37")
df_BLW2 <- read.csv(file ="data/meta_data_compost_BLW_without_errors.csv", sep =";")
factors <- c("batch", "site", "company", "comp.system", "site_ID") # Change to factor
for (i in factors) {
  df_BLW2[, i] <- as.factor(df_BLW2[,i])
}
df_BLW2.r <- df_BLW2 %>% filter(!treatment %in% compost.exclude) # Reduce normal distribution to 37 composts
rm(factors)

# Compost_ID for changing the labeling of the data frame
compost_ID <- read.csv(file ="data/Key_treatment_compostID.csv", sep= ";")
compost_ID$treatment <- as.factor(compost_ID$treatment)
compost_ID$compost_ID <- as.factor(compost_ID$compost_ID)
compost_ID$site_ID <- as.factor(compost_ID$site_ID)
compost_ID <- compost_ID %>% filter(!treatment %in% c("K27", "K37")) # Exclude two error composts
compost_ID$treatment <-droplevels(compost_ID$treatment)

# Design
# Key for sequencing samples
design.BLW2 <- read.csv(file ="data/design_BLW2.csv")
design.BLW2$site <- as.factor(design.BLW2$site)

# Add the new numbers to the design.BLW2
design.BLW2 <- merge(design.BLW2, compost_ID, by ="treatment", ) # Add new labeling
rownames(design.BLW2) <- design.BLW2$ID
design.BLW2$ID <- NULL
design.BLW2 <- design.BLW2[order(rownames(design.BLW2)),]
design.BLW2$batch.y <-NULL
setnames(design.BLW2, "batch.x", "batch")
design.BLW2 = design.BLW2 %>% filter(!treatment %in% compost.exclude) # 37 compost

# Tax file ASV
tax <- read.table(file="Sequencing_data/8.all.ASV_ITSx.ITS2.mothur.wang.taxonomy", sep =";") # GDTB
colnames(tax) <- c("V1", "Phyla", "Class", "Order", "Family", "Genus", "Species", "V8")
tax[, c("ASV", "Kingdom")] <- str_split_fixed(tax[,1], "\t", 2)
tax <- tax[c("ASV", "Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species", "V8")]
rownames(tax) <-tax$ASV
tax<- tax[mixedorder(rownames(tax)), ] 
tax <-tax[,-1]
tax$V8 <- NULL

# Remove percentages
tax$Kingdom <- sub("\\(.*","",tax$Kingdom) # remove the percentage after the kingdom
tax$Phyla <- sub("\\(.*","",tax$Phyla)
tax$Class <- sub("\\(.*","",tax$Class)
tax$Order <- sub("\\(.*","",tax$Order)
tax$Family <- sub("\\(.*","",tax$Family)
tax$Genus <- sub("\\(.*","",tax$Genus)
tax$Species <- sub("\\(.*","",tax$Species)

# Unrarefied ASV table (for indicators species analysis)
asv <- read.table(file="Sequencing_data/9.all.ASV_map.txt", header =T)
colnames(asv) <- sub("_S[0-9]+", "", colnames(asv))
rownames(asv) <- asv$ASV;asv$ASV <- NULL
colnames(asv) <- sub("F_", "", colnames(asv)) # remove F_
asv.NTC <- asv
asv <-asv[, !grepl("^NTC", colnames(asv))]
asv<- asv[rowSums(asv)!=0,] 

# Remove non-bacterial tax.BLW1
tax.F <- tax[tax$Kingdom =="k__Fungi",] 
tax.n.F <- tax[tax$Kingdom !="k__Fungi",]
asv.F <- asv[rownames(tax.F),]
smry_asv_seq(asv.F)

# Summarize non-fungal sequences
asv.n.F <- asv[rownames(tax.n.F),]
colnames(tax.n.F) <- paste0("Tax_", colnames(tax.n.F))
tax.asv.n.F <-merge(tax.n.F, asv.n.F, by =0)
rownames(tax.asv.n.F) <- tax.asv.n.F$Row.names; tax.asv.n.F$Row.names <-NULL
n.F.sum <-tax.asv.n.F %>% 
  group_by(Tax_Kingdom) %>% 
  summarise(across(!starts_with("Tax_"), sum)) %>% as.data.frame()
rownames(n.F.sum) <- n.F.sum$Tax_Kingdom; n.F.sum$Tax_Kingdom <- NULL
n.F.sum %>% rowSums()
rm(n.F.sum, asv.n.F, tax.n.F)

# Only 37 composts
asv.r <- asv.F[, !grepl("^K10|K6|K17", colnames(asv.F))] # Reduce only to 37 composts + 4 peat substrate
asv.r <- asv.r[rowSums(asv.r)!=0,] # delete all ASVs that were only present in the discarded samples
tax.r = tax[rownames(asv.r),] # Reduced samples to 164 samples
smry_asv_seq(asv.r)

# Plot sequence per sample

sq.num <-asv.r %>% colSums() %>% as.data.frame()
colnames(sq.num) <- c("seq.number")

design.BLW2.sqn <- merge(design.BLW2, sq.num, by =0)
rownames(design.BLW2.sqn) <- design.BLW2.sqn$Row.names
design.BLW2.sqn$Row.names <- NULL

ggplot(data = design.BLW2.sqn , aes(x= treatment, y = seq.number))+ geom_boxplot()
#ggsave(filename = "figures/seqnumber_per_samples.png", height =5, width = 11)

ggplot(data = design.BLW2.sqn , aes(x= site, y = seq.number))+ geom_boxplot()
#ggsave(filename = "figures/seqnumber_per_site.png", height =5, width = 8)

design.BLW2.sqn$batch <- as.factor(design.BLW2.sqn$batch)
ggplot(data = design.BLW2.sqn , aes(x= batch, y = seq.number))+ geom_boxplot()
#ggsave(filename = "figures/seqnumber_per_batch.png", height =5, width = 8)

expected_starts <- c(Kingdom = "k", Phyla = "p", Class = "c", Order = "o", Family = "f", Genus = "g", Species = "s")
for (col_name in names(expected_starts)) {
  # Check if the column exists in your dataframe to avoid errors
  if (col_name %in% colnames(tax.r)) {
    # Get the expected start letter for the current column
    expected_letter <- expected_starts[col_name]
    
    # Loop through each row in the current column
    for (row_index in 1:nrow(tax.r)) {
      # Extract the first letter of the current cell
      current_letter <- substr(tax.r[[col_name]][row_index], 1, 1)
      
      # Check if the first letter matches the expected letter
      if (current_letter != expected_letter) {
        # If not, replace the cell value with "unclassified"
        tax.r[[col_name]][row_index] <- "unclassified"
      }}}}
rm(col_name, current_letter, expected_letter, expected_starts, row_index)

#write.table(tax.F, file ="data/tax.F.BLW2.txt") # Taxonomy only fungi
#write.table(tax.r, file ="data/tax.F.r.txt") # Taxonomy only 37 composts
#write.table(asv.F, file = "data/asv.F.BLW2.txt") # unrarefied ASVs only Fungi
#write.table(asv.r, file = "data/asv.F.r.txt") # unrarefied ASVs only Fungi, without 37 composts + peat

# Normalization----------------

# Exclude samples with not enough sequences

# ATTENTION: only run one or the other!!#######
# Without K19
high.reads <-design.BLW2.sqn %>% filter(seq.number>20000) %>% rownames()
asv.r.hr <- asv.r[,high.reads]

# With K19
high.reads <-design.BLW2.sqn %>% filter(seq.number>10000) %>% rownames()
asv.r.hr <- asv.r[,high.reads]

######################

# Rarefy abundance table 100 times
ISS.iters <- mclapply(as.list(1:100), function(x) rrarefy(t(asv.r.hr), min(colSums(asv.r.hr))), # 20396 sequences minimum
                      mc.cores = 1) # 100 rarefied abundance tables

# For alpha diversity calculations & beta diversity based on Ackermann methods
ISS.array <- laply(ISS.iters, as.matrix)
ISS <- apply(ISS.array, 2:3, median) # median of all iterations
ISS.prop = prop.table(ISS, margin = 1) * 100  # get proportions
rm(ISS.array)

# Rarefied and robustly detected ASVs only composts
robdetect <- 2
robdetect <- 3

ISS = ISS %>% t() %>% as.data.frame()
ISS01 = ifelse(ISS > 0, 1, 0) %>% t() # presence/absence
ISS01.meta <-merge(design.BLW2["treatment"], ISS01, by = 0)  # merge with design file
row.names(ISS01.meta)<-ISS01.meta$Row.names;ISS01.meta$Row.names<-NULL
ISS_01_rob <-as.data.frame(aggregate(ISS01.meta[,2:ncol(ISS01.meta)], list(ISS01.meta$treatment), median)) # calculate median

if (robdetect ==3) {
  ISS_01_rob <- ISS_01_rob %>% mutate_all(~ifelse(. == 0.5,0, .)) # Put the 2/4 (0.5) also to zero
} 

if (robdetect ==2) {
  ISS_01_rob <- ISS_01_rob %>% mutate_all(~ifelse(. == 0.5,1, .)) 
}

# Only composts
ISS_01_rob_comp = ISS_01_rob[!grepl("St", ISS_01_rob$Group.1),] # Only select composts
rownames(ISS_01_rob_comp) =ISS_01_rob_comp$Group.1 ; ISS_01_rob_comp$Group.1  <- NULL
select <-colSums(ISS_01_rob_comp) >0 # Only select ASVs that are at least robustly detected in one compost
ISS_01_rob_comp_red <-ISS_01_rob_comp[,select] # Select only robustly detected ASVs (unrarefied data)
sum(colSums(ISS_01_rob_comp_red)== 0) # check if really all the 0 are gone
select = colnames(ISS_01_rob_comp_red)
ISS.rob.comp = ISS[select,] %>% t() %>% as.data.frame()
design.BLW2.r <-design.BLW2[rownames(ISS.rob.comp),] # Subset design matrix
rownames(ISS.rob.comp) == rownames(design.BLW2.r) # Check order
ISS.rob.comp.avg = ISS.rob.comp %>% aggregate(list(design.BLW2.r$treatment), mean) # Calculate the average of the four samples
ISS.rob.comp.avg = ISS.rob.comp.avg[!grepl("Std", ISS.rob.comp.avg$Group.1), ] 
rownames(ISS.rob.comp.avg) <-ISS.rob.comp.avg[,1]; ISS.rob.comp.avg[,1] <- NULL
ISS.rob.comp <- ISS.rob.comp[!grepl("Std",rownames(ISS.rob.comp)),]

# Compost and peat
rownames(ISS_01_rob) =ISS_01_rob$Group.1 ; ISS_01_rob$Group.1  <- NULL
select <-colSums(ISS_01_rob) >0 # Only select ASVs that are at least robustly detected in one compost
ISS_01_rob_red <-ISS_01_rob[,select] # Select only robustly detected ASVs (unrarefied data)
sum(colSums(ISS_01_rob_red)== 0) # check if really all the 0 are gone
select = colnames(ISS_01_rob_red)
ISS.rob = ISS[select,] %>% t() %>% as.data.frame()
rownames(ISS.rob) == rownames(design.BLW2.r) # Check order
ISS.rob.avg = ISS.rob %>% aggregate(list(design.BLW2.r$treatment), mean) # Calculate the average of the four samples
rownames(ISS.rob.avg) <-ISS.rob.avg[,1]; ISS.rob.avg[,1] <- NULL

# Same as ISS.rob.comp.avg but with unrarefied abundances table
asv.t = asv.r.hr %>% t() %>% as.data.frame()
rownames(asv.t) == rownames(design.BLW2.r) # Check order
asv.t.avg = asv.t %>% aggregate(list(design.BLW2.r$treatment), mean) # calculate mean for the four replicates

asv.t.avg.comp = asv.t.avg[!grepl("Std", asv.t.avg$Group.1), ] # discard peat
rownames(asv.t.avg.comp) <- asv.t.avg.comp$Group.1; asv.t.avg.comp$Group.1 <-NULL
asv.rob.comp.avg = asv.t.avg.comp[, colnames(ISS.rob.comp.avg)] # select only robustly detected ASVs

ISS <- ISS[rowSums(ISS) != 0,]
ISS.comp <- ISS %>% select(starts_with("K")) # Only peat
ISS.comp <- ISS.comp[rowSums(ISS.comp) != 0,]
# Summary of different files

smry_asv_seq(ISS) # All rarefied
smry_asv_seq(ISS.comp) # Only composts
smry_asv_seq(ISS.rob %>% t()) # All robustly detected
smry_asv_seq(ISS.rob.comp %>% t())# Only composts
smry_asv_seq(ISS.rob.comp.avg %>% t()) # Average compost

# Save data files without K19
#write.table(ISS, file ="data/ISS.F.BLW2.txt") # all bacterial asvs after rarifying
#write.table(ISS.comp, file= "data/ISS.F.BLW2.comp.txt")
write.table(ISS.rob %>% t(), file = paste0("data/ISS.rob_",robdetect, ".txt")) # for beta diversiy compost and peat with replicates
write.table(ISS.rob.avg %>% t(), file=paste0("data/ISS.rob.avg_", robdetect, ".txt")) # For beta diversity compost and peat without replicates
write.table(ISS.rob.comp %>% t(), file = paste0("data/ISS.rob.comp_", robdetect,".txt")) # For beta diversity, only composts with replicates
write.table(ISS.rob.comp.avg %>% t(), file=paste0("data/ISS.rob.comp.avg_", robdetect, ".txt")) # for Indicator species analysis PCB, ISA, MaAsLin2, for beta diversiy only composts without replicates

# unrarified
write.table(asv.rob.comp.avg, file= paste0("data/asv.rob.comp.avg_", robdetect, ".txt")) # For indictor species analysis Aldex2

# Save data files with 10'000 reads
#write.table(ISS, file =paste0("data/10000_reads/ISS.F.BLW2.txt")) # all bacterial asvs after rarifying
write.table(ISS.comp, file= paste0("data/10000_reads/ISS.F.BLW2.comp.txt"))
write.table(ISS.rob %>% t(), file = paste0("data/10000_reads/ISS.rob_", robdetect,".txt")) # for beta diversiy compost and peat with replicates
write.table(ISS.rob.avg %>% t(), file=paste0("data/10000_reads/ISS.rob.avg_", robdetect, ".txt")) # For beta diversity compost and peat without replicates
write.table(ISS.rob.comp %>% t(), file = paste0("data/10000_reads/ISS.rob.comp_", robdetect,".txt")) # For beta diversity, only composts with replicates
write.table(ISS.rob.comp.avg %>% t(), file=paste0("data/10000_reads/ISS.rob.comp.avg_", robdetect, ".txt")) # for Indicator species analysis PCB, ISA, MaAsLin2, for beta diversiy only composts without replicates

# unrarified
write.table(asv.rob.comp.avg, file=paste0("data/10000_reads/asv.rob.comp.avg_", robdetect, ".txt"))# For indictor species analysis Aldex2

rm(ISS01, ISS01.meta, ISS_01_rob, ISS_01_rob_comp, select, ISS_01_rob_comp_red, ISS_01_rob_red, asv.t,
   asv.t.avg, asv.t.avg.comp, design.BLW2.sqn, sq.num, tax.asv.n.F, i)

# Calculate alpha diversity----------------------
# Observed richness
sobs <- mclapply(ISS.iters, function(x) specnumber(x), mc.cores = 1)
sobs <- apply(laply(sobs, as.matrix), 2, mean)

# Shannon diversity
shannon <- mclapply(ISS.iters, function(x) diversity(x, index = "shannon"))
shannon <- apply(laply(shannon, as.matrix), 2, mean)

# Inverse Simpson diversity
invsimpson <- mclapply(ISS.iters, function(x) diversity(x, index = "invsimpson"))
invsimpson <- apply(laply(invsimpson, as.matrix), 2, mean)

# Pielou's evenness
evenness <- shannon/log10(sobs)

# Bind to data.frame
ISS.alpha <- data.frame(cbind(sobs, evenness, shannon, invsimpson))
rm(sobs, shannon, invsimpson, evenness)
ISS.alpha$evenness = ISS.alpha$shannon/ log(ISS.alpha$sobs) # correct evenness -> put to calculation part
#write.csv(ISS.alpha, file ="data/ISS.alpha.F.BLW2.csv")
#write.csv(ISS.alpha, file ="data/withK19/ISS.alpha.F.BLW2.csv")

# Calculate mean and merge with meta file for correlation analysis without peat
rownames(ISS.alpha) == rownames(design.BLW2.r)
ISS.alpha.merged =merge(design.BLW2.r, ISS.alpha, by =0) # merge with design
ISS.alpha.mean = ISS.alpha.merged %>% group_by(treatment) %>% # Calculate averages for all variables
  dplyr::summarize(mean.sobs= mean(sobs), mean.shannon =mean(shannon),
                   mean.even = mean(evenness), mean.ivsimp = mean(invsimpson)) %>%
  as.data.frame()
asv1.alpha.BLW2 = ISS.alpha.mean[ grep("K", ISS.alpha.mean$treatment),] # Only selected compsots
asv1.alpha.BLW2$td.shannon = exp(asv1.alpha.BLW2$mean.shannon) # Calculate additional variable
df_BLW2.div <-merge(df_BLW2.r, asv1.alpha.BLW2, by = "treatment") # Fuse with meta file
#write.csv(df_BLW2.div, file ="data/df_BLW2.div.csv")
#write.csv(df_BLW2.div, file ="data/withK19/df_BLW2.div.csv")
rm(asv1.alpha.BLW2, ISS.alpha.merged, ISS.alpha.mean)


# Calculate beta diversity--------------------

rownames(ISS.iters[[1]]) == rownames(design.BLW2.r) # Check if the rownames are the same

# All the samples replicates separately

# Bray-curtis
ISS.iters.bray <- mclapply(ISS.iters, function(x) vegdist(x, method = "bray"), mc.cores = 1)
ISS.array.bray <- laply(ISS.iters.bray, as.matrix)
ISS.bray <- as.dist(apply(ISS.array.bray, 2:3, mean))
rm(ISS.iters.bray, ISS.array.bray)
#write.table(ISS.bray %>% as.matrix() %>% as.data.frame(), file ="data/ISS.bray.txt")
#write.table(ISS.bray %>% as.matrix() %>% as.data.frame(), file ="data/withK19/ISS.bray.txt")

# Jaccard
ISS.iters.jac <- mclapply(ISS.iters, function(x) vegdist(x, method = "jaccard"), mc.cores = 1)
ISS.array.jac <- laply(ISS.iters.jac, as.matrix)
ISS.jac <- as.dist(apply(ISS.array.jac, 2:3, mean))
rm(ISS.iters.jac, ISS.array.jac)
#write.table(ISS.jac %>% as.matrix() %>% as.data.frame(), file ="data/ISS.jac.txt")
#write.table(ISS.jac %>% as.matrix() %>% as.data.frame(), file ="data/withK19/ISS.jac.txt")
# Replicates together 

# Bray curtis
ISS.iters.avg.bray <- mclapply(ISS.iters, function(x, method ="bray") avg.vegdist(x), mc.cores = 1)
ISS.array.avg.bray <- laply(ISS.iters.avg.bray, as.matrix)
ISS.bray.avg <- as.dist(apply(ISS.array.avg.bray, 2:3, mean))
rm(ISS.iters.avg.bray, ISS.array.avg.bray)
#write.table(ISS.bray.avg %>% as.matrix() %>% as.data.frame(), file ="data/ISS.bray.avg.txt")
write.table(ISS.bray.avg %>% as.matrix() %>% as.data.frame(), file ="data/withK19/ISS.bray.avg.txt")

# Jaccard
ISS.iters.avg.jac <- mclapply(ISS.iters, function(x, method ="jaccard") avg.vegdist(x), mc.cores = 1)
ISS.array.avg.jac <- laply(ISS.iters.avg.jac, as.matrix)
ISS.jac.avg <- as.dist(apply(ISS.array.avg.jac, 2:3, mean))
rm(ISS.iters.avg.jac, ISS.array.avg.jac)
#write.table(ISS.jac.avg %>% as.matrix() %>% as.data.frame(), file ="data/ISS.jac.avg.txt")
write.table(ISS.jac.avg %>% as.matrix() %>% as.data.frame(), file ="data/withK19/ISS.jac.avg.txt")
#data  <- read.table("ISS.bray.txt") %>% as.matrix() %>% as.dist() # How to read in the data again

# Plotting beta diversity with NTC
asv.NTC 
asv.NTC.bray <-vegdist(t(asv.NTC), method ="bray")
nmds.BLW2 <- metaMDS(asv.NTC.bray)
nmds.BLW2.points = nmds.BLW2$points %>% as.data.frame()
colnames(nmds.BLW2.points) = c("nmds1", "nmds2")
nmds.BLW2.points$treatment = rownames(nmds.BLW2.points)
nmds.BLW2.points$treatment <- sub("_[1-9]", "", nmds.BLW2.points$treatment)
nmds.BLW2$stress # 0.097

ggplot(data = nmds.BLW2.points, aes(x = nmds1, y = nmds2))+
  geom_point(aes(color =treatment))+
  theme_classic(base_size = 14) + theme(legend.position = "right") +
  geom_text(x=5,y=-2,label='stress = 0.097', size=4)+
  ylab("NMDS2")+xlab("NMDS1")+
  coord_fixed(ratio =1)+
  geom_text_repel(aes(x = nmds1 , y = nmds2 , label = treatment), show.legend = FALSE, max.overlaps =20)


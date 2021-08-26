# Identify potential pollen parents for hybrid individuals using RNAseq data mapped to Alloteropsis semialata genome (CDS)
rm(list = ls())
library(data.table)

setwd("/path/to/directory") 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Format table containing only multiallelic sites for which all hybrids are heterozygous and all nonhybrids are homozygous
# 1.1. Read raw table
tb <- fread("variants_only_multiallelic_sites_heter_and_homoz_reformatted.txt", h = T)
head(tb)

# 1.2. Formatting
# Create a unique site name
tb$site <- paste(tb$CHROM, tb$POS, sep="-")

# Rearrange columns
tb$CHROM = tb$site
tb <- tb[,-c(2,3,4)]

# Transpose
x <- dcast(melt(tb, id.vars = "site"), variable ~ site)
colnames(x)[1] <- "sample"
tb <- x
tb <- tb[-1,]
tb[1:20, 1:6]

# Get metadata
meta <- read.csv("metadata.csv", h=T)
head(meta)

# Add photosynthetic type
photo <-  c()
for(i in 1:length(as.character(tb$sample))){
  photo <- c(photo, as.character(subset(meta$photo, meta$accession == as.character(tb$sample)[i])))
}
tb$photo <- photo

# Define hybrids
key1 <- c()
for(i in 1:length(as.character(tb$sample))){
  if(grepl("hybrid", tb$sample[i])) {key1 <- c(key1, "hybrid")} else {key1 <- c(key1, "non_hybrid")}
}
tb$key <- factor(tb$key, levels=unique(c("non_hybrid", "hybrid")))
tb[1:20, 1:6]

# Rearrange columns again and sort 
df <- data.frame(accession = as.character(tb$sample), key = key1, photo = photo)
df[3,2] = "hybrid" # EML11-2 is a hybrid 
df <- data.frame(cbind(df, tb[,2:dim(x)[2]]))
df[1:20, 1:7]
df$photo <- factor(df$photo, levels=unique(c("C3", "C3 x C4", "C4", "C3-C4 x C4", "C3-C4")))
df <- df[with(df, order(df$key, df$photo)),]
df[1:33, 1:6]

# Save it
write.table(df, "variants_only_multiallelic_sites_heter_and_homoz_reformatted_sorted.txt")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 2. Determine the type of cross (c3 x c4, c3-c4 x c4)

## 2.1. Identify markers that can differentiate photosythetic types 
good_int_mark <- c()
good_c3_mark <- c()
good_c4_mark <- c()
good_marker <- c()
options(warn=2) # interrupt loop when there is a warning 

# Loop through markers, save markers that can differentiate photosythetic types in the table "good_marker" 
for(i in 4:dim(df)[2]){
  t <- df[,c(1,2,3,i)]
  print(dim(df)[2]-i)

  # 1) Check if all accession within a photo type have the same genotype
    # 1.1. Check if there is less than 25% missing data; If so, 
    # 1.2. Check if there is more than one variant in this photo type (after excluding missing data), if not, save it 
  
  checkbases.c3 <- c()
  checkbases.int <- c()
  checkbases.c4 <- c()
  
  c3.out = F
  int.out = F
  c4.out = F
  
  # C3
  if(length(grep("-/-", subset(t, t$photo == "C3")[,4]))/length(as.character(subset(t, t$photo == "C3")[,4])) <= 0.25) {
    checkbases.c3 <- unique(unlist(strsplit(subset(t, t$photo == "C3")[,4], "/")))
    if(length(checkbases.c3[!grepl("-", checkbases.c3)]) == 1) {good_c3_mark <- c(good_c3_mark, colnames(t)[4])} else {good_c3_mark=good_c3_mark;
    c3.out = T}
  } else {good_c3_mark=good_c3_mark; c3.out = T}
  
  # C3-C4
  if(length(grep("-/-", subset(t, t$photo == "C3-C4")[,4]))/length(as.character(subset(t, t$photo == "C3-C4")[,4])) <= 0.25) {
    checkbases.int <- unique(unlist(strsplit(subset(t, t$photo == "C3-C4")[,4], "/")))
    if(length(checkbases.int[!grepl("-", checkbases.int)]) == 1) {good_int_mark <- c(good_int_mark, colnames(t)[4])} else {good_int_mark=good_int_mark;
    int.out = T}
  } else {good_int_mark=good_int_mark; int.out = T}
  
  # C4 
  if(length(grep("-/-", subset(t, t$photo == "C4")[,4]))/length(as.character(subset(t, t$photo == "C4")[,4])) <= 0.25) {
    checkbases.c4 <- unique(unlist(strsplit(subset(t, t$photo == "C4")[,4], "/")))
    if(length(checkbases.c4[!grepl("-", checkbases.c4)]) == 1) {good_c4_mark <- c(good_c4_mark, colnames(t)[4])} else {good_c4_mark=good_c4_mark;
    c4.out = T}
  } else {good_c4_mark=good_c4_mark; c4.out = T}
  
  # 2. Check if there are multiple alleles within each photo type first, then 
  # 2.1. If c3, int and c4 have different bases, save the marker in a separate vector
  if(any(c3.out, int.out, c4.out)) {good_marker = good_marker} else {
    if(checkbases.c3[!grepl("-", checkbases.c3)] != checkbases.int[!grepl("-", checkbases.int)] & 
       checkbases.c3[!grepl("-", checkbases.c3)] != checkbases.c4[!grepl("-", checkbases.c4)] &
       checkbases.int[!grepl("-", checkbases.int)] != checkbases.c4[!grepl("-", checkbases.c4)]) {good_marker <- c(good_marker, colnames(t)[4])} else {
         good_marker = good_marker}
  }
}
 
good_marker
length(good_marker)

dim(subset(df, select = good_marker))
gm <- data.frame(cbind(df[,c(1,2,3)], subset(df, select = good_marker)))
gm[1:10, 1:5]

# Save it 
write.csv(gm, "table_selected_markers_diff_phototypes.csv")

## 2.1. Calculate the proportion of times the hybrids match with the expected genotypes
# Import it
#gm <- read.csv("table_selected_markers_diff_phototypes.csv")
head(gm)
gm[1:10, 1:5]

gm$accession <- as.character(gm$accession)

# Go through each hybrid accession and check whether the genotype matches the expected genotype based on the cross type
total.valid.snps <- c()
total.positive.matches <- c()
wrong.matches <- c()
for(i in 1:length(as.character(subset(gm$accession, gm$key == "hybrid")))){
  # create vector to save the total number of snps being analysed (to be able to discard missing data)
  total.snps = 0 
  # create vector to save the total number of positive matches (i.e. genotype is the one expected for that cross)
  positive.matches = 0
  # get the type of cross 
  cross.type <- as.character(subset(gm$photo, gm$accession == as.character(subset(gm$accession, gm$key == "hybrid"))[i]))
  # get the accession name
  acc <- as.character(subset(gm$accession, gm$key == "hybrid"))[i]
  # save names of SNPs not matching
  nomatch <- c()
  
  # go through each SNP and check whether it matches the expected or not 
  for(j in 4:dim(gm)[2]){
    print(paste(acc, dim(gm)[2]-j))
    # get the SNP
    t <- gm[,c(1,2,3,j)]
    t[,4] <- as.character(t[,4])
    # get the bases for each photo type (these are all homozygous and different from each other; see previous loop)
    checkbases.c3 <- unique(unlist(strsplit(subset(t, t$photo == "C3")[,4], split = "/")))
    checkbases.c3 = checkbases.c3[!grepl("-", checkbases.c3)]
    checkbases.int <- unique(unlist(strsplit(subset(t, t$photo == "C3-C4")[,4], "/")))
    checkbases.int = checkbases.int[!grepl("-", checkbases.int)]
    checkbases.c4 <- unique(unlist(strsplit(subset(t, t$photo == "C4")[,4], "/")))  
    checkbases.c4 = checkbases.c4[!grepl("-", checkbases.c4)]
    # get the genotype of the accession i in the SNP j
    gt1 <- unlist(strsplit(subset(t[,4], t$accession == acc), "/"))
    # If it's missing data, skip it, else test whether it matches or not
    if("-" %in% gt1) {total.snps = total.snps} else {total.snps = total.snps+1;
      if(cross.type == "C3 x C4"){
        if(checkbases.c3 %in% gt1 & checkbases.c4 %in% gt1) {positive.matches = positive.matches+1} else {positive.matches=positive.matches; 
        nomatch <- c(nomatch, colnames(t)[4])}
       } else if(cross.type == "C3-C4 x C4"){
          if(checkbases.int %in% gt1 & checkbases.c4 %in% gt1) {positive.matches = positive.matches+1} else {positive.matches=positive.matches;
          nomatch <- c(nomatch, colnames(t)[4])}
          }
    }
  }
  # save total snps and positive matches
  total.valid.snps <- c(total.valid.snps, total.snps)
  total.positive.matches <- c(total.positive.matches, positive.matches)
  temp1 <- data.frame(acc = rep(acc, length(nomatch)), cross = rep(cross.type, length(nomatch)), nomatch = nomatch)
  wrong.matches <- data.frame(rbind(wrong.matches, temp1))
}
result.genotype.phototype <- data.frame(acc = as.character(subset(gm$accession, gm$key == "hybrid")), 
                                        photo = as.character(subset(gm$photo, gm$key == "hybrid")),
                                        total.valid.snps = total.valid.snps, total.positive.matches = total.positive.matches,
                                        prop = total.positive.matches/total.valid.snps)
result.genotype.phototype
wrong.matches

# Save results
write.csv(result.genotype.phototype, "table_selected_markers_diff_phototypes_correct_genotype_matches.csv")
write.csv(wrong.matches, "table_selected_markers_diff_phototypes_correct_genotype_matches_list_markers.csv")

levels(wrong.matches$nomatch)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 3. Determine the potential pollen parents - Part1: find snps that are unique within each photo type ("singletons")

# Keeping only one accession per pop (prevents reducing the number of singletons per potential pollen parent, as closely related accessions - for example, from the same population - will likely have the same genotype) 
df2 <- df
nonhyb
df2 <- df2[!grepl("KWT3_B_1_14_merged", df2$accession),]
df2 <- df2[!grepl("AUS2_B_1_14_merged", df2$accession),]
df2 <- df2[!grepl("parent_A_TWN2-1", df2$accession),]
df2 <- df2[!grepl("PHIL-17_CAGATC_L001", df2$accession),]
df2 <- df2[!grepl("TW10_B_1_14_merged", df2$accession),]
df2 <- df2[!grepl("TAN02-1_ACTTGA_L001", df2$accession),]
df2$accession <- factor(df2$accession) 

df2[1:28, 1:6]
dim(df2)
nonhyb <- as.character(subset(df2$accession, df2$key == "non_hybrid"))
unique.markers <- data.frame()

# Loop through markers 
for(i in 4:dim(df2)[2]){
  # extract marker
  t <- df2[,c(1,2,3,i)]
  print(dim(df2)[2]-i)
  # go through each accession 
  for(k in 1:length(nonhyb)){
    # get photo type
    pt <- as.character(subset(t$photo, t$accession == nonhyb[k]))
    # get all genotypes for that photo type
    gt.all <- subset(t[,4], t$photo == pt)
    # get all genotypes for the other parent 
    if(pt == "C3") {gt.comp <- subset(t[,4], t$photo == "C4")} else if(pt == "C3-C4") {gt.comp <- subset(t[,4], t$photo == "C4")} else if(pt == "C4"){
      gt.comp <- subset(t[,4], t$photo == "C3" | t$photo == "C3-C4")
    }
    
    # get genotype for the accession
    gt.acc <- subset(t[,4], t$accession == nonhyb[k])
    # save the accession and the marker
    log.mark <- data.frame(accession = nonhyb[k], photo = pt, marker = colnames(t)[4], genotype = gt.acc)
    # save the marker if the genotype for accession k is unique within the accessions from the same photo type, and there is no accession in the other 
    #   potential parent with the same genotype
    # (skip if missing data)
    if(gt.acc == "-/-") {unique.markers=unique.markers} else {
      if(length(grep(gt.acc, gt.all)) == 1 & gt.acc %in% gt.comp == F) {unique.markers <- data.frame(rbind(unique.markers, log.mark))} else {
        unique.markers=unique.markers}
    }
  }
}

dim(unique.markers)
tail(unique.markers)
back1 <- unique.markers
head(back1)
unique.markers <- unique.markers[with(unique.markers, order(unique.markers$accession)),]
head(unique.markers)
tail(unique.markers)

# save
write.csv(unique.markers, "table_unique_markers_per_accession_representative_set.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 4. Determine the potential pollen parents - Part2: find positive matches between the hybrid allele and the sample-specific singletons

# Read list of unique markers
uma <- read.csv("table_unique_markers_per_accession_representative_set.csv")
uma <- uma[,-1]
head(uma)
head(subset(uma, uma$accession == "parent_B_Aus1-1"))

# create new df2 with only unique markers
cbind(df2[,c(1,2,3)], df2$ASEM_AUS1_26507.158)
df2uni <- data.frame(cbind(df2[,c(1,2,3)], subset(df2, select = levels(uma$marker))))
df2uni[1:30, 1:6]

length(unique(as.character(subset(uma$marker, uma$photo == "C4"))))

# 4.1. C3 x C4 hybrids (for cases where C3 mom is known)

# for each marker where there is at least one c4 accession with a unique genotype, 
potential.parent.C3xC4 <- data.frame()
options(warn=2) # interrupt loop when there is a warning 
tomarkers <- length(unique(as.character(subset(uma$marker, uma$photo == "C4"))))
for(i in unique(as.character(subset(uma$marker, uma$photo == "C4")))){
  # i = "ASEM_AUS1_18006.2908"
  print(paste(i, tomarkers - grep(i, unique(as.character(subset(uma$marker, uma$photo == "C4")))), sep="       "))
  # identify the "potential father" to be tested, 
  #   i.e. the C4 accession for which marker "i" has a unique homozygous genotype (i.e. there isn't any other C4 with the same genotype)
  acc.uni <- sort(as.character(subset(uma$accession, uma$marker == i & uma$photo == "C4"))) # sort to put it in the same order as in tab1 (below)
  
  # reduce the main dataset to the marker i
  tab1 <- cbind(df2[,c(1,2,3)], gt = as.character(df2[,which(colnames(df2)==i)]))
  
  # get the C3 genotype
  c3.genot <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$photo == "C3")), "/")))
  c3.genot = c3.genot[!grepl("-", c3.genot)]
  
  # T1. check if the c3 genotype is available - if not, or if there are multiple c3, skip it
  if(length(c3.genot) != 1){potential.parent.C3xC4=potential.parent.C3xC4} else{
    
    # get the genotype of the "potential father"
    pot.father.gt <- as.character(tab1[tab1$accession %in% acc.uni,]$gt)
    
    # proportion of missing genotypes within the C4 - gives an idea of how likely this is the actual parent 
    # (i.e. if many c4s are missing, than it is possible that there are other c4s with the same genotype)
    prop.c4.missing <- round(length(grep("-/-", as.character(subset(tab1$gt, tab1$photo == "C4"))))/
                               length(as.character(subset(tab1$gt, tab1$photo == "C4"))),2)
    
    # L2. Go through each hybrid and test if the allele of the potential parent matches with its c4 allele
    for(k in as.character(subset(tab1$accession, tab1$photo == "C3 x C4"))){
      # k = "hybrid_H05_EML11-200_E" 
      # get the hybrid genotype
      gt.test <- unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == k)), "/"))
      
      # T2. check if the c3 genotype is in the hybrid, if not skip it (potentially allele-specific expression in the c3 parents)
      if(c3.genot %in% gt.test == F) {potential.parent.C3xC4=potential.parent.C3xC4} else {
        
        # T3. test if there are multiple c4 potential parents, if so, check one by one to see if its genotype matches with the hybrid
        if(length(pot.father.gt) > 1){
          for(p in 1:length(pot.father.gt)){
            genot1 <- unique(unlist(strsplit(pot.father.gt[p], "/")));
            
            # prepare log, if to be saved
            dad.log <- data.frame(hybrid = k, marker = i,
                                  genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                  c3_allele = c3.genot, potential_father = acc.uni[p], 
                                  potential_father_allele = genot1, c4.missing = prop.c4.missing);
            # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
            if("-" %in% gt.test) {potential.parent.C3xC4=potential.parent.C3xC4} else {
              actual.c4.gt <- gt.test[!grepl(c3.genot, gt.test)]
              if(genot1 == actual.c4.gt) {potential.parent.C3xC4 <- data.frame(rbind(potential.parent.C3xC4, dad.log))}
            }
          }
        } else {
          pot.father.gt <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == acc.uni)), "/")));
          # prepare log, if to be saved
          dad.log <- data.frame(hybrid = k, marker = i,
                                genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                c3_allele = c3.genot, potential_father = acc.uni, potential_father_allele = pot.father.gt, 
                                c4.missing = prop.c4.missing);
          # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
          if("-" %in% gt.test) {potential.parent.C3xC4=potential.parent.C3xC4} else {
            actual.c4.gt <- gt.test[!grepl(c3.genot, gt.test)]
            if(pot.father.gt == actual.c4.gt) {potential.parent.C3xC4 <- data.frame(rbind(potential.parent.C3xC4, dad.log))}
          }
        } # T3 test if there are multiple c4 potential parents
      } # T2 test if the c3 allele is in the hybrid
    } # L1 loop through hybrids
  } # T1 test if c3 genotype is available
} # END
potential.parent.C3xC4

write.csv(potential.parent.C3xC4, "table_potential_C4_parent_for_c3xc4_hybrids_representative_accessions.csv")

head(potential.parent.C3xC4)

c3xc4dad.filtered <- subset(potential.parent.C3xC4, potential.parent.C3xC4$c4.missing < 0.5)
c3xc4dad.filtered.summary <- table(c3xc4dad.filtered[,c(1,5)])
write.csv(c3xc4dad.filtered.summary, "table_potential_C4_parent_for_c3xc4_hybrids_summary_filter50perc_representative_accessions.csv")

# # # # 
# # # # 
# 4.2 C3-C4 x C4 hybrids ((for cases where the C3-C4 mom is known)

# for each marker where there is at least one c4 accession with a unique genotype, 
potential.parent.c3_c4xc4 <- data.frame()
options(warn=2) # interrupt loop when there is a warning 
tomarkers <- length(unique(as.character(subset(uma$marker, uma$photo == "C4"))))
for(i in unique(as.character(subset(uma$marker, uma$photo == "C4")))){
  #i = "ASEM_AUS1_18006.2908"
  print(paste(i, tomarkers - grep(i, unique(as.character(subset(uma$marker, uma$photo == "C4")))), sep="       "))
  # identify the "potential father" to be tested, 
  #   i.e. the C4 accession for which marker "i" has a unique homozygous genotype (i.e. there isn't any other C4 with the same genotype)
  acc.uni <- sort(as.character(subset(uma$accession, uma$marker == i & uma$photo == "C4"))) # sort to put it in the same order as in tab1 (below)
  
  # reduce the main dataset to the marker i
  tab1 <- cbind(df2[,c(1,2,3)], gt = as.character(df2[,which(colnames(df2)==i)]))
  
  # get the C3 genotype
  c3_c4.genot <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$photo == "C3-C4")), "/")))
  c3_c4.genot = c3_c4.genot[!grepl("-", c3_c4.genot)]
  
  # T1. check if the c3 genotype is available - if not, or if there are multiple c3, skip it
  if(length(c3_c4.genot) != 1){potential.parent.c3_c4xc4=potential.parent.c3_c4xc4} else{
    
    # get the genotype of the "potential father"
    pot.father.gt <- as.character(tab1[tab1$accession %in% acc.uni,]$gt)
    
    # proportion of missing genotypes within the C4 - gives an idea of how likely this is the actual parent 
    # (i.e. if many c4s are missing, than it is possible that there are other c4s with the same genotype)
    prop.c4.missing <- round(length(grep("-/-", as.character(subset(tab1$gt, tab1$photo == "C4"))))/
                               length(as.character(subset(tab1$gt, tab1$photo == "C4"))),2)
    
    # L2. Go through each hybrid and test if the allele of the potential parent matches with its c4 allele
    for(k in as.character(subset(tab1$accession, tab1$photo == "C3-C4 x C4"))){
      # k = "hybrid_H07_TAN16-2-3C_X3" 
      # get the hybrid genotype
      gt.test <- unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == k)), "/"))
      
      # T2. check if the c3 genotype is in the hybrid, if not skip it (potentially allele-specific expression in the c3 parents)
      if(c3_c4.genot %in% gt.test == F) {potential.parent.c3_c4xc4=potential.parent.c3_c4xc4} else {
        
        # T3. test if there are multiple c4 potential parents, if so, check one by one to see if its genotype matches with the hybrid
        if(length(pot.father.gt) > 1){
          for(p in 1:length(pot.father.gt)){
            genot1 <- unique(unlist(strsplit(pot.father.gt[p], "/")));
            
            # prepare log, if to be saved
            dad.log <- data.frame(hybrid = k, marker = i,
                                  genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                  c3_allele = c3_c4.genot, potential_father = acc.uni[p], 
                                  potential_father_allele = genot1, c4.missing = prop.c4.missing);
            # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
            if("-" %in% gt.test) {potential.parent.c3_c4xc4=potential.parent.c3_c4xc4} else {
              actual.c4.gt <- gt.test[!grepl(c3_c4.genot, gt.test)]
              if(genot1 == actual.c4.gt) {potential.parent.c3_c4xc4 <- data.frame(rbind(potential.parent.c3_c4xc4, dad.log))}
            }
          }
        } else {
          pot.father.gt <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == acc.uni)), "/")));
          # prepare log, if to be saved
          dad.log <- data.frame(hybrid = k, marker = i,
                                genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                c3_allele = c3_c4.genot, potential_father = acc.uni, potential_father_allele = pot.father.gt, 
                                c4.missing = prop.c4.missing);
          # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
          if("-" %in% gt.test) {potential.parent.c3_c4xc4=potential.parent.c3_c4xc4} else {
            actual.c4.gt <- gt.test[!grepl(c3_c4.genot, gt.test)]
            if(pot.father.gt == actual.c4.gt) {potential.parent.c3_c4xc4 <- data.frame(rbind(potential.parent.c3_c4xc4, dad.log))}
          }
        } # T3 test if there are multiple c4 potential parents
      } # T2 test if the c3 allele is in the hybrid
    } # L1 loop through hybrids
  } # T1 test if c3 genotype is available
} # END
potential.parent.c3_c4xc4

write.csv(potential.parent.c3_c4xc4, "table_potential_C4_parent_for_c3_c4xc4_hybrids_representative_accessions.csv")

head(potential.parent.c3_c4xc4)

c3_c4xc4dad.filtered <- subset(potential.parent.c3_c4xc4, potential.parent.c3_c4xc4$c4.missing < 0.5)
c3_c4xc4dad.filtered.summary <- table(c3_c4xc4dad.filtered[,c(1,5)])
write.csv(c3_c4xc4dad.filtered.summary, "table_potential_C4_parent_for_c3_c4xc4_hybrids_summary_filter50perc_representative_accessions.csv")

# # # # 
# # # # 
# 4.2. C3 x C4 hybrids (for cases where the C4 mom is known)

# for each marker where there is at least one c4 accession with a unique genotype, 
potential.parent.c4xc3 <- data.frame()
options(warn=2) # interrupt loop when there is a warning 
tomarkers <- length(unique(as.character(subset(uma$marker, uma$photo == "C3"))))
for(i in unique(as.character(subset(uma$marker, uma$photo == "C3")))){
  #i = "ASEM_AUS1_32856.378"
  print(paste(i, tomarkers - grep(i, unique(as.character(subset(uma$marker, uma$photo == "C3")))), sep="       "))
  # identify the "potential father" to be tested, 
  #   i.e. the C3 accession for which marker "i" has a unique homozygous genotype (i.e. there isn't any other C3 with the same genotype)
  acc.uni <- sort(as.character(subset(uma$accession, uma$marker == i & uma$photo == "C3"))) # sort to put it in the same order as in tab1 (below)
  
  # reduce the main dataset to the marker i
  tab1 <- cbind(df2[,c(1,2,3)], gt = as.character(df2[,which(colnames(df2)==i)]))
  
  # get the C4 genotype
  c4.genot <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$photo == "C4")), "/")))
  c4.genot = c4.genot[!grepl("-", c4.genot)]
  
  # T1. check if the c3 genotype is available - if not, or if there are multiple c3, skip it
  if(length(c4.genot) != 1){potential.parent.c4xc3=potential.parent.c4xc3} else{
    
    # get the genotype of the "potential father"
    pot.father.gt <- as.character(tab1[tab1$accession %in% acc.uni,]$gt)
    
    # proportion of missing genotypes within the C4 - gives an idea of how likely this is the actual parent 
    # (i.e. if many c4s are missing, than it is possible that there are other c4s with the same genotype)
    prop.c3.missing <- round(length(grep("-/-", as.character(subset(tab1$gt, tab1$photo == "C3"))))/
                               length(as.character(subset(tab1$gt, tab1$photo == "C3"))),2)
    
    # L2. Go through each hybrid and test if the allele of the potential parent matches with its c4 allele
    for(k in c("hybrid_H20_TWN11-1_X1")){
      #k = "hybrid_H20_TWN11-1_X1"
      # get the hybrid genotype
      gt.test <- unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == k)), "/"))
      
      # T2. check if the c4 genotype is in the hybrid, if not skip it (potentially allele-specific expression in the c4 parents)
      if(c4.genot %in% gt.test == F) {potential.parent.c4xc3=potential.parent.c4xc3} else {
        
        # T3. test if there are multiple c3 potential parents, if so, check one by one to see if its genotype matches with the hybrid
        if(length(pot.father.gt) > 1){
          for(p in 1:length(pot.father.gt)){
            genot1 <- unique(unlist(strsplit(pot.father.gt[p], "/")));
            
            # prepare log, if to be saved
            dad.log <- data.frame(hybrid = k, marker = i,
                                  genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                  c4_allele = c4.genot, potential_father = acc.uni[p], 
                                  potential_father_allele = genot1, c3.missing = prop.c3.missing);
            # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
            if("-" %in% gt.test) {potential.parent.c4xc3=potential.parent.c4xc3} else {
              actual.c3.gt <- gt.test[!grepl(c4.genot, gt.test)]
              if(genot1 == actual.c3.gt) {potential.parent.c4xc3 <- data.frame(rbind(potential.parent.c4xc3, dad.log))}
            }
          }
        } else {
          pot.father.gt <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == acc.uni)), "/")));
          # prepare log, if to be saved
          dad.log <- data.frame(hybrid = k, marker = i,
                                genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                c4_allele = c4.genot, potential_father = acc.uni, potential_father_allele = pot.father.gt, 
                                c3.missing = prop.c3.missing);
          # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
          if("-" %in% gt.test) {potential.parent.c4xc3=potential.parent.c4xc3} else {
            actual.c3.gt <- gt.test[!grepl(c4.genot, gt.test)]
            if(pot.father.gt == actual.c3.gt) {potential.parent.c4xc3 <- data.frame(rbind(potential.parent.c4xc3, dad.log))}
          }
        } # T3 test if there are multiple c4 potential parents
      } # T2 test if the c3 allele is in the hybrid
    } # L1 loop through hybrids
  } # T1 test if c3 genotype is available
} # END
potential.parent.c4xc3

write.csv(potential.parent.c4xc3, "table_potential_C3_parent_for_c4xc3_hybrids_representative_accessions.csv")

head(potential.parent.c4xc3)

c4xc3dad.filtered <- subset(potential.parent.c4xc3, potential.parent.c4xc3$c3.missing < 0.5)
c4xc3dad.filtered.summary <- table(c4xc3dad.filtered[,c(1,5)])
write.csv(c4xc3dad.filtered.summary, "table_potential_C3_parent_for_c4xc3_hybrids_summary_filter50perc_representative_accessions.csv")

# # # # 
# # # # 
# 4.4. C3-C4 x C4 hybrids (for cases where the C4 mom is known)

# for each marker where there is at least one c4 accession with a unique genotype, 
potential.parent.c4xc3_c4 <- data.frame()
options(warn=2) # interrupt loop when there is a warning 
tomarkers <- length(unique(as.character(subset(uma$marker, uma$photo == "C3-C4"))))
for(i in unique(as.character(subset(uma$marker, uma$photo == "C3-C4")))){
  #i = "ASEM_AUS1_26507.158"
  print(paste(i, tomarkers - grep(i, unique(as.character(subset(uma$marker, uma$photo == "C3-C4")))), sep="       "))
  # identify the "potential father" to be tested, 
  #   i.e. the C3-C4 accession for which marker "i" has a unique homozygous genotype (i.e. there isn't any other C3-C4 with the same genotype)
  acc.uni <- sort(as.character(subset(uma$accession, uma$marker == i & uma$photo == "C3-C4"))) # sort to put it in the same order as in tab1 (below)
  
  # reduce the main dataset to the marker i
  tab1 <- cbind(df2[,c(1,2,3)], gt = as.character(df2[,which(colnames(df2)==i)]))
  
  # get the C4 genotype
  c4.genot <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$photo == "C4")), "/")))
  c4.genot = c4.genot[!grepl("-", c4.genot)]
  
  # T1. check if the c3 genotype is available - if not, or if there are multiple c3, skip it
  if(length(c4.genot) != 1){potential.parent.c4xc3_c4=potential.parent.c4xc3_c4} else{
    
    # get the genotype of the "potential father"
    pot.father.gt <- as.character(tab1[tab1$accession %in% acc.uni,]$gt)
    
    # proportion of missing genotypes within the C4 - gives an idea of how likely this is the actual parent 
    # (i.e. if many c4s are missing, than it is possible that there are other c4s with the same genotype)
    prop.c3_c4.missing <- round(length(grep("-/-", as.character(subset(tab1$gt, tab1$photo == "C3-C4"))))/
                                  length(as.character(subset(tab1$gt, tab1$photo == "C3-C4"))),2)
    
    # L2. Go through each hybrid and test if the allele of the potential parent matches with its c4 allele
    for(k in c("hybrid_H18_TWN10-2_X7")){
      #k = "hybrid_H18_TWN10-2_X7"
      # get the hybrid genotype
      gt.test <- unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == k)), "/"))
      
      # T2. check if the c4 genotype is in the hybrid, if not skip it (potentially allele-specific expression in the c4 parents)
      if(c4.genot %in% gt.test == F) {potential.parent.c4xc3_c4=potential.parent.c4xc3_c4} else {
        
        # T3. test if there are multiple c4 potential parents, if so, check one by one to see if its genotype matches with the hybrid
        if(length(pot.father.gt) > 1){
          for(p in 1:length(pot.father.gt)){
            genot1 <- unique(unlist(strsplit(pot.father.gt[p], "/")));
            
            # prepare log, if to be saved
            dad.log <- data.frame(hybrid = k, marker = i,
                                  genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                  c4_allele = c4.genot, potential_father = acc.uni[p], 
                                  potential_father_allele = genot1, c3_c4.missing = prop.c3_c4.missing);
            # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
            if("-" %in% gt.test) {potential.parent.c4xc3_c4=potential.parent.c4xc3_c4} else {
              actual.c3_c4.gt <- gt.test[!grepl(c4.genot, gt.test)]
              if(genot1 == actual.c3_c4.gt) {potential.parent.c4xc3_c4 <- data.frame(rbind(potential.parent.c4xc3_c4, dad.log))}
            }
          }
        } else {
          pot.father.gt <- unique(unlist(strsplit(as.character(subset(tab1$gt, tab1$accession == acc.uni)), "/")));
          # prepare log, if to be saved
          dad.log <- data.frame(hybrid = k, marker = i,
                                genotype = as.character(subset(tab1$gt, tab1$accession == k)), 
                                c4_allele = c4.genot, potential_father = acc.uni, potential_father_allele = pot.father.gt, 
                                c3_c4.missing = prop.c3_c4.missing);
          # 1) skip missing data, 2) remove genotype from the c3 parent, 3) test whether the actual c4 genotype matches 
          if("-" %in% gt.test) {potential.parent.c4xc3_c4=potential.parent.c4xc3_c4} else {
            actual.c3_c4.gt <- gt.test[!grepl(c4.genot, gt.test)]
            if(pot.father.gt == actual.c3_c4.gt) {potential.parent.c4xc3_c4 <- data.frame(rbind(potential.parent.c4xc3_c4, dad.log))}
          }
        } # T3 test if there are multiple c4 potential parents
      } # T2 test if the c3 allele is in the hybrid
    } # L1 loop through hybrids
  } # T1 test if c3 genotype is available
} # END
potential.parent.c4xc3_c4

write.csv(potential.parent.c4xc3_c4, "table_potential_C3-C4_parent_for_c4xc3_c4_hybrids_representative_accessions.csv")

head(potential.parent.c4xc3_c4)

c4xc3_c4dad.filtered <- subset(potential.parent.c4xc3_c4, potential.parent.c4xc3_c4$c3_c4.missing < 0.5)
c4xc3_c4dad.filtered.summary <- table(c4xc3_c4dad.filtered[,c(1,5)])
write.csv(c4xc3_c4dad.filtered.summary, "table_potential_C3-C4_parent_for_c4xc3_c4_hybrids_summary_filter50perc_representative_accessions.csv")


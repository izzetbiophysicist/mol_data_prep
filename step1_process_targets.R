#awk -F '\t' '{print $2 "\t" $10 "\t" $7}' BindingDB_All_202312.tsv > BindingDB_filter.tsv
library(readr)
#inp <- read.csv("~/virtual_screening_pipeline/tst", sep='\t')

inp <- read.csv("~/virtual_screening_pipeline/BindingDB_filter.tsv", sep='\t')
inp <- inp[which(inp$IC50..nM. != '' & is.na(inp$IC50..nM.) == FALSE),]

inp$targ_org <- paste(inp$Target.Name, inp$Target.Source.Organism.According.to.Curator.or.DataSource)

inp$IC50..nM. <- as.numeric(inp$IC50..nM.)
inp <- na.omit(inp)

###############
### Clean dataframes
###############


##############
uniq_targets <- data.frame(unique_targets=unique(sort(inp$targ_org)), n=NA)
for(i in 1:nrow(uniq_targets)){
  uniq_targets$n[i] <- length(which(inp$targ_org == uniq_targets$unique_targets[i]))
  #uniq_targets$n[i] <- length(which(inp$UniProt..SwissProt..Primary.ID.of.Target.Chain == uniq_targets$unique_targets[i]))
}

uniq_filter <- uniq_targets[which(uniq_targets$n > 2000),]

######
uniq_filter$directory <- NA
for(i in 1:nrow(uniq_filter)){
  path <- my_string_cleaned <- gsub("[^[:alnum:]]", "", gsub(" ", "", uniq_filter$unique_targets[i]))
  uniq_filter$directory[i] <- path
  dir.create(paste("/home/lucas/virtual_screening_pipeline/datasets/", path, sep=""))
  tmp <- inp[which(inp$targ_org == uniq_filter$unique_targets[i]),]
  tmp2 <- data.frame(smiles=tmp$Ligand.SMILES, IC50=tmp$IC50..nM.)
  write.csv(tmp2, paste("/home/lucas/virtual_screening_pipeline/datasets/", path, '/dataset.csv',sep=""), quote = FALSE, row.names = FALSE)
}

write.csv(uniq_filter,'/home/lucas/virtual_screening_pipeline/directory_map.csv')

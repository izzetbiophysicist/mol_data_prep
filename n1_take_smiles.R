###########
### Process smiles
###########
library(readr)
library(rcdk)

approved <- read.csv("~/ML_integrase_keras/all_over/structure_links.csv")
app_index <- 1:nrow(approved)
bdb <- read.csv("~/ML_integrase_keras/all_over/BindingDB_All.tsv", sep="\t",  comment.char = "")
decoys <- read.csv("~/ML_integrase_keras/all_over/decoy.smiles")

###################
### filtrar BDB
###################
hivin <- c("HIV-1 Integrase", "HIV-1 integrase (IN) A128T", "Integrase", "Human immunodeficiency virus type 1 integrase")


bdb_hiv <- bdb[which(as.character(bdb$Target.Name.Assigned.by.Curator.or.DataSource) %in% hivin),]


#a <-  data.frame(unique(sort(as.character(bdb$Target.Name.Assigned.by.Curator.or.DataSource))))
#length(which(a$unique.sort.as.character.bdb.Target.Name.Assigned.by.Curator.or.DataSource... %in% hivin[3]))
#write.csv(a, '~/ML_integrase_keras/all_over/a.txt', quote = FALSE, row.names = FALSE)
#bdb[which(bdb$BindingDB.Reactant_set_id == "BDBM50495689"),]

#######################
### Filtrar ativos/inativos
#######################

#####
# https://aac.asm.org/content/aac/41/2/385.full.pdf   cutoff < 100 ug
## https://onlinelibrary.wiley.com.sci-hub.se/doi/epdf/10.1002/ptr.2252 same for protease (<100 ug)
### https://www.researchgate.net/post/What_is_the_criteria_of_selecting_molecules_as_actives_and_inactives_based_on_IC50_value
#####

#### Remover os sem informacao
bdb_hiv <- bdb_hiv[which(bdb_hiv$IC50..nM. != ""), ]


####### escrever sinal
bdb_hiv$sign <- "-"
for(i in 1:nrow(bdb_hiv)){
  
  bla <- strsplit(bdb_hiv[i,10], split="")[[1]]
  if(bla[1] == " "){
  bla <- bla[-which(bla == " ")]
  }
  if(bla[1] == ">"){
    bdb_hiv[i,194] <- ">"  
  }
  if(bla[1] == "<"){
    bdb_hiv[i,194] <- "<"  
  }
}

####### escrever valor
bdb_hiv$value <- 0
for(i in 1:nrow(bdb_hiv)){
  
 bdb_hiv[i, 195] <- parse_number(bdb_hiv[i,10])

}

###############
#### Pegar os com IC50 < 100 um
##############

bdb_hiv_active <- bdb_hiv[which((bdb_hiv$sign == "-" || bdb_hiv$sign == "<") & (bdb_hiv$value <= 1000)),] 
bdb_hiv_inactive <- bdb_hiv[which((bdb_hiv$sign == "-" || bdb_hiv$sign == ">") & (bdb_hiv$value > 1000)),] 

####### Adicionar ral, evg, bic, dtg, cab, allini-1, allini-2
### The "known"

known <- c("CC1=NN=C(O1)C(=O)NC(C)(C)C2=NC(C(=O)NCC3=CC=C(F)C=C3)=C(O)C(=O)N2C",
           "COC1=CC=2N(C=C(C(O)=O)C(=O)C2C=C1CC3=CC=CC(Cl)=C3F)C(CO)C(C)C",
           "CC1CCOC2N1C(=O)C3=C(C(=O)C(=CN3C2)C(=O)NCC4=CC=C(C=C4F)F)O",
           "CC1COC2N1C(=O)C3=C(C(=O)C(=CN3C2)C(=O)NCC4=CC=C(C=C4F)F)O",
           "C1=C(C=C(C(=C1F)CNC(=O)C2=CN3C(=C(C2=O)O)C(=O)N4C5CCC(C5)OC4C3)F)F",
           "COC(C(O)=O)C1=C(C)N=C2C=CC(Br)=CC2=C1C3=CC=C(Cl)C=C32",
           "CC1=NC2=CC=C(Br)C=C2C(C3=CC=C(Cl)C=C3)=C1C(OC(C)(C)C)C(O)=O")


active <- c(bdb_hiv_active$Ligand.SMILES, known)
ic50_active <- c(bdb_hiv_active$value, rep(0.000001, rep(length(known))))
                
#act <- unique(sort(active))
#act2 <- active
act <- parse.smiles(active)
for(i in 1:length(act)){
  if(is.null(act[[i]]) == TRUE){
    a <- i 
  }
}
act <- act[-a]
ic50_active <- ic50_active[-a]

act2 <- vector()
for(i in 1:length(act)){
  act2 <- c(act2, get.smiles(act[i][[1]]))
}




# act2 <- vector()
#  for(i in 1:length(act)){
#    act2 <- c(act2, get.smiles(act[i][[1]]))
#  }
write.csv(act2, '~/ML_integrase_keras/all_over/active.smiles', quote = FALSE, row.names = FALSE)

ina <- parse.smiles(bdb_hiv_inactive$Ligand.SMILES)
ic50_inactive <- c(bdb_hiv_inactive$value)


for(i in 1:length(ina)){
  if(is.null(ina[[i]]) == TRUE){
    a <- i 
  }
}
ina <- ina[-a]
ic50_inactive <- ic50_inactive[-a]

inactive <- vector()
for(i in 1:length(ina)){
  inactive <- c(inactive, get.smiles(ina[i][[1]]))
}


inc <- unique(sort(inactive))
write.csv(inactive, '~/ML_integrase_keras/all_over/inactive.smiles', quote = FALSE, row.names = FALSE)

##############################
##############################
####  Filtrar aprovados
##############################
##############################

are_app <- vector()
approved$is_app <- NA
for(i in 1:nrow(approved)){

bla<-  strsplit(approved[i, 4], split=";")[[1]]
if("approved" %in% bla){
  approved[i,18] <- TRUE
  }  
}

approved_true <- approved[which(approved$is_app == TRUE),]
app_index <- app_index[which(approved$is_app == TRUE)]
#app2 <- approved_true$SMILES
app <- parse.smiles(approved_true$SMILES)
app <- app[-c(2174)]
app_index <- app_index[-2174]

app2 <- vector()
for(i in 1:length(app)){
    
    app2 <- c(app2, get.smiles(app[i][[1]]))
}

appnull <- which(app2 == "")
app2[appnull] <- NA
app2 <- na.omit(app2)
app_index <- app_index[-appnull]

#app2 <- app2[-c(10, 379, 387, 795, 1813, 2157)]

#write.csv(app2, '~/ML_integrase_keras/all_over/approved.smiles', quote = FALSE, row.names = FALSE)

#app2[-c(105,383, 392,801,1820)]
#################
#### read decoys
#################
#### reads decoy lists of each compound and concatenates them
###### 


list <- read.table("~/ML_integrase_keras/all_over/dude-decoys/decoys/decoy_number.txt")
list <- list[which(list$V1 != 1),2]
decoys <- vector()
for(a in list){
  dec <- read.csv(paste("~/ML_integrase_keras/all_over/dude-decoys/decoys/", a, sep=""), skip=1, sep="\t")[,1]
  decoys <- c(decoys, dec)
} 
write.csv(decoys, "~/ML_integrase_keras/all_over/decoys.smiles", quote = FALSE, row.names = FALSE)


#write.csv(app2, '~/ML_integrase_keras/all_over/approved.smiles', quote = FALSE, row.names = FALSE)
#write.csv(act2, '~/ML_integrase_keras/all_over/active.smiles', quote = FALSE, row.names = FALSE)
#write.csv(inactive, '~/ML_integrase_keras/all_over/inactive.smiles', quote = FALSE, row.names = FALSE)
#write.csv(decoys, "~/ML_integrase_keras/all_over/decoys.smiles", quote = FALSE, row.names = FALSE)

############################
###### read and process NPA molecules
############################

npa <- read.csv("~/ML_integrase_keras/all_over/natural/NPAtlas_download.tsv", sep="\t")

nat <- c(as.character(npa$compound_smiles))
nat2 <- nat

# nat2 <- vector()
# for(i in 1:length(nat)){
#   nat2 <- c(nat2, get.smiles(nat[i][[1]]))
# }

########
## remove the ones in the 
########

# toremove <- c()
# toremove <- c(toremove, which(nat2 %in% active))
# toremove <- c(toremove, which(nat2 %in% app2))
# toremove <- c(toremove, which(nat2 %in% inactive))
# toremove <- c(toremove, which(nat2 %in% decoys))
# 
# ### redo, now without the repeated smiles
# remove(nat)
# remove(nat2)
# 
# npa <- npa[-toremove,]
# write.csv(npa, "~/ML_integrase_keras/all_over/natural/npa_tmp.csv")
# 
# 
# 
# npa <- read.csv("~/ML_integrase_keras/all_over/natural/npa_tmp.csv", row.names = 1)
# 
# 
# nat <- c(as.character(npa$compound_smiles))
# nat <- parse.smiles(nat)
# 
# nat2 <- vector()
# for(i in 1:length(nat)){
#   nat2 <- c(nat2, get.smiles(nat[i][[1]]))
# }

#write.csv(npa, "~/ML_integrase_keras/all_over/natural/npa_uniq.csv")
write.csv(nat2, "~/ML_integrase_keras/all_over/natural/npa_smiles.csv", quote = FALSE, row.names = FALSE)
write.csv(app2, '~/ML_integrase_keras/all_over/approved.smiles', quote = FALSE, row.names = FALSE)
write.csv(act2, '~/ML_integrase_keras/all_over/active.smiles', quote = FALSE, row.names = FALSE)
write.csv(inactive, '~/ML_integrase_keras/all_over/inactive.smiles', quote = FALSE, row.names = FALSE)
write.csv(decoys, "~/ML_integrase_keras/all_over/decoys.smiles", quote = FALSE, row.names = FALSE)
write.csv(app_index, "~/ML_integrase_keras/all_over/app_index.csv", quote = FALSE, row.names = FALSE)

write.csv(ic50_active, "~/ML_integrase_keras/all_over/active_ic50.csv", quote = FALSE, row.names = FALSE)
write.csv(ic50_inactive, "~/ML_integrase_keras/all_over/inactive_ic50.csv", quote = FALSE, row.names = FALSE)

######################
##### verify uniqueness
######################

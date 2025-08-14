###
# Data wrangling for phosphocite kinase databases
###


library(tidyverse)
library(rentrez) #for taxid
library(readODS) #for dealing with multiple sheets of spreadsheet
library(stringr)


db_col_names<-c("Substrate","Substrate_Acc","Motif_Sequence","position","Kinase","Kinase_Acc","kinase_group","kinase_family","kinase_sub_family","binding_domain","taxonomy_id","species_name","db_source","SIGNOR_Score")

mdf <- data.frame(matrix(ncol = length(db_col_names), nrow = 0))
colnames(mdf) <- db_col_names

mdf[] <- lapply(mdf, as.character)



#-----------prepare phosphoELM----------------------
phosphoELMpath<-"/media/master/Data/PhosphoData/phosphoELM_all_2015-04.csv"
phosphoELM<-read_tsv(phosphoELMpath)

#Extract the Motif Sequence (15 AA) from the protein substrate sequence (PhosphoELM only contains full sequence)
phosphoELM$sequenceshort<-substr(phosphoELM$sequence,phosphoELM$position-7,phosphoELM$position+7)

#Leftpad the sequence if position is at start otherwise right pad 
phosphoELM$sequenceshort <- ifelse(phosphoELM$position < 8,
                  str_pad(phosphoELM$sequenceshort, width = 15, side = "left", pad = "_"),
                  str_pad(phosphoELM$sequenceshort, width = 15, side = "right", pad = "_"))

#Remove entries where kinase is not known
filtered_phosphoELM <- phosphoELM[!is.na(phosphoELM$kinases), ]

# get unique species
species<-unique(filtered_phosphoELM$species)

# Function to get Taxonomy ID from NCBI
get_taxid <- function(species_name) {
  result <- entrez_search(db="taxonomy", term=species_name)
  if (length(result$ids) > 0) {
    return(result$ids[1])  # Return the first matched taxonomic ID
  } else {
    return(NA)
  }
}

#Query NCBI for tax_ids
tax_ids <- sapply(species, get_taxid)

#map taxid to species
tax_map <- data.frame(
  species = species,
  tax_id = tax_ids
)

# Merge by species to have TAXids in the dataset
filtered_phosphoELM_taxed <- merge(filtered_phosphoELM, tax_map, by = "species", all.x = TRUE)  

#rename columns for merging
filtered_phosphoELM_taxed <- filtered_phosphoELM_taxed %>%
  rename(
    Substrate_Acc = acc,
    Motif_Sequence = sequenceshort,
    Kinase=kinases,
    species_name=species,
    taxonomy_id=tax_id
  )

#add source
filtered_phosphoELM_taxed$db_source<-"phosphoELM"


# Add missing columns with NA
filtered_phosphoELM_taxed[, setdiff(names(mdf), names(filtered_phosphoELM_taxed))] <- NA

# Combine PhosphoELM with masterdf
filtered_phosphoELM_taxed[] <- lapply(filtered_phosphoELM_taxed, as.character)
mdf <- bind_rows(mdf, filtered_phosphoELM_taxed %>% select(db_col_names))


unique(mdf$Kinase)

#---------- prepare SIGNOR -----------

signorpath<-"/media/master/Data/PhosphoData/Signor/M_musculus_phosphorylations_24_02_25.tsv"
signor<-read_tsv(signorpath)

#removes repeated (identical) rows
uniq_signor<-distinct(signor)

#get position
uniq_signor$position<-substr(uniq_signor$RESIDUE,4,nchar(uniq_signor$RESIDUE))

#rename columns for merging
uniq_signor <- uniq_signor %>%
  rename(
    Substrate_Acc = IDB,
    Motif_Sequence = SEQUENCE,
    Kinase=ENTITYA,
    taxonomy_id=TAX_ID,
    Kinase_Acc=IDA,
    Substrate=ENTITYB,
    SIGNOR_Score=SCORE
  )

# Add source
uniq_signor$db_source<-"SIGNOR"

# Add missing columns with NA
uniq_signor[, setdiff(names(mdf), names(uniq_signor))] <- NA

# Combine SIGNOR with masterdf
uniq_signor[] <- lapply(uniq_signor, as.character) #need to convert to character (it recognizes taxid as numbers)
mdf <- bind_rows(mdf, uniq_signor %>% select(db_col_names))

head(mdf)
sort(unique(filtered_phosphoELM_taxed$Kinase))


# ------------- prepare PTMsig

Ptmsigpath<-"/media/master/Data/PhosphoData/data_PTMsigDB_all_sites_v2.0.0.ods"
#read sheets in seperately
sheet1 <- read_ods(Ptmsigpath, sheet=1)
sheet2 <-read_ods(Ptmsigpath, sheet=2)
sheet3 <- read_ods(Ptmsigpath, sheet=3)

#sheets correspond to organisms, add that information
sheet1$taxonomy_id<-9606
sheet2$taxonomy_id<-10090
sheet3$taxonomy_id<-10116

Ptmsig<-rbind(sheet1,sheet2,sheet3)

filtered_Ptmsig<-Ptmsig %>% filter(category=="KINASE-iKiP"|category=="KINASE-PSP")

#get kinases from the signature (has the form: "KINASE-PSP_kinasename") problem is it often has multiple names 
filtered_Ptmsig <- filtered_Ptmsig %>% mutate(Kinase = str_extract(signature, "(?<=_).*"))

#extract Substrate acession
filtered_Ptmsig <- filtered_Ptmsig %>% mutate(position= str_extract(site.uniprot,"(?<=;.).*"))

#extract position
filtered_Ptmsig <- filtered_Ptmsig %>% mutate(Substrate_Acc=str_extract(site.uniprot,".+?(?=;)"))


#rename columns for merging
filtered_Ptmsig <- filtered_Ptmsig %>%
  rename(
    Motif_Sequence = site.flanking,
    db_source = category
  )

#add missing columns
filtered_Ptmsig[, setdiff(names(mdf), names(filtered_Ptmsig))] <- NA

View(filtered_Ptmsig)

#merge data frames
filtered_Ptmsig[] <- lapply(filtered_Ptmsig, as.character) 
mdf <- bind_rows(mdf, filtered_Ptmsig %>% select(db_col_names))

mdfcopy<-mdf
mdf<-mdfcopy

#2016 non human non mouse non rat entries / ~60k
sum(na.omit(mdf$taxonomy_id)!="10090"&na.omit(mdf$taxonomy_id)!="10116"&na.omit(mdf$taxonomy_id)!="9606") 


#discard non human/mouse/rat entries
mdf_humr<-mdf[(mdf$taxonomy_id=="10090"|mdf$taxonomy_id=="10116"|mdf$taxonomy_id=="9606")&!is.na(mdf$taxonomy_id),]

sort(unique(mdf_humr$Kinase))





#clean up kinase nomenclature
mdf_humr$Kinase <- toupper(mdf_humr$Kinase)

# remove - characters
mdf_humr$Kinase <- sub("-","",mdf_humr$Kinase)



#replace written out notation
mdf_humr$Kinase <- gsub("_ALPHA", "A", mdf_humr$Kinase)
mdf_humr$Kinase <- gsub("_BETA", "B", mdf_humr$Kinase)  
mdf_humr$Kinase <- gsub("_GAMMA", "C", mdf_humr$Kinase)
mdf_humr$Kinase <- gsub("_DELTA", "D", mdf_humr$Kinase)  
mdf_humr$Kinase <- gsub("_EPSILON", "E", mdf_humr$Kinase)
mdf_humr$Kinase <- gsub("_ETA", "H", mdf_humr$Kinase)
mdf_humr$Kinase <- gsub("_IOTA", "I", mdf_humr$Kinase)
mdf_humr$Kinase <- gsub("_THETA", "T", mdf_humr$Kinase)
mdf_humr$Kinase <- gsub("_ZETA", "Z", mdf_humr$Kinase)

sort(unique(mdf_humr$Kinase))
sort(unique(mdf_humr[!mdf_humr$Kinase %in% kinhubref$`HGNC Name`&!mdf_humr$Kinase %in% kinhubref$xName&!mdf_humr$Kinase %in% kinhubref$`Manning Name`,]$Kinase))

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("_GROUP", mdf_humr$Kinase)]<-sub("_.*", "",mdf_humr$Kinase[grepl("_GROUP", mdf_humr$Kinase)])
#remove the group from kinase as we only know group
mdf_humr$Kinase[grepl("_GROUP", mdf_humr$Kinase)] <- NA

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<AKT\\>", mdf_humr$Kinase)]<-mdf_humr$Kinase[grepl("\\<AKT\\>", mdf_humr$Kinase)]
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<AKT\\>", mdf_humr$Kinase)] <- NA

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<AMPK\\>", mdf_humr$Kinase)]<-mdf_humr$Kinase[grepl("\\<AMPK\\>", mdf_humr$Kinase)]
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<AMPK\\>", mdf_humr$Kinase)] <- NA

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<BMPR1A/1B/2\\>", mdf_humr$Kinase)]<-"STKR"
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<BMPR1A/1B/2\\>", mdf_humr$Kinase)] <- NA


#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<PIM\\>", mdf_humr$Kinase)]<-"PIM"
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<PIM\\>", mdf_humr$Kinase)] <- NA

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<ERK1/2\\>", mdf_humr$Kinase)]<-"MAPK"
mdf_humr$kinase_sub_family[grepl("\\<ERK1/2\\>", mdf_humr$Kinase)]<-"ERK1"
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<ERK1/2\\>", mdf_humr$Kinase)] <- NA

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<JNK\\>", mdf_humr$Kinase)]<-"MAPK"
mdf_humr$kinase_sub_family[grepl("\\<JNK\\>", mdf_humr$Kinase)]<-"JNK"
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<JNK\\>", mdf_humr$Kinase)] <- NA


#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<LATS1/2\\>", mdf_humr$Kinase)]<-"NDR"
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<LATS1/2\\>", mdf_humr$Kinase)] <- NA

#deal with kinases that are not actually known (only family is known)
mdf_humr$kinase_family[grepl("\\<PAK\\>", mdf_humr$Kinase)]<-"STE20"
#remove the group from kinase as we only know  AKT group
mdf_humr$Kinase[grepl("\\<PAK\\>", mdf_humr$Kinase)] <- NA


mdf_humr

sort(unique(mdf_humr$Kinase))


#choose the gene name if multiple names are given (ptmsig sometimes has ProteinName/GeneName)
#checks if the kinase entry is from ptmsig and has one and only one "/" if so it takes the entry after the backslash otherwise takes teh entire name.
mdf_humr$Kinase <- ifelse(grepl("/", mdf_humr$Kinase) & nchar(gsub("[^/]", "", mdf_humr$Kinase)) == 1 & mdf_humr$db_source %in% c("KINASE-iKiP","KINASE-PSP"),
                      sub(".*/", "", mdf_humr$Kinase),  # Extract after "/"
                      mdf_humr$Kinase)

sort(unique(mdf_humr$Kinase))
sort(unique(mdf_humr[!mdf_humr$Kinase %in% kinhubref$`HGNC Name`&!mdf_humr$Kinase %in% kinhubref$xName&!mdf_humr$Kinase %in% kinhubref$`Manning Name`,]$Kinase))


#unify individual kinases
mdf_humr$Kinase[grepl("TITIN KINASE", mdf_humr$Kinase)]<-"TTN"

mdf_humr$Kinase[grepl("LYNA.LYN", mdf_humr$Kinase)]<-"LYN"
mdf_humr$Kinase[grepl("LYNB.LYN", mdf_humr$Kinase)]<-"LYN"

mdf_humr$Kinase[grepl("TROPOMYOSIN KINASE", mdf_humr$Kinase)]<-"NTRK1" 

mdf_humr$Kinase[grepl("2810408M09RIK", mdf_humr$Kinase)]<-"TP53RK" 

mdf_humr$Kinase[grepl("EG3 KINASE", mdf_humr$Kinase)]<-"MELK" 

mdf_humr$Kinase[grepl("AURORA A", mdf_humr$Kinase)]<-"AURKA" 
mdf_humr$Kinase[grepl("AURORA B", mdf_humr$Kinase)]<-"AURKB" 

mdf_humr$Kinase[grepl("CAMKIA", mdf_humr$Kinase)]<-"CAMK1"
mdf_humr$Kinase[grepl("CAMKIIA", mdf_humr$Kinase)]<-"CAMK2A" 
mdf_humr$Kinase[grepl("CAMKIID", mdf_humr$Kinase)]<-"CAMK2D"
mdf_humr$Kinase[grepl("CAMKIV", mdf_humr$Kinase)]<-"CAMK4"
mdf_humr$Kinase[grepl("CAMKKA", mdf_humr$Kinase)]<-"CAMKK1"

mdf_humr$Kinase[grepl("CAK COMPLEX", mdf_humr$Kinase)]<-"CDK7"

# BCR ABL is a fusion gene i dont know how to treat it.
mdf_humr$Kinase[grepl("BCRABL", mdf_humr$Kinase)]<-"ABL1"

mdf_humr$Kinase[mdf_humr$Kinase == "CDC2CCNB1"]<-"CDC2"
mdf_humr$Kinase[mdf_humr$Kinase == "CDC7DBF4"]<-"CDC7"
mdf_humr$Kinase[mdf_humr$Kinase == "CDK14CCND3"]<-"CDK14"
mdf_humr$Kinase[mdf_humr$Kinase == "CDK1CCNB1"]<-"CDK1"
mdf_humr$Kinase[mdf_humr$Kinase == "CYCLINB/CDK1"]<-"CDK1"
mdf_humr$Kinase[mdf_humr$Kinase == "CYCLINA2/CDK2"]<-"CDK2"
mdf_humr$Kinase[mdf_humr$Kinase == "CYCLINE/CDK2"]<-"CDK2"
mdf_humr$Kinase[mdf_humr$Kinase == "CYCLINC/CDK3"]<-"CDK3"
mdf_humr$Kinase[mdf_humr$Kinase == "CYCLINE1/CDK3"]<-"CDK3"
mdf_humr$Kinase[mdf_humr$Kinase == "CYCLIND/CDK4"]<-"CDK4"
mdf_humr$Kinase[grepl("CDK2CCN", mdf_humr$Kinase)]<-"CDK2"
mdf_humr$Kinase[grepl("CDK3CCN", mdf_humr$Kinase)]<-"CDK3"
mdf_humr$Kinase[grepl("CDK4CCN", mdf_humr$Kinase)]<-"CDK4"
mdf_humr$Kinase[grepl("CDK5", mdf_humr$Kinase)]<-"CDK5"
mdf_humr$Kinase[grepl("CDK6", mdf_humr$Kinase)]<-"CDK6"
mdf_humr$Kinase[grepl("CDK7", mdf_humr$Kinase)]<-"CDK7"
mdf_humr$Kinase[grepl("CDK8", mdf_humr$Kinase)]<-"CDK8"
mdf_humr$Kinase[grepl("CDK9", mdf_humr$Kinase)]<-"CDK9"
mdf_humr$Kinase[grepl("CDK12", mdf_humr$Kinase)]<-"CDK12"
mdf_humr$Kinase[grepl("CDK13", mdf_humr$Kinase)]<-"CDK13"
mdf_humr$Kinase[grepl("CDK14", mdf_humr$Kinase)]<-"CDK14"
mdf_humr$Kinase[grepl("CDK16", mdf_humr$Kinase)]<-"CDK16"

mdf_humr$Kinase[mdf_humr$Kinase == "GRK2"]<-"BARK1"
mdf_humr$Kinase[mdf_humr$Kinase == "GRK3"]<-"BARK2"
mdf_humr$Kinase[mdf_humr$Kinase == "GSK3B/AXIN/APC"]<-"GSK3B"


mdf_humr$Kinase[mdf_humr$Kinase == "CGK2"]<-"PRKG2"
mdf_humr$Kinase[mdf_humr$Kinase == "CK2A"]<-"CSNK2A1"
mdf_humr$Kinase[mdf_humr$Kinase == "CK2B"]<-"CSNK2A2"

mdf_humr$Kinase[mdf_humr$Kinase == "MAP3K20"]<-"ZAK"
mdf_humr$Kinase[mdf_humr$Kinase == "MAP3K21"]<-"MLK4"

mdf_humr$Kinase[mdf_humr$Kinase == "MLCK"]<-"MYLK3"

mdf_humr$Kinase[mdf_humr$Kinase == "MTORC1"]<-"MTOR" #unsure about these
mdf_humr$Kinase[mdf_humr$Kinase == "MTORC2"]<-"MTOR"

mdf_humr$Kinase[mdf_humr$Kinase == "NPM1:ALK"]<-"ALK"
mdf_humr$Kinase[mdf_humr$Kinase == "NPMALK/NPM/ALK"]<-"ALK"


mdf_humr$Kinase[mdf_humr$Kinase == "PKG1/CGKI"]<-"PKG1"
mdf_humr$Kinase[mdf_humr$Kinase == "PKG2/CGKII"]<-"PKG2"

mdf_humr$Kinase[grepl("PRKCB", mdf_humr$Kinase)]<-"PKCB"
mdf_humr$Kinase[grepl("RIPK5", mdf_humr$Kinase)]<-"DSTYK"
mdf_humr$Kinase[grepl("RSK5", mdf_humr$Kinase)]<-"MSK1"


mdf_humr_copy<-mdf_humr

mdf_humr

mdf_humr<-mdf_humr[mdf_humr$Kinase!="ACP1"|is.na(mdf_humr$Kinase),] #is a phosphatase?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CALCINEURIN"|is.na(mdf_humr$Kinase),] #not a kinase?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CCDPK"|is.na(mdf_humr$Kinase),] #idk what this is
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CDC14A"|is.na(mdf_humr$Kinase),] #a phosphatase? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CDC14B"|is.na(mdf_humr$Kinase),] #a phosphatase? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CDC25A"|is.na(mdf_humr$Kinase),] #a phosphatase? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CDC25B"|is.na(mdf_humr$Kinase),] #a phosphatase? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CDC25C"|is.na(mdf_humr$Kinase),] #a phosphatase? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CDKN3"|is.na(mdf_humr$Kinase),] #a phosphatase? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CKM COMPLEX"|is.na(mdf_humr$Kinase),] #idk? 
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CRK"|is.na(mdf_humr$Kinase),] #idk wrong organism?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CTDNEP1"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CTDSP1"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CTDSP2"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="CTDSPL"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="DLG1"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[!grepl("DUSP", mdf_humr$Kinase),] #phosphatases
mdf_humr<-mdf_humr[mdf_humr$Kinase!="EYA1"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="EYA3"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="GM20390"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="HRAS"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="IKBKG"|is.na(mdf_humr$Kinase),] #kinsae inhibitor?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="IKKCOMPLEX"|is.na(mdf_humr$Kinase),] #idk?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="IL6ST"|is.na(mdf_humr$Kinase),] #kinase activator subunit
mdf_humr<-mdf_humr[mdf_humr$Kinase!="KRAS"|is.na(mdf_humr$Kinase),] #kras idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="MEK1/2"|is.na(mdf_humr$Kinase),] #idk,  MAP3K1 is aka MEKK1 but not sure if this is what this means
mdf_humr<-mdf_humr[mdf_humr$Kinase!="MNAT1"|is.na(mdf_humr$Kinase),] #idk cdk activating stuff
mdf_humr<-mdf_humr[mdf_humr$Kinase!="P38"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PDKC"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PDXP"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[!grepl("PHLPP", mdf_humr$Kinase),] #phosphatases
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PHPT1"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PIAS4"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PK"|is.na(mdf_humr$Kinase),] #idk not clear
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PKA"|is.na(mdf_humr$Kinase),] #idk not clear
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PKAA"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PKBB"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PLCG1"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PLCG2"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PP1"|is.na(mdf_humr$Kinase),] #organism?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PP2B"|is.na(mdf_humr$Kinase),] #organism?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PP2CA_R1A_BD"|is.na(mdf_humr$Kinase),] #organism?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PP2CA_R1A_R2A"|is.na(mdf_humr$Kinase),] #organism?
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPEF1"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1A"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1B"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1D"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1E"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1F"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1K"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PPM1L"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[!grepl("PPP", mdf_humr$Kinase),]#phosphatases
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PRKAR2B"|is.na(mdf_humr$Kinase),]#kinase regulator
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PTEN"|is.na(mdf_humr$Kinase),]#phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PTEN"|is.na(mdf_humr$Kinase),]
mdf_humr<-mdf_humr[mdf_humr$Kinase!="PTP4A3"|is.na(mdf_humr$Kinase),]#phopsphatase
mdf_humr<-mdf_humr[!grepl("PTPN", mdf_humr$Kinase),]#phosphatase
mdf_humr<-mdf_humr[!grepl("PTPR", mdf_humr$Kinase),]#phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="RALA"|is.na(mdf_humr$Kinase),] #ras related protein
mdf_humr<-mdf_humr[mdf_humr$Kinase!="RHOA"|is.na(mdf_humr$Kinase),] #gtpase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="ROBO"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="RPAP2"|is.na(mdf_humr$Kinase),] #phosphatase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="RPS6K"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="SPHK2"|is.na(mdf_humr$Kinase),]# SPHK2 sphingosine kinase
mdf_humr<-mdf_humr[mdf_humr$Kinase!="TCEB1"|is.na(mdf_humr$Kinase),] #idk
mdf_humr<-mdf_humr[mdf_humr$Kinase!="TFIIH"|is.na(mdf_humr$Kinase),]
mdf_humr<-mdf_humr[mdf_humr$Kinase!="TGM2"|is.na(mdf_humr$Kinase),] #idk




mdf_humr



# some of the entries are e.g. MAP3K10.MLK2 then we only care for the gene name: MAP3K10.
# checks if the kinase entry is from ptmsig and has one and only one "." if so it takes the entry before the . otherwise takes the entire name.
mdf_humr$Kinase <- ifelse(grepl(".", mdf_humr$Kinase) & nchar(gsub("[^.]", "", mdf_humr$Kinase)) == 1 & mdf_humr$db_source %in% c("KINASE-iKiP","KINASE-PSP"),
                          sub("\\..*", "", mdf_humr$Kinase),  # Extract before "."
                          mdf_humr$Kinase)



#Load kinase family information, data from kinhub (human kinome information)
kinhubpath<-"/media/master/Data/PhosphoData/FamilyClassification/HumanKinomeFamilyClassification.csv"
kinhubref<-read_tsv(kinhubpath)

kinhubref$xName<-toupper(kinhubref$xName)
kinhubref$`Manning Name`<-toupper(kinhubref$`Manning Name`)
kinhubref$`HGNC Name`<-toupper(kinhubref$`HGNC Name`)
kinhubref$Group<-toupper(kinhubref$Group)
kinhubref$Family<-toupper(kinhubref$Family)
kinhubref$SubFamily<-toupper(kinhubref$SubFamily)

#unique kinases left
sort(unique(mdf_humr$Kinase))

#kinase things that couldnt get classified but is relevant enough to not toss? 
sort(unique(mdf_humr[!mdf_humr$Kinase %in% kinhubref$`HGNC Name`&!mdf_humr$Kinase %in% kinhubref$xName&!mdf_humr$Kinase %in% kinhubref$`Manning Name`,]$Kinase))
# CSN2KB, FAM20C, NME1, NME2, PFKFB4, PI4K2A, seem fine
# PIK3C3, PIK3CA... all are kinase subunits?
# pkcc mouse kinase?
# pkm rat kinase?
# PRKAB1 PRKAG2 AMPK subunits?

#not that many entries that are still remainign questionable
table(mdf_humr[!mdf_humr$Kinase %in% kinhubref$`HGNC Name`&!mdf_humr$Kinase %in% kinhubref$xName&!mdf_humr$Kinase %in% kinhubref$`Manning Name`,]$Kinase)


# unify family names
mdf_humr$kinase_sub_family[mdf_humr$kinase_family=="AMPK"]<-"AMPK"
mdf_humr$kinase_family[mdf_humr$kinase_family=="AMPK"]<-"CAMKL"
mdf_humr$kinase_family[mdf_humr$kinase_family=="CAMKI"]<-"CAMK1"
mdf_humr$kinase_family[mdf_humr$kinase_family=="CAMKII"]<-"CAMK2"
mdf_humr$kinase_family[mdf_humr$kinase_family=="GSK3"]<-"GSK"
mdf_humr$kinase_sub_family[mdf_humr$kinase_family=="JNK"]<-"JNK"
mdf_humr$kinase_family[mdf_humr$kinase_family=="JNK"]<-"MAPK"
mdf_humr$kinase_family[mdf_humr$kinase_family=="MAP2k"]<-"STE7"
mdf_humr<-mdf_humr[!grepl("MAP3K", mdf_humr$kinase_family),] #map3k spread multiple families
mdf_humr$kinase_sub_family[mdf_humr$kinase_family=="MARK"]<-"MARK"
mdf_humr$kinase_family[mdf_humr$kinase_family=="MARK"]<-"CAMKL"
mdf_humr$kinase_sub_family[mdf_humr$kinase_family=="p38"]<-"p38"
mdf_humr$kinase_family[mdf_humr$kinase_family=="p38"]<-"MAPK"
mdf_humr$kinase_sub_family[mdf_humr$kinase_family=="P70S6K"]<-"RSKp70"
mdf_humr$kinase_family[mdf_humr$kinase_family=="P70S6K"]<-"RSK"
mdf_humr$kinase_family[mdf_humr$kinase_family=="PAK"]<-"STE20"
mdf_humr$kinase_family[mdf_humr$kinase_family=="PKB"]<-"AKT"
mdf_humr$kinase_family[mdf_humr$kinase_family=="PKG/CGK"]<-"PKG"
mdf_humr$kinase_sub_family[mdf_humr$kinase_family=="ROCK"]<-"ROCK"
mdf_humr$kinase_family[mdf_humr$kinase_family=="ROCK"]<-"DMPK"



# kinases that are not accounted for
sort(unique(mdf_humr[!mdf_humr$Kinase %in% kinhubref$`HGNC Name`&!mdf_humr$Kinase %in% kinhubref$xName&!mdf_humr$Kinase %in% kinhubref$`Manning Name`,]$Kinase))
maybetrim<-c("CSNK2B","FAM20C","NME1","NME2","PFKFB4","PI4K2A","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R5","PKCC","PKM","PRKAB1","PRKAG2")

#all kinases
unique(mdf_humr$Kinase)

#kinase counts
table(mdf_humr$Kinase)

# get kinase family, group, subfamily and uniprot id from kinhub dataframe (very ugly and takes forever as it is not using proper methods)
# make sure to not forget about capitalization

for(i in 1:nrow(mdf_humr)){
  for(j in 1:nrow(kinhubref)){
    if (!is.na(mdf_humr$Kinase[i]) && !is.na(kinhubref$xName[j]) && mdf_humr$Kinase[i]==kinhubref$xName[j]){
      mdf_humr$kinase_family[i]<-kinhubref$Family[j]
      mdf_humr$kinase_sub_family[i]<-kinhubref$SubFamily[j]
      mdf_humr$kinase_group[i]<-kinhubref$Group[j]
      mdf_humr$Kinase_Acc[i]<-kinhubref$UniprotID[j]
      break
    }
    if (!is.na(mdf_humr$Kinase[i]) && !is.na(kinhubref$`Manning Name`[j]) && mdf_humr$Kinase[i]==kinhubref$`Manning Name`[j]){
      mdf_humr$kinase_family[i]<-kinhubref$Family[j]
      mdf_humr$kinase_sub_family[i]<-kinhubref$SubFamily[j]
      mdf_humr$kinase_group[i]<-kinhubref$Group[j]
      mdf_humr$Kinase_Acc[i]<-kinhubref$UniprotID[j]
      break
    }
    if (!is.na(mdf_humr$Kinase[i]) && !is.na(kinhubref$`HGNC Name`[j]) && mdf_humr$Kinase[i]==kinhubref$`HGNC Name`[j]){
      mdf_humr$kinase_family[i]<-kinhubref$Family[j]
      mdf_humr$kinase_sub_family[i]<-kinhubref$SubFamily[j]
      mdf_humr$kinase_group[i]<-kinhubref$Group[j]
      mdf_humr$Kinase_Acc[i]<-kinhubref$UniprotID[j]
      break
    }
  }
  print(i)
}

mdf_humr

# get the group from the family (again really ugly code)
for(i in 1:nrow(mdf_humr)){
  for(j in 1:nrow(kinhubref)){
    if (!is.na(mdf_humr$kinase_family[i]) && mdf_humr$kinase_family[i]==kinhubref$Family[j]){
      mdf_humr$kinase_group[i]<-kinhubref$Group[j]
      break
    }
  }
  print(i)
}


#remove entries with unexpected aminoacid entries
filtered_mdf <- mdf_humr[!grepl("X|B|Z", mdf_humr$Motif_Sequence), ]

#filter out entries that couldnt be classified
filtered_mdf <- filtered_mdf[!filtered_mdf$Kinase %in% maybetrim,]

#filter out histidine phosphorylations (only 26)
filtered_mdf<- filtered_mdf[!grepl(".......H.......",filtered_mdf$Motif_Sequence),]

#save annotated df
write_tsv(filtered_mdf,"masterDatabaseall3.tsv")
#filtered_mdf<-read.csv("masterDatabaseall3.tsv",sep="\t")


#count tables
table(filtered_mdf$Kinase)
table(filtered_mdf$kinase_family)
table(filtered_mdf$kinase_group)




#count occurances of AminoAcids
sum(str_count(filtered_mdf$Motif_Sequence,pattern=".......S......."))


#as soon as we group by sequence, we go from 57k observations to 18.8k.
unique(mdf_humr)

length(unique(mdf_humr$Motif_Sequence))
         

filtered_grouped <- filtered_mdf %>% group_by(Motif_Sequence)

#how many kinases per sequence
summarise(filtered_grouped,Kinase_count = n_distinct(Kinase))


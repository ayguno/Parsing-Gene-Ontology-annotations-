###############################################################################
# Author : Ozan Aygun
# Date : 09/27/2016
#
# Purpose : to provide GOCC and other annotations to desired protein/site specific
# data tables in a systematic manner.
#
# Update: 10/03/20116
#         we will update the Uniprot table by using all available entries from UNIPROT
#         and perform two types of merging 1: by UNIPROT entry 2: by Gene Name
###############################################################################

library(dplyr)


# First prepare uniprot query table:

# This ID table is downloaded from UNIPROT FTP site on 10/03/2016:
# FTP site : ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/


ftpdata <-read.table(file = "HUMAN_9606_idmapping.dat",
                          ,sep = "\t", stringsAsFactors = FALSE)

# Over 4 million entries, the first column is Uniprot IDS

UniprotIDs <- data.frame(unique(ftpdata[,1]))
# contains 133996 unique uniprot IDs


# Wrote this into a csv file to use in an Uniprot search
write.table(file = "UniprotHuman10032016.txt", UniprotIDs, row.names = F, quote = F,col.names = F )


# Use this ID file ("UniprotHuman10032016.txt") to perform a Uniprot search and download the matched annotations:
# This proces takes quite some time.

# This query matched 109961 Human entries from Uniprot(Reviewed + Unreviewed)
# This table is downloaded as : uniprot-yourlist%3AM2016100314483A1C7ED25EE8374758DF3FD545FD34EB804.tab 

######################################
# Start here to continue the work
######################################



# Read the downloaded Uniprot table 

setwd("Z:/LabMembers/Ozan/MitoProject")


# We have to delete the FUNCTION column otherwise file is not read properly.

Uniprot <- read.delim(file = "uniprot-yourlist%3AM2016100314483A1C7ED25EE8374758DF3FD545FD34EB804.tab", stringsAsFactors = FALSE, header = TRUE)
#  109961 entries are sucesfully read as we expected.

# This gives a table of 31 annotations for 109961 Unique Uniprot IDs

# Add identifier "Uniprot" in front of the column names of this table.

colnames(Uniprot) <-paste("UNIPROT",colnames(Uniprot), sep = "_")
Mitochondrial_Evidence_UNIPROT=apply(apply(Uniprot,2, grepl,pattern="mitoch|Mitoch|MITOCH"),1,any) 
Uniprot$Mitochondrial_Evidence_UNIPROT <- Mitochondrial_Evidence_UNIPROT

# Compile the gene name column

Uniprot$UNIPROT_Gene.name <- sapply(Uniprot$UNIPROT_Gene.names,function(x){
        temp<-regexpr(" ",x)
        return(substr(x,start = 1, stop = temp-1))
})

w <- which(Uniprot$UNIPROT_Gene.name == "")
Uniprot$UNIPROT_Gene.name[w] <- Uniprot$UNIPROT_Gene.names[w] 


# Read the next reference table, Calvo's list:

Calvo <- read.delim(file="Sarah_Calvo_Entrez2GO_human.gene_manual_curation.txt", stringsAsFactors = F)

# This is a table consisting of 19224 Genes and 34 features

# Add identifier "CALVO" in front of the column names of this table.

colnames(Calvo) <-paste("CALVO",colnames(Calvo), sep = "_")
Mitochondrial_Evidence_CALVO=apply(apply(Calvo,2, grepl,pattern="mitoch|Mitoch|MITOCH"),1,any) 
Calvo$Mitochondrial_Evidence_CALVO <- Mitochondrial_Evidence_CALVO

# Read the last reference table, mitomatrix gold+ list

MitoGoldP <- read.delim(file = "MitomatrixGoldplus_final proteome_and_combination1.txt", stringsAsFactors = FALSE)
# This table contains information for 521 geneIDs
# Only first 3 columns have information;

MitoGoldP <- MitoGoldP[,1:3]
MitoGoldP <- mutate(MitoGoldP, annotation = "Mitochondria")

# Add identifier "MITOGOLDP" in front of the column names of this table.

colnames(MitoGoldP) <- paste("MITOGOLDP", colnames(MitoGoldP), sep = "_")
Mitochondrial_Evidence_MITOGOLDP=apply(apply(MitoGoldP,2, grepl,pattern="mitoch|Mitoch|MITOCH"),1,any) 
MitoGoldP$Mitochondrial_Evidence_MITOGOLDP <- Mitochondrial_Evidence_MITOGOLDP

# Now we can match the identifiers to the desired data sets and annotate them
# for the presence of mitochondrial evidence

# Read the desired tables:

CarrLab_Biotin_SiteCentric_data <- read.csv2(file = "//bennett/msdatasm/Namrata/Apex/20150326_BiotinEnrichment_Rep1_3_480_497_AND_227_Research_QE_35ions/Rep01/proteinPeptideComparisonColumnsExport.VM.1.ssv",
                                              stringsAsFactors = F)



CarrLab_SA_ProteinCentric_data <- read.csv2(file = "//bennett/msdatasm/Namrata/APEX/20160229_Rep1_3_SA_SILAC_Enrichment/Rep1/proteinProteinCentricColumnsExport.3.ssv",
                                            stringsAsFactors = F)



Rhee_etal_SA_ProteinCentric_data <- read.csv2(file = "//bennett/msdatasm/Namrata/APEX/20160218_Rheeetal_MM_data_re_search/Rep01/proteinProteinCentricColumnsExport_fiveDigits.3.ssv",
                                              stringsAsFactors = F)


# All of these tables contain "geneSymbol" and "accession_numbers" columns. Use them to match
# reference tables.

paths <-data.frame(CarrLab_Biotin_SiteCentric_data = "//bennett/msdatasm/Namrata/Apex/20150326_BiotinEnrichment_Rep1_3_480_497_AND_227_Research_QE_35ions/Rep01/proteinPeptideComparisonColumnsExport.VM.1.ssv",
                     CarrLab_SA_ProteinCentric_data="//bennett/msdatasm/Namrata/APEX/20160229_Rep1_3_SA_SILAC_Enrichment/Rep1/proteinProteinCentricColumnsExport.3.ssv",
                     Rhee_etal_SA_ProteinCentric_data= "//bennett/msdatasm/Namrata/APEX/20160218_Rheeetal_MM_data_re_search/Rep01/proteinProteinCentricColumnsExport_fiveDigits.3.ssv" )

for( i in seq_along(colnames(paths))){
        
        sourcedata <- read.csv2(file= as.character(paths[1,i]), stringsAsFactors = FALSE)
        
        gene_protein <- unique(sourcedata %>% select(accession_number, geneSymbol))
        colnames(gene_protein) <- paste(colnames(paths)[i],colnames(gene_protein), sep = "_")
        
        w<-which(gene_protein[,1]== "")
        if(length(w) >0) {gene_protein[w,1] <- "Not Available"}
        w<-which(gene_protein[,2]== "")
        if(length(w) >0) {gene_protein[w,2] <- "Not Available"}
        
        # merge Uniprot evidence by Gene Name
        
        temp <- merge(gene_protein,Uniprot, by.x = paste(colnames(paths)[i],"geneSymbol",sep = "_"), 
                              by.y = "UNIPROT_Gene.name", all.x = TRUE, all.y = FALSE, sort = FALSE)
        
        # add MERGED_BY_GENE_NAME suffix to these columns
        
        colnames(temp) <-paste(colnames(temp),"_MERGED_BY_GENE_NAME",sep = "")
        
        temp <- subset(temp, !duplicated(temp[,1]))
        
        gene_protein <- merge(gene_protein,temp, by.x =paste(colnames(paths)[i],
                                                             "geneSymbol",sep = "_"),
                              by.y = colnames(temp)[1],
                              all.x = TRUE, all.y = FALSE, sort = FALSE)

        
        # merge Uniprot evidence by Uniprot ID
        
        gene_protein <- merge(gene_protein,Uniprot, by.x = paste(colnames(paths)[i],"accession_number",sep = "_"), 
                              by.y = "UNIPROT_Entry", all.x = TRUE, all.y = FALSE, sort = FALSE)
        
        
        # merge Calvo evidence
        
        gene_protein <- merge(gene_protein,Calvo, by.x = paste(colnames(paths)[i],"geneSymbol",sep = "_"), 
                              by.y = "CALVO_Symbol", all.x = TRUE, all.y = FALSE, sort = FALSE)
        
        # merge mitomatrix gold+ list evidence
        
        
        gene_protein <- merge(gene_protein,MitoGoldP, by.x = paste(colnames(paths)[i],"geneSymbol",sep = "_"), 
                              by.y = "MITOGOLDP_Gene.Symbol", all.x = TRUE, all.y = FALSE, sort = FALSE)
        
        
        # compile the data
        ID <- gene_protein[,1:2]
        w <-which(grepl("^Mitochondrial_Evidence", colnames(gene_protein)))
        Evidence <- gene_protein[,w]
        w <- c(1,2,w)
        rest <- gene_protein[,-w]
        
        
        final_Evidence <- data.frame(ID,Evidence)
        Mitochondrial_Evidence_ANY <- apply(Evidence,1,any,na.rm=TRUE)
        
        final_Evidence <- data.frame(final_Evidence,Mitochondrial_Evidence_ANY,rest)
                                     
                                    
        glued <- merge(sourcedata,final_Evidence,by.x = "accession_number", 
                       by.y = paste(colnames(paths)[i],"accession_number",sep = "_"),
                       all.x = TRUE, all.y = FALSE,sort = FALSE)
        
        
        
        
        setwd(dirname(as.character(paths[1,i])))
        if(!exists("GOCC_annotation")){dir.create("GOCC_annotation")}
        
        write.table(final_Evidence,file = paste("./GOCC_annotation/",basename(as.character(paths[1,i])), "_Annotation_Only.txt", sep = ""), sep = "\t", row.names = FALSE)
        write.table(glued,file = paste("./GOCC_annotation/",basename(as.character(paths[1,i])), "_Annotation_data_glued.txt", sep = ""), sep = "\t", row.names = FALSE)
        
        print(paste("Completed: ", dirname(as.character(paths[1,i]))))
        print(paste("Completed: ", colnames(paths)[i]))
        
}







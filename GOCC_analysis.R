###############################################################################
# Author : Ozan Aygun
# Date : 09/27/2016
#
# Purpose : to provide GOCC and other annotations to desired protein/site specific
# data tables in a systematic manner.
#
###############################################################################

library(dplyr)


# First prepare uniprot query table:


specMilldata <-read.table(file = "//bennett/seqdb/UniProt.human.20141017.RNFISnr.150contams",
                          ,sep = "\t", stringsAsFactors = FALSE)
x<-grepl(">",specMilldata[,1])

specMillUniprot <- data.frame(specMilldata[x,])

x <- grepl("HUMAN",specMillUniprot[,1]) 

specMillUniprot <-as.character(specMillUniprot[x,])

x<-regexpr("\\|",specMillUniprot)

specMillUniprot <- substr(specMillUniprot,(x+1),length(specMillUniprot)) 

x<-regexpr("\\|",specMillUniprot)

specMillUniprot <- unique(substr(specMillUniprot,1,x-1)) 


# We have 32976 unique, Human uniprot IDs from Spectrum Mill database
specMillUniprot <- data.frame(UniprotHumanSpectrumMill=specMillUniprot)

# Wrote this into a csv file to use in an Uniprot search
write.csv(file = "UniprotHumanSpectrumMill.csv", specMillUniprot )


# Use this IDs to perform a Uniprot search and download the matched annotations:
# This proces takes quite some time.



######################################
# Start here to continue the work
######################################



# Read the downloaded table 

setwd("Z:/LabMembers/Ozan/MitoProject")


Uniprot <- read.delim(file = "uniprot-HUMAN-09272016.tab", stringsAsFactors = FALSE)


# This gives a table of 31 annotations for 22684 Unique Uniprot IDs

# Add identifier "Uniprot" in front of the column names of this table.

colnames(Uniprot) <-paste("UNIPROT",colnames(Uniprot), sep = "_")
Mitochondrial_Evidence_UNIPROT=apply(apply(Uniprot,2, grepl,pattern="mitoch|Mitoch|MITOCH"),1,any) 
Uniprot$Mitochondrial_Evidence_UNIPROT <- Mitochondrial_Evidence_UNIPROT


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
        
        # merge Uniprot evidence
        
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







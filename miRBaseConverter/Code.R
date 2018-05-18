
#Scenario 1: Constructing Gene Regulatory Network with mRNAs, TFs and miRNAs
## 1. Processing for miRTarbase 7.0 
library("miRBaseConverter")
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/miRTarbase_hsa_7.rda"))
miRNANames=miRTarbase_hsa_7$miRNA
version=checkMiRNAVersion(miRNANames)
##version="miRBase v21"
Accessions=miRNA_NameToAccession(miRNANames,version = version)
miRTarbase_hsa_7=cbind(Accessions,miRTarbase_hsa_7)

## 2. Processing for miRecords 4 
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/miRecords_v4.rda"))
index=which(miRecords_v4$`Target gene_species_scientific` =="Homo sapiens" 
            & miRecords_v4$miRNA_species=="Homo sapiens")
miRecords_hsa_v4=miRecords_v4[index,]
miRNANames=miRecords_hsa_v4$miRNA_mature_ID
version=checkMiRNAVersion(miRNANames)
##version="miRBase v16"
Accessions=miRNA_NameToAccession(miRNANames,version = version)
miRNANames=miRNA_AccessionToName(Accessions$Accession,targetVersion="v21")
colnames(miRNANames)=c("Accession","miRNAName_v21")
## Retrive the ENTREZ ID for target genes
geneSymbols=miRecords_hsa_v4$`Target gene_name`
library(HGNChelper)
res=checkGeneSymbols(geneSymbols)
library("org.Hs.eg.db")
geneSymbols=select(org.Hs.eg.db, keys=res$Suggested.Symbol, 
                  columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
miRecords_hsa_v4=cbind(miRNANames,geneSymbols,miRecords_hsa_v4)
##remove the missing miRNAs and target genes.
index=which(is.na(miRecords_hsa_v4$Accession) | is.na(miRecords_hsa_v4$ENTREZID))
miRecords_hsa_v4=miRecords_hsa_v4[-index,]

## 3. Processing for  transmir_v1.2
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/transmir_v1_2.rda"))
index=which(transmir_v1_2$organism=="human")
transmir_hsa_v1_2=transmir_v1_2[index,]
miRNANames=transmir_hsa_v1_2$mir
miRNANames=paste0("hsa-",miRNANames)
version=checkMiRNAVersion(miRNANames)
##version="miRBase v17"
##Convert the precursor names to the mature miRNA names
miRNANames=miRNA_PrecursorToMature(miRNANames,version="v21")
## Retrieve the miRNA accessions for Mature1( tails with "5p")
miRNANames=miRNANames$Mature1
Accessions=miRNA_NameToAccession(miRNANames,version="v21")
transmir_hsa_v1_2=cbind(Accessions,transmir_hsa_v1_2)
## remove the missing mature miRNA
index=which(is.na(transmir_hsa_v1_2$Accession))
transmir_hsa_v1_2=transmir_hsa_v1_2[-index,]

## 4. Processing for ENCODE_TF_miRNA
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/ENCODE_TF_miRNA.rda"))
##select the TF-miRNA interactions
index=which(ENCODE_TF_miRNA$Interections=="(TF-miRNA)")
ENCODE_TF_miRNA=ENCODE_TF_miRNA[index,]
## Retrieve the miRNA accessions for miRNA names
miRNANames=ENCODE_TF_miRNA$miRNA
version=checkMiRNAVersion(miRNANames)
##version="miRBase v17"
##convert the miRBase version from 17 to 21
Accessions=miRNA_NameToAccession(miRNANames,version = version)
miRNANames=miRNA_AccessionToName(Accessions$Accession,targetVersion = "v21")
colnames(miRNANames)=c("Accession","miRNAName_v21")
##Retrieve the gene ENTREZ ID for gene names
geneSymbols=ENCODE_TF_miRNA$TF
library(HGNChelper)
res=checkGeneSymbols(geneSymbols)
library("org.Hs.eg.db")
geneSymbols=select(org.Hs.eg.db, keys=res$Suggested.Symbol, 
                 columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
##Combine the TF ENTREZ IDs and miRNA accessions
ENCODE_TF_miRNA=cbind(geneSymbols,miRNANames)


## 5. Processing for TRED and ENCODE Chip-seq databases
library("HGNChelper")
library("org.Hs.eg.db")
##TRED
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/TRED.rda"))
TRED=unlist2(TRED)
TRED=data.frame("TF_Symbol"=names(TRED),"targetGene_ENTREZID"=TRED,stringsAsFactors = FALSE)
TF_Symbol=TRED$TF_Symbol
res=checkGeneSymbols(TF_Symbol)
TF=select(org.Hs.eg.db, keys=res$Suggested.Symbol, 
                 columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
colnames(TF)=c("TF_Symbol","TF_ENTREZID")

targetGene_ENTREZID=TRED$targetGene_ENTREZID
targetGene=select(org.Hs.eg.db, keys=targetGene_ENTREZID, 
                 columns=c("SYMBOL","ENTREZID"), keytype="ENTREZID")
colnames(targetGene)=c("targetGene_ENTREZID","targetGene_Symbol")
TRED=cbind(TF,targetGene)
##Remove the TFs with missing ENTREZID
index=which(is.na(TRED$TF_ENTREZID))
TRED=TRED[-index,]

##ENCODE Chip-seq
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/ENCODE_TF_mRNA.rda"))
ENCODE_TF_mRNA=unlist2(ENCODE_TF_mRNA)
ENCODE_TF_mRNA=data.frame("TF_Symbol"=names(ENCODE_TF_mRNA),
                          "targetGene_ENTREZID"=ENCODE_TF_mRNA,stringsAsFactors = FALSE)
TF_Symbol=ENCODE_TF_mRNA$TF_Symbol
res=checkGeneSymbols(TF_Symbol)
TF=select(org.Hs.eg.db, keys=res$Suggested.Symbol, 
          columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
colnames(TF)=c("TF_Symbol","TF_ENTREZID")

targetGene_ENTREZID=ENCODE_TF_mRNA$targetGene_ENTREZID
targetGene=select(org.Hs.eg.db, keys=targetGene_ENTREZID, 
                  columns=c("SYMBOL","ENTREZID"), keytype="ENTREZID")
colnames(targetGene)=c("targetGene_ENTREZID","targetGene_Symbol")
ENCODE_TF_mRNA=cbind(TF,targetGene)
##Remove the TFs with missing ENTREZID
index=which(is.na(ENCODE_TF_mRNA$TF_ENTREZID))
ENCODE_TF_mRNA=ENCODE_TF_mRNA[-index,]


##Case study 1. Constructing miRNA-TF-mRNA network based on the interaction databases

####Extract all the regulation relationships from each of the interaction databases
########miRTarbase_hsa_7
interaction1=cbind("regulatorID"=miRTarbase_hsa_7$Accession,
                   "regulatorSymbol"=miRTarbase_hsa_7$miRNAName_v21,
                   "targetID"=miRTarbase_hsa_7$`Target Gene (Entrez Gene ID)`,
                   "targetSymbol"=miRTarbase_hsa_7$`Target Gene`)
########miRecords_hsa_v4
interaction2=cbind("regulatorID"=miRecords_hsa_v4$Accession,
                   "regulatorSymbol"=miRecords_hsa_v4$miRNAName_v21,
                   "targetID"=miRecords_hsa_v4$ENTREZID,
                   "targetSymbol"=miRecords_hsa_v4$SYMBOL)

########transmir_hsa_v1_2
interaction3=cbind("regulatorID"=transmir_hsa_v1_2$entrezid,
                   "regulatorSymbol"=transmir_hsa_v1_2$gene,
                   "targetID"=transmir_hsa_v1_2$Accession,
                   "targetSymbol"=transmir_hsa_v1_2$miRNAName_v21)

########ENCODE_TF_miRNA
interaction4=cbind("regulatorID"=ENCODE_TF_miRNA$ENTREZID,
                   "regulatorSymbol"=ENCODE_TF_miRNA$SYMBOL,
                   "targetID"=ENCODE_TF_miRNA$Accession,
                   "targetSymbol"=ENCODE_TF_miRNA$miRNAName_v21)
########TRED
interaction5=cbind("regulatorID"=TRED$TF_ENTREZID,
                   "regulatorSymbol"=TRED$TF_Symbol,
                   "targetID"=TRED$targetGene_ENTREZID,
                   "targetSymbol"=TRED$targetGene_Symbol)
########ENCODE_TF_mRNA
interaction6=cbind("regulatorID"=ENCODE_TF_mRNA$TF_ENTREZID,
                   "regulatorSymbol"=ENCODE_TF_mRNA$TF_Symbol,
                   "targetID"=ENCODE_TF_mRNA$targetGene_ENTREZID,
                   "targetSymbol"=ENCODE_TF_mRNA$targetGene_Symbol)

####Combine all the regulation relationships to construct a comprehensive miRNA-TF-mRNA network 
network=rbind(interaction1,interaction2,interaction3,interaction4,interaction5,interaction6)
network=unique(network)
## Write the network in Data table(one to one) format
write.csv(network,file="network.csv", row.names = FALSE)

##Transform the data table network to sparse matrix network
####a. Extract all the molecules
all_molecules=unique(c(network[,2] ,network[,4]))
####b. Generate a symmetric matrix with all "0".
network_matrix=matrix(0,nrow=length(all_molecules),ncol=length(all_molecules))
####c. Fill the regulatory network using the interaction databases
####Make sure the adequate computing resources for this step, 
####otherwise the calculation might collapse. 
index1=match(network[,2],all_molecules )
index2=match(network[,4],all_molecules)
network_matrix[cbind(index1,index2)]=1
sum(network_matrix)
save(network_matrix,all_molecules,file="network_matrix.rda")


#Case study 2: Predicting miRNA Targets by Integrating Different Data Sources

##Convert the miRNA names for TCGA BRCA miRNA expression dataset
load(url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/BRCA200.rda"))
miRNANames=colnames(BRCA_miRNA)
version=checkMiRNAVersion(miRNANames)##miRBase version 16
####Convert the miRBase version 16 to version 21
Accessions=miRNA_NameToAccession(miRNANames,version = version)
BRCA_miRNANames=miRNA_AccessionToName(Accessions$Accession,targetVersion = "v21")
####Remove the dead miRNAs
index=which(is.na(BRCA_miRNANames$TargetName))
BRCA_miRNANames=BRCA_miRNANames[-index,]
BRCA_miRNA=BRCA_miRNA[,-index]
####Convert the precursor to mature names
BRCA_miRNANames=miRNA_PrecursorToMature(BRCA_miRNANames[,2],version = "v21")
####Study the mature miRNAs with tails "5p"
colnames(BRCA_miRNA)=BRCA_miRNANames$Mature1
####Merge the miRNAs with the same name
BRCA_miRNA=t(apply(BRCA_miRNA, 1, function(x) tapply(x, colnames(BRCA_miRNA), mean)))
####Save the BRCA data as .csv file for downstream analysis
BRCA=cbind(BRCA_miRNA,BRCA_mRNA)
write.csv(BRCA,file = "BRCA_412miR_8000mR.csv", row.names = FALSE)


## Process the miRNA Annotations of TargetScan7.0 and miRTabase 7
TargetScan_7=read.csv(
  url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/TargetScan_7.0.csv"),
                      as.is =TRUE)
miRNANames=TargetScan_7$mir
version=checkMiRNAVersion(miRNANames)  ##miRBase version 17
####convert the miRNA names to the miRBase version 21
miRNANames=miRNAVersionConvert(miRNANames,targetVersion = "v21")
TargetScan_7$mir=miRNANames$TargetName
index=which(is.na(TargetScan_7$mir))
TargetScan_7=TargetScan_7[-index,]
write.csv(TargetScan_7,file = "TargetScan7_miRv21.csv", row.names = FALSE)

####Extract the interactions from miRTarbase_7
load(
  url("https://taoshengxu.github.io/ExperimentData/miRBaseConverter/miRTarbase_hsa_7.rda"))
miRTarbase7=miRTarbase_hsa_7[,c(2,4)]
write.csv(miRTarbase7,file = "miRTarbase7.csv", row.names = FALSE)

####Predict miRNA Targets
library("miRLAB")
Pridict_targets=Pearson("./BRCA_412miR_8000mR.csv", c(1:412), c(413:8412)
                        , targetbinding="./TargetScan7_miRv21.csv")
#Extract the hsa-miR-224-5p as the example
miRTop100=bRank(Pridict_targets, 151,100, TRUE)

####Validate the prediction using miRTarBase 7.0
miRTop100Confirmed = Validation(miRTop100, "miRTarbase7.csv")
miRTop100Confirmed=miRTop100Confirmed[[1]]
#### 8 predicted miRNA-mRNAs interactions are confirmed in miRTarbase 7.0.
## miRNA mRNA Correlation
## hsa-miR-224-5p BCL2 -0.45568915
## hsa-miR-224-5p FOSB -0.18087381
## hsa-miR-224-5p DIO1 -0.17142808
## hsa-miR-224-5p SLC12A5 -0.12755148
## hsa-miR-224-5p ZKSCAN8 -0.09462372
## hsa-miR-224-5p EYA4 -0.08509764
## hsa-miR-224-5p PTX3 0.25557907
## hsa-miR-224-5p EFNA3 0.22785473

####Report the sequence of the miRNA of interest
#Query the sequence of hsa-miR-224-5p and hsa-miR-224-3p
Accessions=miRNA_NameToAccession(c("hsa-miR-224-5p","hsa-miR-224-3p"),version = "v22")
getMiRNASequence(Accessions$Accession, targetVersion = "v22")
##      Accession         miRNASequence_v22
## 1 MIMAT0000281 UCAAGUCACUAGUGGUUCCGUUUAG
## 2 MIMAT0009198   AAAAUGGUGCCCUAGUGACUACA




getwd()
library(biomaRt)
ensembl <- useEnsembl(
    biomart = "genes",
    dataset = "rnorvegicus_gene_ensembl",
    mirror = "asia"
)
??biomart
filters = listFilters(ensembl)
filters[1:5,]
filters
View(filters)
library(org.Rn.eg.db)
ensembl_ID=toTable(org.Rn.egENSEMBLTRANS)
head(ensembl_ID)
utr3=getSequence(
    id=ensembl_ID$trans_id,
    type = "ensembl_transcript_id",
    seqType = "3utr",
    mart = ensembl
)
outfile <- file("rno-3utr.fa", "w")
??file
for (i in 1:nrow(utr3)) {
  h = paste(c(">", utr3[i,2]), collapse="")
  writeLines(h, outfile)
  writeLines(utr3[i,1], outfile)
}

close(outfile)
head(outfile)
outfile

head(toTable(org.Rn.egENSEMBL2EG))
head(toTable(org.Rn.egENSEMBLTRANS))
head(toTable(org.Rn.egSYMBOL2EG))
head(toTable(org.Rn.egGENENAME))
symbol=toTable(org.Rn.egSYMBOL2EG)
library(tidyverse)
library(dplyr)
bt=dplyr::inner_join(symbol,ensembl_ID,by = "gene_id")
head(bt)
View(bt)
ENSRNOT00000029565=getSequence(
    id = "ENSRNOT00000029565",
    type = "ensembl_transcript_id",
    seqType = "3utr",
    mart = ensembl
)

getSequence(id = "ENSRNOT00000096679",type = "ensembl_transcript_id",seqType = "3utr",mart = ensembl)

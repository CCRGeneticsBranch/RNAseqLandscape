library(dplyr)
library(jsonlite)
## functions
flattenDF <- function(df=NULL){
  flatdf <- as.data.frame( do.call(rbind, apply(df, 1, function(x){
    each_row <- lapply(x, function(y){
      if(length(y) > 1 ) { 
        return_char = paste(y, collapse = ";")
        return(return_char)
      } else { 
        return(y)
        }
    })
    return(each_row)
  })
  )) %>% mutate_all(as.character)
  return(flatdf)
}

othersJson <- fromJSON("./other.json")
otherHGNCTable <- othersJson$response$docs %>% 
        dplyr::select(one_of("symbol","ensembl_gene_id","gene_family",
                             "name","refseq_accession"))
otherHGNCTableFlat<- flattenDF(df=otherHGNCTable)
write.table(otherHGNCTableFlat, "otherHGNCTableFlat.txt", sep="\t", row.names = FALSE)
saveRDS(otherHGNCTableFlat, "otherHGNCTableFlat.RDS")



pcJson <- fromJSON("./protein-coding_gene.json")
pcHGNCTable <- pcJson$response$docs %>% 
  dplyr::select(one_of("symbol","ensembl_gene_id","gene_family",
                       "name","refseq_accession"))
pcHGNCTableFlat<- flattenDF(df=pcHGNCTable)
write.table(pcHGNCTableFlat, "pcHGNCTableFlat.txt", sep="\t", row.names = FALSE)
saveRDS(pcHGNCTableFlat, "pcHGNCTableFlat.RDS")

pc.other.HGNCTableFlat <-rbind(pcHGNCTableFlat,otherHGNCTableFlat)
saveRDS(pc.other.HGNCTableFlat, "PC.other.HGNCTableFlat.rds")
write.table(pc.other.HGNCTableFlat, "PC.other.HGNCTableFlat.txt", sep="\t", row.names = FALSE)


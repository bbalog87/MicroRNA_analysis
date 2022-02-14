library(data.table)
library(dplyr)


### Import all samples
mirna_files<-list.files(pattern="*.csv")
mirna_files_name<-gsub(pattern = "\\.csv$", "", mirna_files) ## remove extension
mirna_data <-lapply(mirna_files, fread)


### Add sampleName for each entry
for (i in 1:length(mirna_files_name)){
  
  mirna_data[[i]]$SampleName<-mirna_files_name[i]
}


## Drop empty columns, and rename the remaining ones
mirna_NA<- c("V3", "V4","V10","V11","V12")
mirna_COLNAMES<-c("provisional_id", "mirDeep2_Score", "totalReads_count",
                  "matureReads_count", "loopReads_count","starReads_count",
                  "significance_Pvalue","NCBI_blast", "consensus_mature_sequence",
                  "consensus_star_sequence","consensus_precursor_sequence",
                  "precursor_coordinate", "sampleName")



## merge all samples in a master file and  keep only high confident candidate miRNAs
mirna_masterFile<-do.call("rbind", mirna_data) %>% 
  select(-c("V3", "V4","V10","V11","V12")) %>% rename_all(~mirna_COLNAMES) %>% 
  filter(log2(mirDeep2_Score)>2 & matureReads_count>10)
mirna_masterFile<- mirna_masterFile %>% filter(significance_Pvalue=="yes")



mirna_masterFile$MappingRate<-nchar(mirna_masterFile$consensus_mature_sequence)/nchar(mirna_masterFile$consensus_star_sequence)

## Get occurrence in distinct samples
for (i in 1:nrow(mirna_masterFile)) {
  OCC<-mirna_masterFile %>% 
    filter(consensus_mature_sequence==consensus_mature_sequence[i]) %>% 
    distinct(sampleName) %>% nrow()
  mirna_masterFile$DistinctOccurence[i]<-OCC
} 
## keep candidates occurring in at least 2 distinct samples
mirna_masterFile2<-mirna_masterFile %>% filter(DistinctOccurence>=2)
# mapping rate
mirna_candidates_filtered<-mirna_masterFile2[!duplicated(mirna_masterFile2$consensus_mature_sequence),]


### filtering
for (i in 1:nrow(mirna_masterFile_filtered_MV_1)) {
  OCC<-mirna_masterFile_filtered_MV_1 %>% 
    filter(consensus_mature_sequence==consensus_mature_sequence[i]) %>% 
    distinct(sampleName) %>% nrow()
  mirna_masterFile_filtered_MV_1$SamplesOccMature[i]<-OCC
} 

write.csv(mirna_masterFile_filtered_MV_1,
          "mirna_masterFile_filtered_MV_1_SamplesCount.txt",
          row.names = F)



mirna.merged.final<-merge(mirna_masterFile_filtered_MV_1, 
                          mirna_masterFile, by=c("provisional_id","consensus_mature_sequence",
                                                 "consensus_precursor_sequence","precursor_coordinate",
                                                 "sampleName","consensus_star_sequence"))


for (i in 1:nrow(mirna.merged.final)){
  mirna.merged.final$SampleListPrecursor[i]<-as.vector(mirna.merged.final %>% 
    filter(consensus_precursor_sequence==consensus_precursor_sequence[i]) %>% 
      distinct(sampleName))
} 

df2 <- apply(mirna.merged.final,2,as.character)

write.csv(df2,"mirna.merged.final_sampleNames.txt",row.names = F)


df2<-df2 %>%
  group_by(precursor_coordinate) %>%                           
  summarise(totalStarReads = sum(starReads_count), totalMatureReads = sum(matureReads_count)) %>%
  filter(DistinctOccurence>1)


# write final masterf ile

write.csv(df2, "list_provisionalID-PrecursorSequence.txt", row.names = F)

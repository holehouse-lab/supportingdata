# 25-09-2024

# Note: retrieve of long tax lists from NCBI can generate a known error (Error: Bad Request (HTTP 400) regardless of the library used: https://github.com/ropensci/taxa/issues/202 - The error is not related to a specific species identifier nor quantity of previously retrieved species - Random stuff! Work around with a loop function

#install.packages("taxize")

library("taxize")

#Sys.setenv(ENTREZ_KEY = "f4675e56feaad5b8b7bd18b0b4507934e108") # if needed

#############################
# Retrieve taxa information #
#############################

# load data file
listTaxoIDs <- read.table("UPID_to_TaxID.tsv")

# select the column with IDs
list_of_tax <- listTaxoIDs$V2
list_of_tax <- list_of_tax[1:100]

# The loop bellow will retrieve info from ncbi based on a list of numeric identifiers. As errors are expected (Error: Bad Request (HTTP 400), once we get it we wait and try again. At the end merge the retrieved list and proceed. 
#classes_all = NULL

for (i in 1:length(list_of_tax)) {
  
  print(i)
  
  classes_i <- try(classification(list_of_tax[i], db = "ncbi"))
  if (class(classes_i)=="try-error") {
    Sys.sleep(10)
    classes_i <- try(classification(list_of_tax[i], db = "ncbi"))}
  classes_all <- c(classes_all, classes_i)
  
}

###########################################################################
# Extract family and species information from the generated dataset above #
###########################################################################

# From the bulk data above, transform that in a dataframe with spp & family information
library(dplyr)
#current_list <- classes_all[3427:3430]

current_list <- classes_all

df_species_family = NULL

pb = txtProgressBar(min = 0, max = length(current_list), initial = 0)

for (i in 1:length(current_list)) {
  
  #print(current_list[i])
  
  setTxtProgressBar(pb,i)
  
  # if this error was saved in the list, skip it to the next item in the loop
  #if("Error : Bad Request (HTTP 400)\n" %in% current_list[[i]]) next
  
  # as there other errors stored in the list and these have length == 1, skip anything with that length.
  if(length(current_list[[i]]) == 1) next
  
  vector_family <- bind_rows(current_list[i]) %>% 
    filter(rank == 'family')
  
  vector_species <- bind_rows(current_list[i]) %>% 
    filter(rank == 'species')
  
  # if the vector family is empty, populate it to avoid rbind error 
  if (nrow(vector_family) == 0) {
    vector_family <- data.frame(name="empty",
                                rank="empty",
                                id="empty")
    
    
  }
  
  df_species_family_current <- cbind(vector_species,vector_family)
  df_species_family <- rbind(df_species_family,df_species_family_current)
  
}
close(pb)

# rename the columns
colnames(df_species_family) <- c("spp", "rank_spp", "id_spp", "family", "rank_family", "id_family")

##################################################################
# Add metadata info to the table generated in the previous loop. #
##################################################################

# load metadata
lineages_DF <- read.csv("lineages_R.csv", quote="", header = TRUE)
listTaxoIDs$V2 <- as.character(listTaxoIDs$V2)
#listTaxoIDs[12569,]

# Add the information about the NCBI ids to the table generated
merged_DF <- merge(x = listTaxoIDs, y = df_species_family, by.x = "V2", by.y = "id_spp", all = FALSE)

# Add the information about intrinsic disorder
merged2_DF <- merge(x = merged_DF, y = lineages_DF, by.x = "V1", by.y = "Proteome_Id", all = FALSE)
merged2_DF$Average_fraction_per_protein <- as.numeric(merged2_DF$Average_fraction_per_protein) 

# remove virus from the list
merged2_noVirus_DF <- merged2_DF[rowSums(merged2_DF == "virus")==0, , drop = FALSE]

# remove duplicates
merged2_noVirus_unique_DF <- merged2_noVirus_DF %>% distinct()


######################
# Summary per family #
######################

# Get a mean of the intrinsic disorder measure per family, and represent each family by the first occurring species. More: https://stackoverflow.com/questions/40630269/dplyr-summarize-by-string
summary_DF <- merged2_noVirus_unique_DF %>%
  group_by(family) %>%
  summarize(mean_TFAP = mean(Total_fraction_across_proteome), sd_TFAP=sd(Total_fraction_across_proteome),
            spp = first(spp))

write.csv(summary_DF, "MeanPerFamily_TFAP.csv", row.names = FALSE)

summary2_DF <- merged2_noVirus_unique_DF %>%
  group_by(family) %>%
  summarize(mean_AFPP = mean(Average_fraction_per_protein), sd_AFPP=sd(Average_fraction_per_protein),
            spp = first(spp))

write.csv(summary2_DF, "MeanPerFamily_AFPP.csv", row.names = FALSE)


####################################################
# Create a tree with one representative per family #
####################################################

# Add NCBI identifier to the summary df so the tree can be retrieved
merged2_noVirus_unique_DF <- merged2_noVirus_unique_DF %>% select(V1, V2,spp,Organism)
summary_taxID_DF <- left_join(x = summary_DF, y = merged2_noVirus_unique_DF, by = "spp", keep = FALSE)

# remove duplicates
unique_DF <- summary_taxID_DF %>% distinct()

# select the column with IDs
list_of_tax_short <- unique_DF$V2
#list_of_tax <- list_of_tax[1:100]

# Retrieve info from ncbi based on a list of numeric identifiers. As errors are expected (Error: Bad Request (HTTP 400), once we get it we wait and try again. 
Sys.setenv(ENTREZ_KEY = "f4675e56feaad5b8b7bd18b0b4507934e108") # if needed
taxize_class2 <- classification(list_of_tax_short, db = "ncbi")

# convert to tree
taxize_tree_short <- class2tree(taxize_class2, check = TRUE)

# change name of file
entryBaseName="list1480"

# plot
fileNameToSavePDF <- paste0(entryBaseName, "_", "SpeciesTree.pdf")
pdf(file=fileNameToSavePDF, width=100, height=500)
plot(taxize_tree_short)
dev.off()

#Pull the phylo object from the classtree for plotting
tax_phylo <- taxize_tree_short$phylo

library("ape")

#Now that the tree works lets write it into a newick file in case
write.tree(tax_phylo, file = "1480Families_tree")


# end

#write.csv(do.call("rbind",list(unlist(classes_all))),"testme.csv")
#test_read <- list(read.csv("testme.csv"))

#readLines("testme.csv")

#write.csv(classes_all,"testme2.csv")

#read.csv()

# merge lists of results before plotting
taxize_tree <- class2tree(classes_all, check = TRUE)

# change name of file
entryBaseName="list30"

# plot
fileNameToSavePDF <- paste0(entryBaseName, "_", "SpeciesTree.pdf")
pdf(file=fileNameToSavePDF, width=10, height=24)
plot(taxize_tree)
dev.off()


# end

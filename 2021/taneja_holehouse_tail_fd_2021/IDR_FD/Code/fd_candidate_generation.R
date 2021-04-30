library(dplyr)
library(tidyr)

pdb_data = read.csv('../Processed_Data/PDB_Level_Info_iter0_6828.csv')
idr_data = read.csv('../Processed_Data/idr_metadata.csv')

colnames(pdb_data)[1] = 'Uniprot_ID'

idr_data_within_X = idr_data %>% filter(abs(FD_IDR_terminal_dist) <= 10) %>% as.data.frame()

idr_data_within_X_filtered = idr_data_within_X %>% left_join(pdb_data, by = 'Uniprot_ID') %>% as.data.frame() %>% filter(Method %in% c("X-ray diffraction", "Solution NMR")) %>% as.data.frame()

unique(idr_data_within_X_filtered$Uniprot_ID)

for (i in 1:nrow(idr_data_within_X_filtered))
{
  if (grepl(';',idr_data_within_X_filtered[i,'Residues_Covered']))
  {
    semicolon_loc =  gregexpr(pattern =';',idr_data_within_X_filtered[i,'Residues_Covered'])[[1]][1]
    idr_data_within_X_filtered[i,'Residues_Covered'] = substr(idr_data_within_X_filtered[i,'Residues_Covered'],1,semicolon_loc-1)
  }
}

idr_data_within_X_filtered = idr_data_within_X_filtered %>% separate(Residues_Covered, c('Start_Coverage_FD_PDB', 'End_Coverage_FD_PDB'), sep = '-')

idr_data_within_X_filtered$Start_Coverage_FD_PDB = as.numeric(idr_data_within_X_filtered$Start_Coverage_FD_PDB)
idr_data_within_X_filtered$End_Coverage_FD_PDB = as.numeric(idr_data_within_X_filtered$End_Coverage_FD_PDB)

idr_data_within_X_filtered$start_coverage_error = abs(idr_data_within_X_filtered$Start_Coverage_FD - idr_data_within_X_filtered$Start_Coverage_FD_PDB)
idr_data_within_X_filtered$end_coverage_error = abs(idr_data_within_X_filtered$End_Coverage_FD - idr_data_within_X_filtered$End_Coverage_FD_PDB)

idr_data_within_X_filtered$total_coverage_error = idr_data_within_X_filtered$start_coverage_error + idr_data_within_X_filtered$end_coverage_error

#there are multiple pdb entries per protein; each has different coverage
#Start_Coverage_FD and End_Coverage_FD refer to the cumulative coverage for a protein across all structures
#we want to pick the structure with the 'best' coverage
idr_data_within_X_filtered = idr_data_within_X_filtered %>% group_by(Uniprot_ID) %>% filter(total_coverage_error == min(total_coverage_error)) %>% filter(start_coverage_error == min(start_coverage_error)) %>% filter(end_coverage_error == min(end_coverage_error)) %>% as.data.frame()

idr_data_within_X_at_pdb_level = data.frame()
##
for (i in 1:nrow(idr_data_within_X_filtered))
{
  n_term_status = idr_data_within_X_filtered[i,'N_Terminal_Status']
  c_term_status = idr_data_within_X_filtered[i,'C_Terminal_Status']
  start_idr = idr_data_within_X_filtered[i,'Start_Coverage_IDR']
  end_idr = idr_data_within_X_filtered[i,'End_Coverage_IDR']
  
  start_fd_pdb = idr_data_within_X_filtered[i,'Start_Coverage_FD_PDB']
  end_fd_pdb = idr_data_within_X_filtered[i,'End_Coverage_FD_PDB']
  fd_len_pdb = end_fd_pdb - start_fd_pdb
  
  if (n_term_status) {
    if (abs(start_fd_pdb-end_idr) < 10) {
      fd_idr_terminal_dist_pdb = start_fd_pdb-end_idr
      idr_data_within_X_at_pdb_level = rbind(idr_data_within_X_at_pdb_level, cbind(idr_data_within_X_filtered[i,], data.frame('FD_IDR_terminal_dist_PDB' = fd_idr_terminal_dist_pdb, 'FD_Length_PDB' = fd_len_pdb)))
    }
  }
  
  if (c_term_status) {
    if (abs(end_fd_pdb-start_idr) < 10) {
      fd_idr_terminal_dist_pdb = start_idr-end_fd_pdb
      idr_data_within_X_at_pdb_level = rbind(idr_data_within_X_at_pdb_level, cbind(idr_data_within_X_filtered[i,], data.frame('FD_IDR_terminal_dist_PDB' = fd_idr_terminal_dist_pdb, 'FD_Length_PDB' = fd_len_pdb)))
    }
  }
  
}

#get unique entries
idr_data_within_X_at_pdb_level_unique = idr_data_within_X_at_pdb_level %>% group_by(Uniprot_ID, N_Terminal_Status, C_Terminal_Status) %>% filter(nchar(Chain) == min(nchar(Chain))) %>% filter(PDB_ID == min(PDB_ID)) %>% as.data.frame()

idr_data_within_X_at_pdb_level_unique$FD_coverage_percentage = signif(idr_data_within_X_at_pdb_level_unique$FD_Length_PDB/(idr_data_within_X_at_pdb_level_unique$Sequence_Length - idr_data_within_X_at_pdb_level_unique$IDR_Length),2)

write.csv(idr_data_within_X_at_pdb_level_unique, '../Processed_Data/fd_candidate.csv', row.names = FALSE)

library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(stringr)

#this script is for paper

#setup parallel backend to use many processors
cores=detectCores()
registerDoParallel(cores-1)


fastaFile = readAAStringSet("../Raw_Data/uniprot_human_proteome_reviewed_2020_05_clean.fasta")
pdb_data = read.csv('../Processed_Data/PDB_Level_Info_iter0_6828.csv')

seq_name = names(fastaFile)
sequence = paste(fastaFile)

uniprot_id = unlist(lapply(seq_name, function(x) { return(strsplit(x, '|',fixed=TRUE)[[1]][2]) }))

fasta_df = data.frame(uniprot_id, sequence)
colnames(fasta_df) = c('Uniprot_ID', 'Sequence')
fasta_df$Sequence = as.character(fasta_df$Sequence)
fasta_df = fasta_df %>% group_by(Uniprot_ID) %>% mutate(Sequence_Length = nchar(Sequence)) %>% as.data.frame()

idr_df = read.table(file = '../Raw_Data/DISORDER_R3_human_proteome_2019_10_clean_no_comma_shephard_domains.tsv', sep = '\t', header = FALSE)
colnames(idr_df) = c('Uniprot_ID', 'Start_Coverage', 'End_Coverage', 'Domain')
idr_df = idr_df[,-c(which(colnames(idr_df) == 'Domain'))]
colnames(idr_df)[2:3] = paste(colnames(idr_df)[2:3],'IDR',sep='_')

idr_df = inner_join(idr_df, fasta_df[,c('Uniprot_ID', 'Sequence_Length')], by = 'Uniprot_ID')


fd_mutually_exclusive_df = read.csv(file = '../Processed_Data/uniprot_human_proteome_folded_domain_mutally_exclusive_residues.csv', header = TRUE)
colnames(fd_mutually_exclusive_df)[2:3] = paste(colnames(fd_mutually_exclusive_df)[2:3],'FD',sep='_')

#determine whether an IDR is part of a cterminal or nterminal region 
nterminal_domain = c() 
cterminal_domain = c()
nontail_domain = c()

for (i in 1:nrow(idr_df))
{
  nterminal_status = 0
  cterminal_status = 0
  nontail_status = 0 
  
  curr_start_coverage = idr_df[i,'Start_Coverage_IDR']
  curr_end_coverage = idr_df[i,'End_Coverage_IDR']
  
  coverage_len = curr_end_coverage - curr_start_coverage
  
  curr_sequence_len = idr_df[i,'Sequence_Length']
  
  if ( (curr_start_coverage == 1) & (coverage_len < curr_sequence_len) ) {
    nterminal_status = 1
  }
  
  else if ( (curr_end_coverage == curr_sequence_len) & (coverage_len < curr_sequence_len) ) {
    cterminal_status = 1
  }
  
  else if (coverage_len < curr_sequence_len) {
    nontail_status = 1
  }
  
  
  nterminal_domain = c(nterminal_domain, nterminal_status)
  cterminal_domain = c(cterminal_domain, cterminal_status)
  nontail_domain = c(nontail_domain, nontail_status)
}

idr_df$N_Terminal_Status = nterminal_domain
idr_df$C_Terminal_Status = cterminal_domain
idr_df$Non_Tail_Status = nontail_domain

fd_idr_terminal_dist = c()
start_coverage_fd = c()
end_coverage_fd = c() 

for (i in 1:nrow(idr_df))
{
  curr_nterminal_status = idr_df[i,'N_Terminal_Status']
  curr_cterminal_status = idr_df[i,'C_Terminal_Status']
  curr_uniprot_id = idr_df[i,'Uniprot_ID']
  
  if (curr_nterminal_status | curr_cterminal_status)
  {
    curr_start_coverage = idr_df[i,'Start_Coverage_IDR']
    curr_end_coverage = idr_df[i,'End_Coverage_IDR']
    
    relevant_fd_df = fd_mutually_exclusive_df %>% filter(Uniprot_ID == curr_uniprot_id) %>% as.data.frame()
    
    if (curr_nterminal_status) {
      anchor_residue = curr_end_coverage
      relevant_fd_df = relevant_fd_df %>% mutate(dist = Start_Coverage_FD - anchor_residue) %>% as.data.frame()
    } else {
      anchor_residue = curr_start_coverage
      relevant_fd_df = relevant_fd_df %>% mutate(dist = anchor_residue - End_Coverage_FD) %>% as.data.frame()
    }
    
    min_dist_df = relevant_fd_df %>% filter(abs(dist) == min(abs(dist))) %>% as.data.frame()
    fd_idr_terminal_dist = c(fd_idr_terminal_dist, min_dist_df[1,'dist'])
    start_coverage_fd = c(start_coverage_fd, min_dist_df[1,'Start_Coverage_FD'])
    end_coverage_fd = c(end_coverage_fd, min_dist_df[1,'End_Coverage_FD'])
  }
  else
  {
    fd_idr_terminal_dist = c(fd_idr_terminal_dist, NA)
    start_coverage_fd = c(start_coverage_fd, NA)
    end_coverage_fd = c(end_coverage_fd, NA)
  }
}


idr_df$FD_IDR_terminal_dist = fd_idr_terminal_dist
idr_df$Start_Coverage_FD = start_coverage_fd
idr_df$End_Coverage_FD = end_coverage_fd

idr_df = idr_df %>% mutate(IDR_Length = End_Coverage_IDR - Start_Coverage_IDR, FD_Length = End_Coverage_FD - Start_Coverage_FD) %>% as.data.frame()

#remove entries with NA values in dist -- this occurs because there is no PDB derived folded domain for that protein 
idr_df_filtered = idr_df %>% filter((!(is.na(FD_IDR_terminal_dist)) & ((N_Terminal_Status == 1) | (C_Terminal_Status == 1))) | (Non_Tail_Status == 1)) %>% as.data.frame()


#n terminal or c terminal
nterminal_idr_df = idr_df_filtered %>% filter(N_Terminal_Status == 1) %>% as.data.frame()
cterminal_idr_df = idr_df_filtered %>% filter(C_Terminal_Status == 1) %>% as.data.frame()

write.csv(idr_df_filtered, paste('../Processed_Data/idr_metadata.csv',sep='/'),row.names=FALSE)





###############################
#Calculate Charge distribution for each IDR

calc_charge_sequence = function(sequence)
{
  charge_seq = c()
  for (i in 1:nchar(sequence))
  {
    curr_aa = substr(sequence,i,i)
    if (curr_aa %in% c('R','K'))
    {
      charge_seq = c(charge_seq, 1)
    }
    else if (curr_aa %in% c('E','D'))
    {
      charge_seq = c(charge_seq, -1)
    }
    else
    {
      charge_seq = c(charge_seq, 0)
    }
  }
  return(charge_seq)
}


idr_charge_df = list()

for (i in 1:nrow(idr_df_filtered))
{
  curr_nterminal_status = idr_df_filtered[i,'N_Terminal_Status']
  curr_cterminal_status = idr_df_filtered[i,'C_Terminal_Status']
  curr_nontail_status = idr_df_filtered[i,'Non_Tail_Status']
  
  curr_uniprot_id = idr_df_filtered[i,'Uniprot_ID']
  curr_start_coverage = idr_df_filtered[i,'Start_Coverage_IDR']
  curr_end_coverage = idr_df_filtered[i,'End_Coverage_IDR']
  curr_dist = idr_df_filtered[i,'FD_IDR_terminal_dist']
  
  curr_sequence = (fasta_df %>% filter(uniprot_id == curr_uniprot_id) %>% select(Sequence))[1,]
  
  curr_idr_sequence = substr(curr_sequence, curr_start_coverage, curr_end_coverage)
  curr_charge_sequence = calc_charge_sequence(curr_idr_sequence)
  
  unique_id = paste(curr_uniprot_id, curr_start_coverage, curr_end_coverage, sep = '_')
  
  if (curr_nterminal_status | curr_cterminal_status | curr_nontail_status)
  {
    idr_charge_df[[unique_id]] = curr_charge_sequence
  }
  else
  {
    print('error')
  }
}



get_kappa_subset_seq = function(seq, start, end)
{
  if (start <= end) {
    return(calculate_kappa(substr(seq,start,end)))
  } else {
    return(calculate_kappa(substr(seq,end,start)))
  }
}


#for (i in 1:nrow(idr_df_filtered))
idr_charge_stats = foreach(i=1:nrow(idr_df_filtered), .combine=rbind) %dopar% 
{
    print(i)
    curr_nterminal_status = idr_df_filtered[i,'N_Terminal_Status']
    curr_cterminal_status = idr_df_filtered[i,'C_Terminal_Status']
    curr_nontail_status = idr_df_filtered[i,'Non_Tail_Status']
    
    curr_uniprot_id = idr_df_filtered[i,'Uniprot_ID']
    curr_sequence = (fasta_df %>% filter(uniprot_id == curr_uniprot_id) %>% select(Sequence))[1,]
    curr_start_coverage = idr_df_filtered[i,'Start_Coverage_IDR']
    curr_end_coverage = idr_df_filtered[i,'End_Coverage_IDR']
    curr_dist = idr_df_filtered[i,'FD_IDR_terminal_dist']
    
    unique_id = paste(curr_uniprot_id, curr_start_coverage, curr_end_coverage, sep = '_')
    
    curr_charge_sequence = idr_charge_df[[unique_id]]
    
    num_residues = length(curr_charge_sequence)
    num_positive = length(which(curr_charge_sequence == 1))
    num_negative = length(which(curr_charge_sequence == -1))
    fp = num_positive/num_residues
    fn = num_negative/num_residues
    ncpr = (num_positive - num_negative)/num_residues
    tcpr = (num_positive + num_negative)/num_residues

    ret = data.frame(Uniprot_ID = curr_uniprot_id, unique_id = unique_id,
                     N_Terminal_Status = curr_nterminal_status,
                     C_Terminal_Status = curr_cterminal_status,
                     Non_Tail_Status = curr_nontail_status,
                     FD_IDR_terminal_dist = curr_dist,
                     start_pos = curr_start_coverage, 
                     end_pos = curr_end_coverage,
                     num_residues = num_residues,
                     fp = fp,
                     fn = fn,
                     ncpr = ncpr,
                     tcpr = tcpr)
    return(ret)
}


write.csv(idr_charge_stats, paste('../Processed_Data','idr_charge_distribution_stats_all.csv',sep='/'),row.names=FALSE)













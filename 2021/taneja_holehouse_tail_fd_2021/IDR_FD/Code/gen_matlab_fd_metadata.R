library(bio3d)
library(dplyr)
#match pdb sequence number based on fasta sequence 
align_fasta_pdb_seq = function(curr_uniprot_id, pdb_data, curr_chain)
{
  pdb_seq = pdbseq(pdb_data, atom.select(pdb_data, elety="CA", chain=curr_chain))
  
  pdb_seq_string = c()
  
  #get pdb seq string for chain 
  for (i in 1:length(pdb_seq))
  {
    if (i > 1) {
      prev_res_num = names(pdb_seq)[i-1]
      curr_res_num = names(pdb_seq)[i]
      if (curr_res_num != prev_res_num) {
        pdb_seq_string = c(pdb_seq_string, pdb_seq[i])
      }
    } else {
      pdb_seq_string = c(pdb_seq[i])
    }
  }
    
  pdb_seq_string = paste(pdb_seq_string,collapse="")
  
  fasta_sequence = (fasta_df %>% filter(uniprot_id == curr_uniprot_id) %>% select(Sequence))[1,]

  start_loc_list = c()
  end_loc_list = c()
  
  pdb_start_loc_list = c()
  pdb_end_loc_list = c()
  
  pdb_end_pos = 0
  
  #get actual locations of residues 
  while (pdb_end_pos < (nchar(pdb_seq_string)-10))
  {
    search_length = min(nchar(pdb_seq_string)-pdb_end_pos, 10)
    
    pdb_start_pos = pdb_end_pos+1
    
    start_loc_found = FALSE
    for (i in pdb_start_pos:(nchar(pdb_seq_string)-10))
    {
      match_start = gregexpr(substr(pdb_seq_string, i, i+search_length), fasta_sequence)[[1]]
      start_loc = gregexpr(substr(pdb_seq_string, i, i+search_length), fasta_sequence)[[1]][length(match_start)]
      pdb_start_pos = i 
      
      if (start_loc != -1)
      {
        start_loc_found = TRUE
        start_loc_list = c(start_loc_list, start_loc)
        pdb_start_loc_list = c(pdb_start_loc_list, pdb_start_pos)
        
        for (j in rev(pdb_start_pos+1:nchar(pdb_seq_string)))
        {
          end_loc = gregexpr(substr(pdb_seq_string, i, j), fasta_sequence)[[1]]
          if (end_loc != -1)
          {
            match_len = attr(gregexpr(substr(pdb_seq_string, i, j), fasta_sequence)[[1]], 'match.length')-1
            end_loc = start_loc + match_len
            end_loc_list = c(end_loc_list, end_loc)
            pdb_end_pos = pdb_start_pos + match_len  
            pdb_end_loc_list = c(pdb_end_loc_list, pdb_end_pos)
            break 
          }
        }
        break 
      }
    }
    if (start_loc_found == FALSE) {
      pdb_end_pos = pdb_end_pos+1
    }
  }
  
  if (length(pdb_start_loc_list) > 0)
  {
    pdb_dict = list()
    for (i in 1:length(pdb_start_loc_list)) 
    {
      pdb_start = pdb_start_loc_list[i]
      pdb_end = pdb_end_loc_list[i]
      
      idx = 0 
      for (j in pdb_start:pdb_end) {
        pdb_dict[[j]] = start_loc_list[i] + idx
        idx = idx+1
      }    
      
    }
    
    atom_data = pdb_data$atom
    relevant_atom_data = atom_data %>% filter(chain == curr_chain) %>% as.data.frame()
    relevant_atom_data$resno_updated = relevant_atom_data$resno
    
    prev_resno = curr_resno = relevant_atom_data[1,'resno']
    idx = 1
    
    for (i in 1:nrow(relevant_atom_data))
    {
      curr_resno = relevant_atom_data[i,'resno']
      
      if (i > 1) {
        prev_resno = relevant_atom_data[i-1,'resno']
      }
      
      if (curr_resno != prev_resno) {
        idx = idx+1
      } 
      
      if(is.null(pdb_dict[idx][[1]]))
      {
        relevant_atom_data[i,'resno_updated'] = NA
      }
      else
      {
        relevant_atom_data[i,'resno_updated'] = pdb_dict[idx][[1]]
      }
      
    }
    
    return(relevant_atom_data)
  }
  else
  {
    return(data.frame())
  }

}




fastaFile = readAAStringSet("./Raw_Data/uniprot_human_proteome_reviewed_2020_05_clean.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
uniprot_id = unlist(lapply(seq_name, function(x) { return(strsplit(x, '|',fixed=TRUE)[[1]][2]) }))
fasta_df = data.frame(uniprot_id, sequence)
colnames(fasta_df) = c('Uniprot_ID', 'Sequence')
fasta_df$Sequence = as.character(fasta_df$Sequence)
fasta_df = fasta_df %>% group_by(Uniprot_ID) %>% mutate(Sequence_Length = nchar(Sequence)) %>% as.data.frame()

fd_candidate_data = read.csv('./Processed_Data/fd_candidate.csv')
idr_charge_stats = read.csv('./Processed_Data/idr_charge_distribution_stats_all.csv')
fd_charge_stats = read.csv('./Processed_Data/fd_charge_distribution_stats_all.csv')

fd_candidate_data = fd_candidate_data %>% arrange(PDB_ID) %>% as.data.frame()
#fd_candidate_data = fd_candidate_data %>% filter(PDB_ID %in% c('5LN3', '5OCH', '6R7I', '6U62') == FALSE) %>% as.data.frame() #pdb ID with multiple chains 



matlab_metadata = data.frame()


for (i in 1:nrow(fd_candidate_data))
{
  print(i)
  curr_uniprot_id = fd_candidate_data[i,'Uniprot_ID']
  curr_pdb_id = tolower(fd_candidate_data[i,'PDB_ID'])
  n_terminal = fd_candidate_data[i,'N_Terminal_Status']
  c_terminal = fd_candidate_data[i,'C_Terminal_Status']
  fd_coverage_percentage = fd_candidate_data[i,'FD_coverage_percentage']
  curr_chain_all = fd_candidate_data[i,'Chain']
  start_cov_idr = fd_candidate_data[i,'Start_Coverage_IDR']
  end_cov_idr = fd_candidate_data[i,'End_Coverage_IDR']
  start_cov_fd = fd_candidate_data[i,'Start_Coverage_FD']
  end_cov_fd = fd_candidate_data[i,'End_Coverage_FD']
  chain_list = strsplit(curr_chain_all,',')[[1]]
  
  for (curr_chain in chain_list)
  {
    curr_chain = gsub(" ", "", curr_chain, fixed = TRUE)
    idr_charge_data = idr_charge_stats %>% filter(Uniprot_ID == curr_uniprot_id) %>% filter(start_pos == start_cov_idr) %>% filter(end_pos == end_cov_idr) %>% as.data.frame()
    
    fp_idr = signif(idr_charge_data[,'fp'],2)
    fn_idr = signif(idr_charge_data[,'fn'],2)
    ncpr_idr = signif(idr_charge_data[,'ncpr'],2)
    tcpr_idr = signif(idr_charge_data[,'tcpr'],2)
    idr_len = idr_charge_data[,'num_residues']
    
    fd_charge_data = fd_charge_stats %>% filter(Uniprot_ID == curr_uniprot_id) %>% filter(start_pos == start_cov_fd) %>% filter(end_pos == end_cov_fd) %>% as.data.frame()
    
    if (nrow(fd_charge_data) == 0) {
      print('0 ROWS IN FD_CHARGE_DATA')
    }
    
    fp_fd = signif(fd_charge_data[,'fp'],2)
    fn_fd = signif(fd_charge_data[,'fn'],2)
    ncpr_fd = signif(fd_charge_data[,'ncpr'],2)
    tcpr_fd = signif(fd_charge_data[,'tcpr'],2)
    fd_len = fd_charge_data[,'num_residues']
    
    tail = '' 
    
    if (n_terminal == 1) { 
      tail = 'N_terminal'
      idr_start = fd_candidate_data[i,'End_Coverage_IDR']
      fd_start_reported = fd_candidate_data[i,'Start_Coverage_FD']
    } else {
      tail = 'C_terminal'
      idr_start = fd_candidate_data[i,'Start_Coverage_IDR']
      fd_start_reported = fd_candidate_data[i,'End_Coverage_FD']
    }
    
    pdb_filename = paste0('./PDB_Structure_FD_candidates_within10_FD/PDB_Files/',paste0(curr_pdb_id,'_',curr_chain,'.pdb'))
    

    pdb_data = read.pdb(pdb_filename)
    
    
    if ((length(unique(pdb_data$atom$chain)) == 1) && length(unique(pdb_data$atom$resno)) > 10) 
    {
  
      relevant_atom_data = align_fasta_pdb_seq(curr_uniprot_id, pdb_data, curr_chain)

      if (nrow(relevant_atom_data) > 0)
      {
        if (n_terminal == 1) { 
          fd_actual_start = min(relevant_atom_data$resno_updated, na.rm=TRUE)
          fd_pdb_res_num = relevant_atom_data %>% filter(resno_updated == fd_actual_start) %>% select(resno) %>% as.data.frame()
          fd_pdb_res_num = fd_pdb_res_num[1,]
        } else {
          fd_actual_start = max(relevant_atom_data$resno_updated, na.rm=TRUE)
          fd_pdb_res_num = relevant_atom_data %>% filter(resno_updated == fd_actual_start) %>% select(resno) %>% as.data.frame()
          fd_pdb_res_num = fd_pdb_res_num[1,]
        }
    
        #this condition is not satisfied due to error in sequence number on PDB 
        if (abs(idr_start - fd_actual_start) <= 10)
        {
          relevant_atom_data_nearest_residue = relevant_atom_data %>% filter(resno_updated == fd_actual_start) %>% filter(elety == 'CA') %>% as.data.frame()
          
          if (nrow(relevant_atom_data_nearest_residue) == 0) {
            relevant_atom_data_nearest_residue = relevant_atom_data %>% filter(resno_updated == fd_actual_start) %>% as.data.frame()
            relevant_atom_data_nearest_residue = relevant_atom_data_nearest_residue[1,]
          }
          
          #x_loc = relevant_atom_data_nearest_residue[1,'x']
          #y_loc = relevant_atom_data_nearest_residue[1,'y']
          #z_loc = relevant_atom_data_nearest_residue[1,'z']
    
          matlab_metadata = bind_rows(matlab_metadata, data.frame(Uniprot_ID = curr_uniprot_id,
                                                                  PDB_ID = curr_pdb_id,
                                                                  chain = curr_chain,
                                                                  fp_idr = fp_idr,
                                                                  fn_idr = fn_idr,
                                                                  ncpr_idr = ncpr_idr,
                                                                  tcpr_idr = tcpr_idr,
                                                                  idr_len = idr_len,
                                                                  fp_fd = fp_fd,
                                                                  fn_fd = fn_fd,
                                                                  ncpr_fd = ncpr_fd,
                                                                  tcpr_fd = tcpr_fd,
                                                                  fd_len = fd_len,
                                                                  tail = tail,
                                                                  idr_start = idr_start,
                                                                  fd_start_real = fd_actual_start,
                                                                  fd_start_reported = fd_start_reported,
                                                                  fd_start_pdb = fd_pdb_res_num,
                                                                  fd_coverage_percentage = fd_coverage_percentage))
          
        }
        else 
        {
          print(paste0(i,' excluded'))
        }
      }
    }
  }
}
  

matlab_metadata_unique = matlab_metadata %>% arrange(PDB_ID) %>% as.data.frame()

matlab_metadata_unique$FD_IDR_terminal_dist_real = abs(matlab_metadata_unique$idr_start - matlab_metadata_unique$fd_start_real)
matlab_metadata_unique$FD_start_error = abs(matlab_metadata_unique$fd_start_real - matlab_metadata_unique$fd_start_reported)
matlab_metadata_unique = matlab_metadata_unique %>% group_by(PDB_ID) %>% filter(FD_IDR_terminal_dist_real == min(FD_IDR_terminal_dist_real)) %>% as.data.frame()
matlab_metadata_unique = matlab_metadata_unique %>% group_by(PDB_ID) %>% filter(chain == min(chain)) %>% as.data.frame()

write.csv(matlab_metadata_unique, paste('./Processed_Data','matlab_metadata.csv',sep='/'),row.names=FALSE)


matlab_metadata_unique %>% filter(tail == 'N_terminal') %>% filter((idr_start - fd_start_real) > 5)

###tails removed for
####2vp7_A -- removed 333 - 339
####2wa0_A -- removed -14 to 0
####3fhr_A -- removed 33 to 38
####6r80_A -- removed 901 to 908



matlab_metadata_unique %>% filter(tail == 'C_terminal') %>% filter((idr_start - fd_start_real) < -5)
###tails removed for
####1eig_A -- removed last 5 residues
####1wym_A -- removed last 10 residues
####2do1_A -- removed last 7 residues



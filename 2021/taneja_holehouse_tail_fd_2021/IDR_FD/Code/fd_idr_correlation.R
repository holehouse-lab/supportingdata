library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gghighlight)

##gen idr fd data##

idr_fd_metadata = read.csv('../Processed_Data/matlab_metadata.csv')
fd_surface_charge = read.csv('../Processed_Data/fd_net_charge.csv') %>% select(PDB_ID, charge) %>% as.data.frame()
colnames(fd_surface_charge)[2] = 'fd_net_charge'
idr_fd_metadata = left_join(idr_fd_metadata, fd_surface_charge) %>% as.data.frame()

fd_surface_characteristic_dir = '../Processed_Data/FD_Surface_Characteristics'
fd_surface_characteristic_files = list.files(fd_surface_characteristic_dir)

fd_surface_characteristic_ref_dir = '../Processed_Data/FD_Surface_Characteristics_Random_Residue'
fd_surface_characteristic_ref_files = list.files(fd_surface_characteristic_ref_dir)

q_p_df = data.frame()
q_p_ref_distribution_df = data.frame()

for (file in fd_surface_characteristic_files)
{
  pdb_name = strsplit(file, '_')[[1]][1]
  
  if (grepl('q_f_radius.csv',file))
  {
    curr_pdb_q_f = read.csv(paste(fd_surface_characteristic_dir,file,sep='/'), header=FALSE)
    colnames(curr_pdb_q_f) = c('r', 'patchiness', 'q', 'N', 'total', 'count')
    curr_pdb_q_f$PDB_ID = pdb_name
    q_p_df = bind_rows(q_p_df, curr_pdb_q_f)
  }
}


for (file in fd_surface_characteristic_ref_files)
{
  pdb_name = strsplit(file, '_')[[1]][1]
  
  if (grepl('q_f_radius.csv',file))
  {
    curr_pdb_q_f = read.csv(paste(fd_surface_characteristic_ref_dir,file,sep='/'),header=FALSE)
    underscore_loc = gregexpr('_', file)[[1]]
    start_loc = underscore_loc[3]+1
    end_loc = underscore_loc[6]
    type_str = substr(file,start_loc,end_loc-1)
    colnames(curr_pdb_q_f) = c('r', 'patchiness', 'q', 'N', 'total', 'count')
    curr_pdb_q_f$PDB_ID = pdb_name
    curr_pdb_q_f$type = type_str
    q_p_ref_distribution_df = bind_rows(q_p_ref_distribution_df, curr_pdb_q_f)
  }
}

q_p_df$PDB_ID = as.factor(q_p_df$PDB_ID)
q_p_df$type = 'IDR_nearest'
q_p_ref_distribution_df$PDB_ID = as.factor(q_p_ref_distribution_df$PDB_ID)

q_p_df = left_join(q_p_df, idr_fd_metadata %>% select(-c('fp', 'fn', 'ncpr', 'tcpr')) %>% as.data.frame(), by = c('PDB_ID', 'type')) %>% as.data.frame()
q_p_ref_distribution_df = left_join(q_p_ref_distribution_df, idr_fd_metadata %>% select(-c('fp', 'fn', 'ncpr', 'tcpr')) %>% as.data.frame(), by = c('PDB_ID', 'type')) %>% as.data.frame()

q_p_df$FD_IDR_terminal_dist_real_net = (q_p_df$idr_start - q_p_df$fd_start_real)
q_p_ref_distribution_df$FD_IDR_terminal_dist_real_net = (q_p_ref_distribution_df$idr_start - q_p_ref_distribution_df$fd_start_real)

q_p_df$PDB_ID_w_type = paste(q_p_df$PDB_ID, q_p_df$type, q_p_df$tail, sep = '_')
q_p_ref_distribution_df$PDB_ID_w_type = paste(q_p_ref_distribution_df$PDB_ID, q_p_ref_distribution_df$type, q_p_ref_distribution_df$tail, sep = '_')





q_p_df_max_r = q_p_df %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()
q_p_ref_distribution_df_max_r = q_p_ref_distribution_df %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()

q_p_df_max_r = q_p_df_max_r[,c('r', 'patchiness', 'q', 'PDB_ID_w_type')]
colnames(q_p_df_max_r)[1:3] = c('r_max', 'patchiness_max', 'q_max')


q_p_ref_distribution_df_max_r = q_p_ref_distribution_df_max_r[,c('r', 'patchiness', 'q', 'PDB_ID_w_type')]
colnames(q_p_ref_distribution_df_max_r)[1:3] = c('r_max', 'patchiness_max', 'q_max')



##gen idr stats##

fcr_data = read.csv('../../IDR_only/Output_Data/Window/fcr_window.csv')
ncpr_data = read.csv('../../IDR_only/Output_Data/Window/ncpr_window.csv')
hydro_data = read.csv('../../IDR_only/Output_Data/Window/hydro_window.csv')
idr_len_data = read.csv('../../IDR_only/Output_Data/Window/idr_len_window.csv')
uniprot_data = read.csv('../../IDR_only/Output_Data/Window/uniprot_id_window.csv')

nterm_idr_window_data = cbind(as.numeric(fcr_data[,1]), as.numeric(ncpr_data[,1]), as.numeric(hydro_data[,1]), as.numeric(idr_len_data[,1]), uniprot_data[,1])
cterm_idr_window_data = cbind(as.numeric(fcr_data[,2]), as.numeric(ncpr_data[,2]), as.numeric(hydro_data[,2]), as.numeric(idr_len_data[,2]), uniprot_data[,2])

colnames(nterm_idr_window_data) = c('fcr', 'ncpr', 'hydropathy', 'idr_len', 'Uniprot_ID')
colnames(cterm_idr_window_data) = c('fcr', 'ncpr', 'hydropathy', 'idr_len', 'Uniprot_ID')

nterm_idr_window_data = as.data.frame(nterm_idr_window_data)
cterm_idr_window_data = as.data.frame(cterm_idr_window_data)

nterm_idr_window_data$fcr = as.numeric(nterm_idr_window_data$fcr)
nterm_idr_window_data$ncpr = as.numeric(nterm_idr_window_data$ncpr)
nterm_idr_window_data$hydropathy = as.numeric(nterm_idr_window_data$hydropathy)
nterm_idr_window_data$idr_len = as.numeric(nterm_idr_window_data$idr_len)

cterm_idr_window_data$fcr = as.numeric(cterm_idr_window_data$fcr)
cterm_idr_window_data$ncpr = as.numeric(cterm_idr_window_data$ncpr)
cterm_idr_window_data$hydropathy = as.numeric(cterm_idr_window_data$hydropathy)
cterm_idr_window_data$idr_len = as.numeric(cterm_idr_window_data$idr_len)



fd_df = read.csv(file = '../Raw_Data/uniprot_human_proteome_folded_domain.csv', header = TRUE)
colnames(fd_df)[1] = 'Uniprot_ID'
fd_df = fd_df[,-c(which(colnames(fd_df) == 'FD'))]
fd_df$PDB_ID = tolower(fd_df$PDB_ID)

nterm_idr_window_data = left_join(nterm_idr_window_data, fd_df[,c('Uniprot_ID', 'PDB_ID')])
cterm_idr_window_data = left_join(cterm_idr_window_data, fd_df[,c('Uniprot_ID', 'PDB_ID')])


q_p_df_nterm = q_p_df %>% filter(tail == 'N_terminal') %>% as.data.frame()
idr_fd_data_nterm = left_join(q_p_df_nterm, nterm_idr_window_data, by = 'PDB_ID')

q_p_df_cterm = q_p_df %>% filter(tail == 'C_terminal') %>% as.data.frame()
idr_fd_data_cterm = left_join(q_p_df_cterm, cterm_idr_window_data, by = 'PDB_ID')

all_idr_fd_data = rbind(idr_fd_data_nterm, idr_fd_data_cterm)



q_p_ref_df_nterm = q_p_ref_distribution_df %>% filter(tail == 'N_terminal') %>% as.data.frame()
idr_fd_data_nterm_ref = left_join(q_p_ref_df_nterm, nterm_idr_window_data, by = 'PDB_ID')

q_p_ref_df_cterm = q_p_ref_distribution_df %>% filter(tail == 'C_terminal') %>% as.data.frame()
idr_fd_data_cterm_ref = left_join(q_p_ref_df_cterm, cterm_idr_window_data, by = 'PDB_ID')

all_idr_fd_ref_data = rbind(idr_fd_data_nterm_ref, idr_fd_data_cterm_ref)


create_correlation_matrix  = function(data)
{
  correlation_df = data.frame()
  #create correlation matrix
  for (r_val in c(2,3,4,5,6,7,8,9,10))
  {
    for (len in c(2,3,4,5,6,7,8,9,10))
    {
      temp = data %>% filter(r == r_val) %>% filter(idr_len == len) %>% as.data.frame()
      fcr_correlation = cor.test(temp$q, temp$fcr, method = 'spearman')
      ncpr_correlation = cor.test(temp$q, temp$ncpr, method = 'spearman')
      hydropathy_correlation = cor.test(temp$q, temp$hydropathy, method = 'spearman')
      
      fcr_correlation_val = round(as.numeric(fcr_correlation$estimate),2)
      ncpr_correlation_val = round(as.numeric(ncpr_correlation$estimate),2)
      hydropathy_correlation_val = round(as.numeric(hydropathy_correlation$estimate),2)
      
      fcr_correlation_pval = round(as.numeric(fcr_correlation$p.value),2)
      ncpr_correlation_pval = round(as.numeric(ncpr_correlation$p.value),2)
      hydropathy_correlation_pval = round(as.numeric(hydropathy_correlation$p.value),2)
      
      correlation_df = bind_rows(correlation_df, data.frame(r = r_val, idr_len = len,
                                                            fcr_correlation = fcr_correlation_val,
                                                            ncpr_correlation = ncpr_correlation_val,
                                                            hydropathy_correlation = hydropathy_correlation_val,
                                                            fcr_correlation_pval = fcr_correlation_pval,
                                                            ncpr_correlation_pval = ncpr_correlation_pval,
                                                            hydropathy_correlation_pval = hydropathy_correlation_pval))
      
    }
  }
  
  return(correlation_df)
}
  
correlation_matrix = create_correlation_matrix(all_idr_fd_data)

correlation_matrix$r = as.factor(correlation_matrix$r )
correlation_matrix$idr_len = as.factor(correlation_matrix$idr_len)

plot_save_dir = '../Output/Surface_Charge_Characteristic_Plots/FD_IDR_Correlation/Local_Correlation_Matrix/N_CTerm'


print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = fcr_correlation)) + geom_text(aes(label = fcr_correlation)) + 
        scale_fill_gradient(limits = c(-.15,.15), name = expression(rho~'('~FCR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From FD:IDR Junction") +
        theme_classic() + 
          theme(text=element_text(size=20), legend.title=element_text(size=20)))
  ggsave(paste(plot_save_dir, paste0('fcr_correlation.tiff'), sep = '/'), dpi=300) 



print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = fcr_correlation_pval)) + geom_text(aes(label = fcr_correlation_pval)) + 
        scale_fill_gradient(limits = c(0,1), name = expression(rho[pval]~'('~FCR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From FD:IDR Junction") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('fcr_correlation_pval.tiff'), sep = '/'), dpi=300) 


print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = ncpr_correlation)) + geom_text(aes(label = ncpr_correlation)) + 
        scale_fill_gradient(limits = c(-.15,.15), name = expression(rho~'('~NCPR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From FD:IDR Junction") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('ncpr_correlation.tiff'), sep = '/'), dpi=300) 



print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = ncpr_correlation_pval)) + geom_text(aes(label = ncpr_correlation_pval)) + 
        scale_fill_gradient(limits = c(0,1), name = expression(rho[pval]~'('~NCPR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From FD:IDR Junction") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('ncpr_correlation_pval.tiff'), sep = '/'), dpi=300) 



print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = hydropathy_correlation)) + geom_text(aes(label = hydropathy_correlation)) + 
        scale_fill_gradient(limits = c(-.15,.15), name = expression(rho~'('~hydro[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From FD:IDR Junction") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('hydropathy_correlation.tiff'), sep = '/'), dpi=300) 




print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = hydropathy_correlation_pval)) + geom_text(aes(label = hydropathy_correlation_pval)) + 
        scale_fill_gradient(limits = c(0,1), name = expression(rho[pval]~'('~hydro[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From FD:IDR Junction") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('hydropathy_correlation_pval.tiff'), sep = '/'), dpi=300) 



correlation_matrix = create_correlation_matrix(all_idr_fd_ref_data)

correlation_matrix$r = as.factor(correlation_matrix$r )
correlation_matrix$idr_len = as.factor(correlation_matrix$idr_len)

plot_save_dir = './Output/Surface_Charge_Characteristic_Plots/FD_IDR_Correlation/Local_Correlation_Matrix/Random_Residue'


print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = fcr_correlation)) + geom_text(aes(label = fcr_correlation)) + 
        scale_fill_gradient(limits = c(-.15,.15), name = expression(rho~'('~FCR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From Random Residue") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('fcr_correlation.tiff'), sep = '/'), dpi=300) 



print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = fcr_correlation_pval)) + geom_text(aes(label = fcr_correlation_pval)) + 
        scale_fill_gradient(limits = c(0,1), name = expression(rho[pval]~'('~FCR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From Random Residue") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('fcr_correlation_pval.tiff'), sep = '/'), dpi=300) 



print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = ncpr_correlation)) + geom_text(aes(label = ncpr_correlation)) + 
        scale_fill_gradient(limits = c(-.15,.15), name = expression(rho~'('~NCPR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From Random Residue") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('ncpr_correlation.tiff'), sep = '/'), dpi=300) 



print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = ncpr_correlation_pval)) + geom_text(aes(label = ncpr_correlation_pval)) + 
        scale_fill_gradient(limits = c(0,1), name = expression(rho[pval]~'('~NCPR[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From Random Residue") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('ncpr_correlation_pval.tiff'), sep = '/'), dpi=300) 


print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = hydropathy_correlation)) + geom_text(aes(label = hydropathy_correlation)) + 
        scale_fill_gradient(limits = c(-.15,.15), name = expression(rho~'('~hydro[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From Random Residue") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('hydropathy_correlation.tiff'), sep = '/'), dpi=300) 


print(ggplot(correlation_matrix, aes(r, idr_len)) + 
        geom_tile(aes(fill = hydropathy_correlation_pval)) + geom_text(aes(label = hydropathy_correlation_pval)) + 
        scale_fill_gradient(limits = c(0,1), name = expression(rho[pval]~'('~hydro[1:x]~','~hat(phi)(r)~')'), low = "white", high = "red") +
        labs(x = "Number of IDR Residues", y = "Distance From Random Residue") +
        theme_classic() + 
        theme(text=element_text(size=20), legend.title=element_text(size=20)))
ggsave(paste(plot_save_dir, paste0('hydropathy_correlation_pval.tiff'), sep = '/'), dpi=300) 



##global correlation##



gen_corr_scatterplot = function(df)
{
  save_dir = paste0('./Output/Surface_Charge_Characteristic_Plots/FD_IDR_Correlation/','Global_Correlation')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  p = (ggscatter(df, x = "q", y = "ncpr",
                 add = "reg.line",  
                 xlab = expression(hat(phi)), ylab = 'NCPR',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman") + font("xy.text", size = 25) + font("xlab", size = 25) + font("ylab", size = 25) 
  print(p)
  ggsave(paste(save_dir, paste0('q_ncpr_corr.tiff'), sep = '/'))
  
  p = (ggscatter(df, x = "q", y = "fcr",
                 add = "reg.line", 
                 xlab = expression(hat(phi)), ylab = 'FCR',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman") + font("xy.text", size = 25) + font("xlab", size = 25) + font("ylab", size = 25) 
  print(p)
  ggsave(paste(save_dir, paste0('q_fcr_corr.tiff'), sep = '/'),dpi = 300)
  
  p = (ggscatter(df, x = "fd_net_charge", y = "ncpr",
                 add = "reg.line", 
                 xlab = 'FD Net Charge', ylab = 'NCPR',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman") + font("xy.text", size = 25) + font("xlab", size = 25) + font("ylab", size = 25) 
  print(p)
  
  ggsave(paste(save_dir, paste0('fd_net_charge_ncpr_corr.tiff'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(df, x = "fd_net_charge", y = "fcr",
                 add = "reg.line",  
                 xlab = 'FD Net Charge', ylab = 'FCR',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman") + font("xy.text", size = 25) + font("xlab", size = 25) + font("ylab", size = 25) 
  print(p)
  ggsave(paste(save_dir, paste0('fd_net_charge_fcr_corr.tiff'), sep = '/'),dpi = 300)
  
}

all_idr_fd_data_max_r = all_idr_fd_data %>% group_by(PDB_ID_w_type) %>% filter(r == max(r)) %>% filter(idr_len == max(idr_len)) %>% as.data.frame()
gen_corr_scatterplot(all_idr_fd_data_max_r)
  


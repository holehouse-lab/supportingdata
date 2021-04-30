library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gghighlight)


idr_fd_metadata = read.csv('../Processed_Data/matlab_metadata.csv')
fd_surface_charge = read.csv('../Output/Surface_Charge_Characteristic_Plots/Net_Charge/net_charge.csv') %>% select(PDB_ID, charge)
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

q_p_df = left_join(q_p_df, idr_fd_metadata, by = c('PDB_ID', 'type')) %>% as.data.frame()
q_p_ref_distribution_df = left_join(q_p_ref_distribution_df, idr_fd_metadata, by = c('PDB_ID', 'type')) %>% as.data.frame()

q_p_df$FD_IDR_terminal_dist_real_net = (q_p_df$idr_start - q_p_df$fd_start_real)
q_p_ref_distribution_df$FD_IDR_terminal_dist_real_net = (q_p_ref_distribution_df$idr_start - q_p_ref_distribution_df$fd_start_real)

q_p_df$PDB_ID_w_type = paste(q_p_df$PDB_ID, q_p_df$type, q_p_df$tail, sep = '_')
q_p_ref_distribution_df$PDB_ID_w_type = paste(q_p_ref_distribution_df$PDB_ID, q_p_ref_distribution_df$type, q_p_ref_distribution_df$tail, sep = '_')


gen_surface_charge_characteristic_plots = function(df, df_type, y_val,y_val_str)
{
  get_plot_title = function(df){
    
    num_proteins = length(unique(df$PDB_ID))
    num_instances = length(unique(df$PDB_ID_w_type))
    
    if (df_type == 'Reference_Distribution') {
      title = paste0("Num Proteins = ",num_proteins, ' ; ', "Num Instances = ", num_instances)
    } else {
      title = paste0("Num Proteins = ",num_proteins)
    }
    
    return(title)
  }
  
  if (df_type == 'Reference_Distribution') {
    x_val_str = 'Distance From Random Residue (r)'
  } else {
    x_val_str = 'Distance From FD:IDR Junction (r)'
  }
  
  
  
  if (y_val == 'patchiness') {
    ylim_vec = c(0,1)
  } else {
    ylim_vec = c(-30,30)
  }
  xlim_vec = c(0,100)
  
  parent_dir = paste('../Output/Surface_Charge_Characteristic_Plots/Q_P_Characterization',df_type, y_val_str, sep = '/')
  save_dir = paste0(parent_dir,'/All')
  
  df_min_r = df %>% group_by(PDB_ID_w_type) %>% slice(which.min(r)) %>% as.data.frame()
  df_min_r = df_min_r[,c('r', 'patchiness', 'q', 'PDB_ID_w_type')]
  colnames(df_min_r)[1:3] = c('r_min', 'patchiness_min', 'q_min')
  df = left_join(df, df_min_r, by = 'PDB_ID_w_type') %>% as.data.frame()
  df$qmin_lt0 = as.integer(df$q_min < 0)
  df$qmin_lt0 = factor(df$qmin_lt0)
  df$pmin_lt50 = as.integer(df$patchiness_min < .5)
  df$pmin_lt50 = factor(df$pmin_lt50)
  
  if (y_val == 'patchiness') {
    color_str = 'pmin_lt50'
  } else {
    color_str = 'qmin_lt0'
  }
  
  
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  y_val_mean = paste('mean', y_val, sep = '_')
  y_val_se_upper = paste('se_upper', y_val, sep = '_')
  y_val_se_lower = paste('se_lower', y_val, sep = '_')
  
  mean_df_f_r = df %>% group_by(r) %>% summarise(mean_patchiness = mean(patchiness), mean_q = mean(q)) %>% as.data.frame()
  mean_df_max_r = df %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(r), mean_patchiness = mean(patchiness), mean_q = mean(q), sd_q = sd(q), sd_patchiness = sd(patchiness), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_q = mean_q - 1.96*sd_q, 
                                           se_upper_q = mean_q + 1.96*sd_q,
                                           se_lower_patchiness = mean_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = mean_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df)
  
  print(ggplot(df, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  
  print(ggplot(df, aes_string(x='r', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'r', y = y_val_mean), color = 'red') 
          #geom_point(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean), color = 'white') +
          #geom_errorbar(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean, ymin=y_val_se_lower, ymax=y_val_se_upper), color = 'white')
        )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
  
  df_coverage_gte_.75 = df %>% filter(fd_coverage_percentage > .75) %>% as.data.frame()
  
  mean_df_f_r = df_coverage_gte_.75 %>% group_by(r) %>% summarise(mean_patchiness = mean(patchiness), mean_q = mean(q)) %>% as.data.frame()
  mean_df_max_r = df_coverage_gte_.75 %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(r), mean_patchiness = mean(patchiness), mean_q = mean(q), sd_q = sd(q), sd_patchiness = sd(patchiness), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_q = mean_q - 1.96*sd_q, 
                                           se_upper_q = mean_q + 1.96*sd_q,
                                           se_lower_patchiness = mean_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = mean_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  save_dir = paste0(parent_dir,'/Coverage_Gte_X')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df_coverage_gte_.75)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='r', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'r', y = y_val_mean), color = 'red') 
          #geom_point(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean), color = 'white') +
          #geom_errorbar(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean, ymin=y_val_se_lower, ymax=y_val_se_upper), color = 'white')
  )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
  
  relevant_pdb= df %>% group_by(PDB_ID_w_type) %>% slice(which.min(r)) %>% filter(q < 0) %>% as.data.frame()
  df_qrmin_neg = df %>% filter(PDB_ID_w_type %in% relevant_pdb$PDB_ID_w_type) %>% as.data.frame()
  
  mean_df_f_r = df_qrmin_neg %>% group_by(r) %>% summarise(mean_patchiness = mean(patchiness), mean_q = mean(q)) %>% as.data.frame()
  mean_df_max_r = df_qrmin_neg %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(r), mean_patchiness = mean(patchiness), mean_q = mean(q), sd_q = sd(q), sd_patchiness = sd(patchiness), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_q = mean_q - 1.96*sd_q, 
                                           se_upper_q = mean_q + 1.96*sd_q,
                                           se_lower_patchiness = mean_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = mean_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  
  save_dir = paste0(parent_dir,'/QRmin_Neg')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  
  plot_title = get_plot_title(df_qrmin_neg)

  print(ggplot(df_qrmin_neg, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_qrmin_neg, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_qrmin_neg, aes_string(x='r', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'r', y = y_val_mean), color = 'red') 
          #geom_point(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean), color = 'white') +
          #geom_errorbar(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean, ymin=y_val_se_lower, ymax=y_val_se_upper), color = 'white')
  )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
  
  
  relevant_pdb= df %>% group_by(PDB_ID_w_type) %>% slice(which.min(r)) %>% filter(q > 0) %>% as.data.frame()
  df_qrmin_pos = df %>% filter(PDB_ID_w_type %in% relevant_pdb$PDB_ID_w_type) %>% as.data.frame()
  
  
  mean_df_f_r = df_qrmin_pos %>% group_by(r) %>% summarise(mean_patchiness = mean(patchiness), mean_q = mean(q)) %>% as.data.frame()
  mean_df_max_r = df_qrmin_pos %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(r), mean_patchiness = mean(patchiness), mean_q = mean(q), sd_q = sd(q), sd_patchiness = sd(patchiness), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_q = mean_q - 1.96*sd_q, 
                                           se_upper_q = mean_q + 1.96*sd_q,
                                           se_lower_patchiness = mean_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = mean_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  
  save_dir = paste0(parent_dir,'/QRmin_Pos')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df_qrmin_pos)
  
  print(ggplot(df_qrmin_pos, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_qrmin_pos, aes_string(x='r', y=y_val, group='PDB_ID_w_type', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_qrmin_pos, aes_string(x='r', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'r', y = y_val_mean), color = 'red')
          #geom_point(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean), color = 'white') +
          #geom_errorbar(data = mean_df_max_r, aes_string(x = 'mean_r', y = y_val_mean, ymin=y_val_se_lower, ymax=y_val_se_upper), color = 'white')
  )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
}

gen_surface_charge_characteristic_plots(q_p_df, 'IDR_Origin', 'q','Mean Potential')
gen_surface_charge_characteristic_plots(q_p_df, 'IDR_Origin', 'patchiness','Patchiness')


gen_surface_charge_characteristic_plots(q_p_ref_distribution_df, 'Reference_Distribution', 'patchiness','Patchiness')
gen_surface_charge_characteristic_plots(q_p_ref_distribution_df, 'Reference_Distribution', 'q','Mean Potential')



q_p_df_max_r = q_p_df %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()
q_p_ref_distribution_df_max_r = q_p_ref_distribution_df %>% group_by(PDB_ID_w_type) %>% slice(which.max(r)) %>% as.data.frame()

q_p_df_max_r %>% summarise(mean_q = mean(q), mean_patchiness = mean(patchiness)) %>% as.data.frame()
q_p_df_max_r %>% summarise(quantile_q = quantile(q, c(0.25, 0.5, 0.75)), quantile_patchiness = quantile(patchiness, c(0.25, 0.5, 0.75))) %>% as.data.frame()

q_p_df_max_r = q_p_df_max_r[,c('r', 'patchiness', 'q', 'PDB_ID_w_type')]
colnames(q_p_df_max_r)[1:3] = c('r_max', 'patchiness_max', 'q_max')


q_p_ref_distribution_df_max_r = q_p_ref_distribution_df_max_r[,c('r', 'patchiness', 'q', 'PDB_ID_w_type')]
colnames(q_p_ref_distribution_df_max_r)[1:3] = c('r_max', 'patchiness_max', 'q_max')



##rm histogram##
plot_save_dir = '../Output/Surface_Charge_Characteristic_Plots/Rm'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

plot_name = 'rm.tiff'
print(ggplot(q_p_df_max_r, aes(x = r_max)) +
        geom_histogram() +
        geom_vline(aes(xintercept = mean(r_max)), linetype = "dashed") +
        labs(x = expression(r[m]), title=paste0("N = ", nrow(q_p_df_max_r)), x ="NCPR", y = "Count") +
        theme_classic() + 
        theme(text=element_text(size=20)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)
##

#rX q vs rm q#
plot_save_dir = '../Output/Surface_Charge_Characteristic_Plots/RXq_v_Rmq'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)


q_p_df_w_max_r = left_join(q_p_df, q_p_df_max_r, by = 'PDB_ID_w_type') %>% as.data.frame()

q_p_df_w_max_r_at_r5 = q_p_df_w_max_r %>% filter(r == 5) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 
plot_name = 'qdiff_ecdf_r5.tiff'
print(ggplot(q_p_df_w_max_r_at_r5, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(5)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r5)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)



q_p_df_w_max_r_at_r10 = q_p_df_w_max_r %>% filter(r == 10) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 
plot_name = 'qdiff_hist_r10.tiff'
print(ggplot(q_p_df_w_max_r_at_r10, aes(x = potential_diff)) +
        geom_histogram() +
        geom_vline(aes(xintercept = mean(potential_diff)), linetype = "dashed") +
        labs(x = expression(~hat(phi)(10)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r10)), y = "Count") +
        theme_classic() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)


plot_name = 'qdiff_ecdf_r10.tiff'
print(ggplot(q_p_df_w_max_r_at_r10, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(10)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r10)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)



q_p_df_w_max_r_at_r15 = q_p_df_w_max_r %>% filter(r == 15) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 

plot_name = 'qdiff_hist_r15.tiff'
print(ggplot(q_p_df_w_max_r_at_r15, aes(x = potential_diff)) +
        geom_histogram() +
        geom_vline(aes(xintercept = mean(potential_diff)), linetype = "dashed") +
        labs(x = expression(~hat(phi)(15)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r15)), y = "Count") +
        theme_classic() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)


plot_name = 'qdiff_ecdf_r15.tiff'
print(ggplot(q_p_df_w_max_r_at_r15, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(15)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r15)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)



q_p_df_w_max_r_at_r20 = q_p_df_w_max_r %>% filter(r == 20) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 

plot_name = 'qdiff_hist_r20.tiff'
print(ggplot(q_p_df_w_max_r_at_r20, aes(x = potential_diff)) +
        geom_histogram() +
        geom_vline(aes(xintercept = mean(potential_diff)), linetype = "dashed") +
        labs(x = expression(~hat(phi)(20)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r20)), y = "Count") +
        theme_classic() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)


plot_name = 'qdiff_ecdf_r20.tiff'
print(ggplot(q_p_df_w_max_r_at_r20, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(20)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_df_w_max_r_at_r20)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)


####rX q vs rm q from reference point###

q_p_ref_distribution_df_w_max_r = left_join(q_p_ref_distribution_df, q_p_ref_distribution_df_max_r, by = 'PDB_ID_w_type') %>% as.data.frame()
plot_save_dir = '../Output/Surface_Charge_Characteristic_Plots/RXq_v_Rmq/Reference_Point'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)


q_p_ref_distribution_df_w_max_r_at_r5 = q_p_ref_distribution_df_w_max_r %>% filter(r == 5) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 


plot_name = 'qdiff_ecdf_r5.tiff'
print(ggplot(q_p_ref_distribution_df_w_max_r_at_r5, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(5)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_ref_distribution_df_w_max_r_at_r5)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)




q_p_ref_distribution_df_w_max_r_at_r10 = q_p_ref_distribution_df_w_max_r %>% filter(r == 10) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 


plot_name = 'qdiff_ecdf_r10.tiff'
print(ggplot(q_p_ref_distribution_df_w_max_r_at_r10, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(10)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_ref_distribution_df_w_max_r_at_r10)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)



q_p_ref_distribution_df_w_max_r_at_r15 = q_p_ref_distribution_df_w_max_r %>% filter(r == 15) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 

plot_name = 'qdiff_ecdf_r15.tiff'
print(ggplot(q_p_ref_distribution_df_w_max_r_at_r15, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(15)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_ref_distribution_df_w_max_r_at_r15)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)




q_p_ref_distribution_df_w_max_r_at_r20 = q_p_ref_distribution_df_w_max_r %>% filter(r == 20) %>% mutate(potential_diff = q-q_max, abs_potential_diff = abs(q-q_max)) %>% as.data.frame() 

plot_name = 'qdiff_ecdf_r20.tiff'
print(ggplot(q_p_ref_distribution_df_w_max_r_at_r20, aes(x = potential_diff)) +
        stat_ecdf(size = 2) +
        scale_x_continuous(breaks = seq(-10,10,2)) +
        coord_cartesian(xlim = c(-10, 10)) + 
        labs(x = expression(~hat(phi)(20)~-~hat(phi)(r[m])), title=paste0("N = ", nrow(q_p_ref_distribution_df_w_max_r_at_r20)), y = "Percentage") +
        theme_bw() + 
        theme(text=element_text(size=25)))
ggsave(paste(plot_save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)





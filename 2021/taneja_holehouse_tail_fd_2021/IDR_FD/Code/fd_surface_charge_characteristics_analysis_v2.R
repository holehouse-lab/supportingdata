library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gghighlight)
#library(plyr)
#detach("package:plyr", unload=TRUE)

idr_fd_metadata = read.csv('./Processed_Data/matlab_metadata.csv')
phi_patchiness_input = read.csv('./Processed_Data/phi_patchiness_input.csv')

phi_patchiness_input = phi_patchiness_input %>% select(Uniprot_ID, PDB_ID, chain, reference_point, residue) %>% as.data.frame()

####add FD net charge and dipole moment 
#fd_surface_charge = read.csv('./Output/Surface_Charge_Characteristic_Plots/Net_Charge/net_charge.csv') %>% select(PDB_ID, charge)
#colnames(fd_surface_charge)[2] = 'fd_net_charge'
#idr_fd_metadata = left_join(idr_fd_metadata, fd_surface_charge) %>% as.data.frame()

fd_surface_characteristic_dir = './Processed_Data/Phi_Patchiness_Data'
fd_surface_characteristic_files = list.files(fd_surface_characteristic_dir)

phi_patchiness_f_r = data.frame()
phi_patchiness_f_surf_dist = data.frame()

phi_patchiness_patch_loc = data.frame()
phi_patchiness_patch_dist_distribution = data.frame()

for (i in 1:length(fd_surface_characteristic_files))
{
  if (i %% 100 == 0) { print (i) }
  
  file = fd_surface_characteristic_files[i]
  curr_pdb = strsplit(file, '_')[[1]][1]
  curr_chain = strsplit(file, '_')[[1]][2]
  curr_residue = strsplit(file, '_')[[1]][4]
  
  tmp = phi_patchiness_input %>% filter(PDB_ID == curr_pdb) %>% filter(chain == curr_chain) %>% filter(residue == curr_residue) %>% as.data.frame()
  ref_point = tmp[1,'reference_point']
  
  if (grepl('patchiness_q_f_radius.csv',file))
  {
    curr_phi_patchiness_f_r = read.csv(paste(fd_surface_characteristic_dir,file,sep='/'))
    curr_phi_patchiness_f_r = curr_phi_patchiness_f_r %>% select(-num_points) %>% as.data.frame()
    curr_phi_patchiness_f_r$PDB_ID = curr_pdb
    curr_phi_patchiness_f_r$chain = curr_chain
    curr_phi_patchiness_f_r$ref_point = ref_point
    curr_phi_patchiness_f_r$residue_ref = curr_residue
    phi_patchiness_f_r = bind_rows(phi_patchiness_f_r, curr_phi_patchiness_f_r)
  }
  
  if (grepl('patchiness_q_f_surface_distance.csv',file))
  {
    curr_phi_patchiness_f_surface_dist = read.csv(paste(fd_surface_characteristic_dir,file,sep='/'))
    curr_phi_patchiness_f_surface_dist = curr_phi_patchiness_f_surface_dist %>% select(-num_points) %>% as.data.frame()
    curr_phi_patchiness_f_surface_dist$PDB_ID = curr_pdb
    curr_phi_patchiness_f_surface_dist$chain = curr_chain
    curr_phi_patchiness_f_surface_dist$ref_point = ref_point
    curr_phi_patchiness_f_surface_dist$residue_ref = curr_residue
    phi_patchiness_f_surf_dist = bind_rows(phi_patchiness_f_surf_dist, curr_phi_patchiness_f_surface_dist)
  }
  
  if (grepl('patch_loc.csv',file))
  {
    curr_patch_loc = read.csv(paste(fd_surface_characteristic_dir,file,sep='/'))
    if ('patch_centroid_idr_dist' %in% colnames(curr_patch_loc)) {
      print(file)
    }
    curr_patch_loc$PDB_ID = curr_pdb
    curr_patch_loc$chain = curr_chain
    curr_patch_loc$ref_point = ref_point
    curr_patch_loc$residue_ref = curr_residue
    phi_patchiness_patch_loc = bind_rows(phi_patchiness_patch_loc, curr_patch_loc)
  }
  if (grepl('patch_dist',file))
  {
    curr_patch_dist_distribution = read.csv(paste(fd_surface_characteristic_dir,file,sep='/'))
    curr_patch_dist_distribution$PDB_ID = curr_pdb
    curr_patch_dist_distribution$chain = curr_chain
    curr_patch_dist_distribution$ref_point = ref_point
    curr_patch_dist_distribution$residue_ref = curr_residue
    phi_patchiness_patch_dist_distribution = bind_rows(phi_patchiness_patch_dist_distribution, curr_patch_dist_distribution)
  }
  
}


phi_patchiness_f_r = left_join(phi_patchiness_f_r, idr_fd_metadata %>% select(Uniprot_ID, PDB_ID, fp_idr, fn_idr, ncpr_idr, tcpr_idr, idr_len, fp_fd, fn_fd, ncpr_fd, tcpr_fd, fd_len, fd_coverage_percentage) %>% as.data.frame())
phi_patchiness_f_surf_dist = left_join(phi_patchiness_f_surf_dist, idr_fd_metadata %>% select(Uniprot_ID, PDB_ID, fp_idr, fn_idr, ncpr_idr, tcpr_idr, idr_len, fp_fd, fn_fd, ncpr_fd, tcpr_fd, fd_len, fd_coverage_percentage) %>% as.data.frame())
phi_patchiness_patch_loc = left_join(phi_patchiness_patch_loc, idr_fd_metadata %>% select(Uniprot_ID, PDB_ID, fp_idr, fn_idr, ncpr_idr, tcpr_idr, idr_len, fp_fd, fn_fd, ncpr_fd, tcpr_fd, fd_len, fd_coverage_percentage) %>% as.data.frame())
phi_patchiness_patch_dist_distribution = left_join(phi_patchiness_patch_dist_distribution, idr_fd_metadata %>% select(Uniprot_ID, PDB_ID, fp_idr, fn_idr, ncpr_idr, tcpr_idr, idr_len, fp_fd, fn_fd, ncpr_fd, tcpr_fd, fd_len, fd_coverage_percentage) %>% as.data.frame())

phi_patchiness_patch_dist_distribution$ideal_chain_ee_dist = sqrt(phi_patchiness_patch_dist_distribution$idr_len)*3.5
patch_accessibility_df = phi_patchiness_patch_dist_distribution %>% group_by(PDB_ID, residue_ref, patch_centroid_x, patch_centroid_y, patch_centroid_z) %>% summarise(patch_accessibility_percentage = sum(distance_to_sampled_patch_node <  ideal_chain_ee_dist)/n()) %>% as.data.frame()

phi_patchiness_f_r$PDB_ID_w_residue_ref = paste(phi_patchiness_f_r$PDB_ID, phi_patchiness_f_r$residue_ref, sep = '_')
phi_patchiness_f_surf_dist$PDB_ID_w_residue_ref = paste(phi_patchiness_f_surf_dist$PDB_ID, phi_patchiness_f_surf_dist$residue_ref, sep = '_')
phi_patchiness_patch_loc$PDB_ID_w_residue_ref = paste(phi_patchiness_patch_loc$PDB_ID, phi_patchiness_patch_loc$residue_ref, sep = '_')

phi_patchiness_patch_data_all = left_join(phi_patchiness_patch_loc, patch_accessibility_df) %>% filter(patch_size_percentage >= .05) %>% as.data.frame()


gen_surface_charge_characteristic_plots_f_r = function(df, df_type, y_val,y_val_str)
{
  
  if (df_type == 'Reference_Distribution') {
    df = df %>% filter(ref_point == 'Random') %>% as.data.frame()
  } else {
    df = df %>% filter(ref_point != 'Random') %>% as.data.frame()
  }
  
  get_plot_title = function(df){
    
    num_proteins = length(unique(df$PDB_ID))
    num_instances = length(unique(df$PDB_ID_w_residue_ref))
    
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

  
  if (y_val == 'mean_patchiness') {
    ylim_vec = c(0,1)
  } else if (y_val == 'mean_phi') {
    ylim_vec = c(-30,30)
  }
  
  xlim_vec = c(0,100)
  
  parent_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Phi_Patchiness_F_R',df_type, y_val_str, sep = '/')
  save_dir = paste0(parent_dir,'/All')
  
  
  df_min_r = df %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.min(radius)) %>% as.data.frame()
  df_min_r = df_min_r[,c('radius', 'mean_patchiness', 'mean_phi', 'PDB_ID_w_residue_ref')]
  colnames(df_min_r)[1:3] = c('radius_min', 'mean_patchiness_min', 'mean_phi_min')
  df = left_join(df, df_min_r, by = 'PDB_ID_w_residue_ref') %>% as.data.frame()
  df$qmin_lt0 = as.integer(df$mean_phi_min < 0)
  df$qmin_lt0 = factor(df$qmin_lt0)
  df$pmin_lt50 = as.integer(df$mean_patchiness_min < .5)
  df$pmin_lt50 = factor(df$pmin_lt50)
  
  if (y_val == 'mean_patchiness') {
    color_str = 'pmin_lt50'
  } else if (y_val == 'mean_phi') { 
    color_str = 'qmin_lt0'
  }
  
  
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  mean_df_f_r = df %>% group_by(radius) %>% summarise(mean_patchiness = mean(mean_patchiness), mean_phi = mean(mean_phi)) %>% as.data.frame()
  
  mean_df_max_r = df %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(radius)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(radius), avg_patchiness = mean(mean_patchiness), avg_phi = mean(mean_phi), sd_patchiness = sd(mean_patchiness), sd_phi = sd(mean_phi), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_phi = avg_phi - 1.96*sd_phi, 
                                           se_upper_phi = avg_phi + 1.96*sd_phi,
                                           se_lower_patchiness = avg_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = avg_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df)
  
  print(ggplot(df, aes_string(x='radius', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df, aes_string(x='radius', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  
  print(ggplot(df, aes_string(x='radius', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'radius', y = y_val), color = 'red') 
        )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
  
  df_coverage_gte_.75 = df %>% filter(fd_coverage_percentage > .75) %>% as.data.frame()
  
  mean_df_f_r = df_coverage_gte_.75 %>% group_by(radius) %>% summarise(mean_patchiness = mean(mean_patchiness), mean_phi = mean(mean_phi)) %>% as.data.frame()
  
  mean_df_max_r = df_coverage_gte_.75 %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(radius)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(radius), avg_patchiness = mean(mean_patchiness), avg_phi = mean(mean_phi), sd_patchiness = sd(mean_patchiness), sd_phi = sd(mean_phi), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_phi = avg_phi - 1.96*sd_phi, 
                                           se_upper_phi = avg_phi + 1.96*sd_phi,
                                           se_lower_patchiness = avg_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = avg_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  save_dir = paste0(parent_dir,'/Coverage_Gte_X')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df_coverage_gte_.75)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='radius', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='radius', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='radius', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'radius', y = y_val), color = 'red') 
  )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
}



gen_surface_charge_characteristic_plots_f_surface_dist = function(df, df_type, y_val,y_val_str)
{
  
  if (df_type == 'Reference_Distribution') {
    df = df %>% filter(ref_point == 'Random') %>% as.data.frame()
  } else {
    df = df %>% filter(ref_point != 'Random') %>% as.data.frame()
  }
  
  get_plot_title = function(df){
    
    num_proteins = length(unique(df$PDB_ID))
    num_instances = length(unique(df$PDB_ID_w_residue_ref))
    
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
  
  
  if (y_val == 'mean_patchiness') {
    ylim_vec = c(0,1)
  } else if (y_val == 'mean_phi') {
    ylim_vec = c(-30,30)
  }
  
  xlim_vec = c(0,100)
  
  parent_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Phi_Patchiness_F_Surface_Dist',df_type, y_val_str, sep = '/')
  save_dir = paste0(parent_dir,'/All')

  
  df_min_r = df %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.min(surface_dist)) %>% as.data.frame()
  df_min_r = df_min_r[,c('surface_dist', 'mean_patchiness', 'mean_phi', 'PDB_ID_w_residue_ref')]
  colnames(df_min_r)[1:3] = c('surface_dist_min', 'mean_patchiness_min', 'mean_phi_min')
  df = left_join(df, df_min_r, by = 'PDB_ID_w_residue_ref') %>% as.data.frame()
  df$qmin_lt0 = as.integer(df$mean_phi_min < 0)
  df$qmin_lt0 = factor(df$qmin_lt0)
  df$pmin_lt50 = as.integer(df$mean_patchiness_min < .5)
  df$pmin_lt50 = factor(df$pmin_lt50)
  
  if (y_val == 'mean_patchiness') {
    color_str = 'pmin_lt50'
  } else if (y_val == 'mean_phi') { 
    color_str = 'qmin_lt0'
  }
  
  
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  mean_df_f_r = df %>% group_by(surface_dist) %>% summarise(mean_patchiness = mean(mean_patchiness), mean_phi = mean(mean_phi)) %>% as.data.frame()
  
  mean_df_max_r = df %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(surface_dist)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(surface_dist), avg_patchiness = mean(mean_patchiness), avg_phi = mean(mean_phi), sd_patchiness = sd(mean_patchiness), sd_phi = sd(mean_phi), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_phi = avg_phi - 1.96*sd_phi, 
                                           se_upper_phi = avg_phi + 1.96*sd_phi,
                                           se_lower_patchiness = avg_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = avg_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df)
  
  print(ggplot(df, aes_string(x='surface_dist', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df, aes_string(x='surface_dist', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  
  print(ggplot(df, aes_string(x='surface_dist', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'surface_dist', y = y_val), color = 'red') 
  )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
  
  df_coverage_gte_.75 = df %>% filter(fd_coverage_percentage > .75) %>% as.data.frame()
  
  mean_df_f_r = df_coverage_gte_.75 %>% group_by(surface_dist) %>% summarise(mean_patchiness = mean(mean_patchiness), mean_phi = mean(mean_phi)) %>% as.data.frame()
  
  mean_df_max_r = df_coverage_gte_.75 %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(surface_dist)) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% summarise(mean_r = mean(surface_dist), avg_patchiness = mean(mean_patchiness), avg_phi = mean(mean_phi), sd_patchiness = sd(mean_patchiness), sd_phi = sd(mean_phi), n = n()) %>% as.data.frame()
  mean_df_max_r = mean_df_max_r %>% mutate(se_lower_phi = avg_phi - 1.96*sd_phi, 
                                           se_upper_phi = avg_phi + 1.96*sd_phi,
                                           se_lower_patchiness = avg_patchiness - 1.96*sd_patchiness,
                                           se_upper_patchiness = avg_patchiness + 1.96*sd_patchiness) %>% as.data.frame()
  
  save_dir = paste0(parent_dir,'/Coverage_Gte_X')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(mean_df_max_r, paste(save_dir, 'mean_vals.csv', sep = '/'), row.names=FALSE)
  
  plot_title = get_plot_title(df_coverage_gte_.75)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='surface_dist', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = c(.5,.2), plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='surface_dist', y=y_val, group='PDB_ID_w_residue_ref', color = color_str)) +
          geom_line(aes_string(color=color_str), alpha=0.1)+
          geom_point(aes_string(color=color_str), alpha=0.1, size = .2) +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme_classic() + geom_smooth(aes_string(group=color_str), se=FALSE, size =2) + 
          theme(legend.position = 'none', plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title))
  ggsave(paste(save_dir, paste0(y_val,'_f_r_scatter_wo_legend.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df_coverage_gte_.75, aes_string(x='surface_dist', y=y_val)) +
          geom_bin2d() +
          scale_fill_continuous(type = "viridis") +
          theme_classic() +
          coord_cartesian(ylim=ylim_vec, xlim = xlim_vec) + 
          theme(legend.position = "none", plot.title = element_text(size=15), text=element_text(size=25)) +
          labs(x = x_val_str, y = y_val_str, title=plot_title) +
          geom_point(data = mean_df_f_r, aes_string(x = 'surface_dist', y = y_val), color = 'red') 
  )
  ggsave(paste(save_dir, paste0(y_val,'_f_r_heatmap.tiff'), sep = '/'),dpi = 300)
  
}





gen_surface_charge_characteristic_plots_f_r(phi_patchiness_f_r, 'IDR_Origin', 'mean_phi','Mean Potential')
gen_surface_charge_characteristic_plots_f_surface_dist(phi_patchiness_f_surf_dist, 'IDR_Origin', 'mean_phi','Mean Potential')

gen_surface_charge_characteristic_plots_f_r(phi_patchiness_f_r, 'IDR_Origin', 'mean_patchiness','Patchiness')
gen_surface_charge_characteristic_plots_f_surface_dist(phi_patchiness_f_surf_dist, 'IDR_Origin', 'mean_patchiness','Patchiness')

gen_surface_charge_characteristic_plots_f_r(phi_patchiness_f_r, 'Reference_Distribution', 'mean_phi','Mean Potential')
gen_surface_charge_characteristic_plots_f_surface_dist(phi_patchiness_f_surf_dist, 'Reference_Distribution', 'mean_phi','Mean Potential')

gen_surface_charge_characteristic_plots_f_r(phi_patchiness_f_r, 'Reference_Distribution', 'mean_patchiness','Patchiness')
gen_surface_charge_characteristic_plots_f_surface_dist(phi_patchiness_f_surf_dist, 'Reference_Distribution', 'mean_patchiness','Patchiness')




patchiness_loc_plots_individual_protein = function(df, df_type)
{
  if (df_type == 'Reference_Distribution') {
    df = df %>% filter(ref_point == 'Random') %>% as.data.frame()
  } else {
    df = df %>% filter(ref_point != 'Random') %>% as.data.frame()
  }
  
  if (df_type == 'Reference_Distribution') {
    x_val_str = 'Distance From Random Residue (r)'
  } else {
    x_val_str = 'Distance From FD:IDR Junction (r)'
  }
  
  
  
  unique_instances = unique(df$PDB_ID_w_residue_ref) 
  

  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Patch_Loc_Individual_Protein/Euclidean_Distance',df_type, sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Patch_Loc_Individual_Protein/Surface_Distance',df_type, sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  df$patch_size_percentage_scaled = df$patch_size_percentage*100
  
  for (curr_PDB_ID_w_residue_ref in unique_instances)
  {
    curr_df = df %>% filter(PDB_ID_w_residue_ref == curr_PDB_ID_w_residue_ref)
    
    if (nrow(curr_df) == 1)
    {
      min_ylim = curr_df$patch_mean_q - 5
      max_ylim = curr_df$patch_mean_q + 5
    }
    else
    {
      if (min(curr_df$patch_mean_q) < 0)
      {
        min_ylim = min(curr_df$patch_mean_q)*2
      }
      else
      {
        min_ylim = min(curr_df$patch_mean_q)/2
      }
      
      if (max(curr_df$patch_mean_q) > 0)
      {
        max_ylim = max(curr_df$patch_mean_q)*2
      }
      else
      {
        max_ylim = max(curr_df$patch_mean_q)/2
      }
    }
    
    ref_point = curr_df[1,'ref_point']
    residue_ref = curr_df[1, 'residue_ref']
    chain = curr_df[1, 'chain']
    
    plot_name = paste(curr_df[1,'PDB_ID_w_residue_ref'], chain, ref_point, sep = '_')
    
    max_radius = phi_patchiness_f_r %>% filter(PDB_ID_w_residue_ref == curr_PDB_ID_w_residue_ref) %>% slice(which.max(radius)) %>% select(radius) %>% as.data.frame()
    max_radius = max_radius[1,]
    xlim_vec_radius = c(0,max_radius)
    
    max_surface_dist = phi_patchiness_f_surf_dist %>% filter(PDB_ID_w_residue_ref == curr_PDB_ID_w_residue_ref) %>% slice(which.max(surface_dist)) %>% select(surface_dist) %>% as.data.frame()
    max_surface_dist = max_surface_dist[1,]
    xlim_vec_surface_dist = c(0,max_surface_dist)
    
    title = paste0('PDB: ', curr_df[1,'PDB_ID'], ', Ref. Point: ', ref_point)

    
    save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Patch_Loc_Individual_Protein/Surface_Distance',df_type, sep = '/')
    print(ggplot(curr_df, aes(x=patch_centroid_idr_surf_dist, y=patch_mean_q)) +
            geom_point(aes(color=patch_accessibility_percentage, size = patch_size_percentage_scaled), alpha = 0.5) +
            scale_size_identity(breaks = c(5, 10, 20, 30, 40, 50), guide = "legend") +
            coord_cartesian(xlim = xlim_vec_surface_dist, ylim = c(min_ylim, max_ylim)) + 
            theme_classic() + 
            scale_color_gradientn(limits = c(0,1), colours = rainbow(5)) +
            theme(text=element_text(size=25), legend.title = element_blank(), legend.key.size = unit(.5, "cm")) +
            labs(x = 'Surface Distance To Patch Centroid', y = expression(~hat(phi)), title=title))
    ggsave(paste(save_dir, paste0(plot_name, '.tiff'), sep = '/'),dpi = 300)
  }
    
}

patchiness_loc_plots_individual_protein(phi_patchiness_patch_data_all, 'IDR_Origin')
patchiness_loc_plots_individual_protein(phi_patchiness_patch_data_all, 'Reference_Distribution')


patchiness_loc_plots_all_protein = function(df, df_type)
{
  if (df_type == 'Reference_Distribution') {
    df = df %>% filter(ref_point == 'Random') %>% as.data.frame()
  } else {
    df = df %>% filter(ref_point != 'Random') %>% as.data.frame()
  }
  
  df = df %>% filter(patch_size_percentage >= .1) %>% as.data.frame()
  
  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Patch_Loc_All_Protein/Euclidean_Distance',df_type, sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  print(ggplot(df, aes(PDB_ID, patch_centroid_idr_euclid_dist))+
      geom_point(aes(color=patch_mean_q), alpha = 0.5) +
      scale_color_gradient2(midpoint=0, low="red", mid="white", high="blue") +
      theme_linedraw()+ 
      labs(x = 'Protein', y = 'Euclidean Distance To Patch Centroid', title = 'Mean Potential') +
      theme(legend.title = element_blank(), text=element_text(size=10), axis.text.x = element_blank())+
      coord_fixed())
  ggsave(paste(save_dir, paste0('color_f_phi', '.tiff'), sep = '/'),dpi = 300)
  

  print(ggplot(df, aes(PDB_ID, patch_centroid_idr_euclid_dist))+
          geom_point(aes(size=patch_size_percentage), alpha = 0.5) +
          scale_radius() +
          theme_linedraw()+ 
          labs(x = 'Protein', y = 'Euclidean Distance To Patch Centroid', title = 'Patch Size %') +
          theme(legend.title = element_blank(), text=element_text(size=10), axis.text.x = element_blank())+
          coord_fixed())
  ggsave(paste(save_dir, paste0('color_f_size', '.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df, aes(PDB_ID, patch_centroid_idr_euclid_dist))+
          geom_point(aes(size=patch_accessibility_percentage), alpha = 0.5) +
          scale_radius() +
          theme_linedraw()+ 
          labs(x = 'Protein', y = 'Euclidean Distance To Patch Centroid', title = 'Patch Accessibility %') +
          theme(legend.title = element_blank(), text=element_text(size=10), axis.text.x = element_blank())+
          coord_fixed())
  ggsave(paste(save_dir, paste0('color_f_accessibility', '.tiff'), sep = '/'),dpi = 300)
  
  ################
  
  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Patch_Loc_All_Protein/Surface_Distance',df_type, sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  print(ggplot(df, aes(PDB_ID, patch_centroid_idr_surf_dist))+
          geom_point(aes(color=patch_mean_q), alpha = 0.5) +
          scale_color_gradient2(midpoint=0, low="red", mid="white", high="blue") +
          theme_linedraw()+ 
          labs(x = 'Protein', y = 'Surface Distance To Patch Centroid', title = 'Mean Potential') +
          theme(legend.title = element_blank(), text=element_text(size=10), axis.text.x = element_blank())+
          coord_fixed())
  ggsave(paste(save_dir, paste0('color_f_phi', '.tiff'), sep = '/'),dpi = 300)
  
  
  print(ggplot(df, aes(PDB_ID, patch_centroid_idr_surf_dist))+
          geom_point(aes(size=patch_size_percentage), alpha = 0.5) +
          scale_radius() +
          theme_linedraw()+ 
          labs(x = 'Protein', y = 'Surface Distance To Patch Centroid', title = 'Patch Size %') +
          theme(legend.title = element_blank(), text=element_text(size=10), axis.text.x = element_blank())+
          coord_fixed())
  ggsave(paste(save_dir, paste0('color_f_size', '.tiff'), sep = '/'),dpi = 300)
  
  print(ggplot(df, aes(PDB_ID, patch_centroid_idr_surf_dist))+
          geom_point(aes(size=patch_accessibility_percentage), alpha = 0.5) +
          scale_radius() +
          theme_linedraw()+ 
          labs(x = 'Protein', y = 'Surface Distance To Patch Centroid', title = 'Patch Accessibility %') +
          theme(legend.title = element_blank(), text=element_text(size=10), axis.text.x = element_blank())+
          coord_fixed())
  ggsave(paste(save_dir, paste0('color_f_accessibility', '.tiff'), sep = '/'),dpi = 300)
  
}


patchiness_loc_plots_all_protein(phi_patchiness_patch_data_all, 'IDR_Origin')



gen_histogram_plots = function()
{
  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Histograms', sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  phi_patchiness_max_surf_dist = phi_patchiness_f_surf_dist %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(surface_dist)) %>% as.data.frame()
  phi_patchiness_max_surf_dist = phi_patchiness_max_surf_dist %>% filter(ref_point != 'Random') %>% as.data.frame()
  
  print(mean(phi_patchiness_max_surf_dist$mean_phi))
  
  print(ggplot(phi_patchiness_max_surf_dist, aes(mean_phi))+
          geom_histogram(position="identity", color = "black", fill = "grey") +
          theme_classic() +
          labs(x = expression(~hat(phi)), y = 'Count') +
          geom_vline(aes(xintercept=mean(mean_phi)), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('mean_potential', '.tiff'), sep = '/'),dpi = 300)
  
  print(mean(phi_patchiness_max_surf_dist$mean_patchiness))
  
  print(ggplot(phi_patchiness_max_surf_dist, aes(mean_patchiness))+
          geom_histogram(position="identity", color = "black", fill = "grey") +
          theme_classic() +
          labs(x = 'Patchiness', y = 'Count') +
          geom_vline(aes(xintercept=mean(mean_patchiness)), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('mean_patchiness', '.tiff'), sep = '/'),dpi = 300)
  
  print(mean(phi_patchiness_max_surf_dist$ncpr_fd))

  print(ggplot(phi_patchiness_max_surf_dist, aes(ncpr_fd))+
          geom_histogram(position="identity", color = "black", fill = "grey") +
          theme_classic() +
          labs(x = 'FD NCPR', y = 'Count') +
          geom_vline(aes(xintercept=mean(ncpr_fd)), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('fd_ncpr', '.tiff'), sep = '/'),dpi = 300)
  
  
  phi_patchiness_patch_data_all_idrorign = phi_patchiness_patch_data_all %>% filter(ref_point != 'Random') %>% as.data.frame()
  
  print(mean(phi_patchiness_patch_data_all_idrorign$patch_size_percentage))
  
  print(ggplot(phi_patchiness_patch_data_all_idrorign, aes(patch_size_percentage))+
          geom_histogram(position="identity", color = "black", fill = "grey") +
          theme_classic() +
          labs(x = 'Patch Size Fraction', y = 'Count') +
          geom_vline(aes(xintercept=mean(patch_size_percentage)), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_size_percentage', '.tiff'), sep = '/'),dpi = 300)
  
  
  phi_patchiness_patch_data_all_idrorign[phi_patchiness_patch_data_all_idrorign$patch_mean_q > 0, 'patch_type'] = 'Positive'
  phi_patchiness_patch_data_all_idrorign[phi_patchiness_patch_data_all_idrorign$patch_mean_q < 0, 'patch_type'] = 'Negative'
  
  mu_patch_q = phi_patchiness_patch_data_all_idrorign %>% group_by(patch_type) %>% summarise(grp.mean=mean(patch_mean_q)) %>% as.data.frame()
  
  colnames(phi_patchiness_patch_data_all_idrorign)[which(colnames(phi_patchiness_patch_data_all_idrorign) == 'patch_type')] = 'Patch Type'
  
  print(ggplot(phi_patchiness_patch_data_all_idrorign, aes(x = patch_mean_q, fill = `Patch Type`, color = `Patch Type`))+
          geom_histogram(position="identity", alpha=.5) +
          theme_classic() +
          labs(x = expression(~hat(phi)), y = 'Count') +
          geom_vline(data=mu_patch_q, aes(xintercept=grp.mean, color = patch_type), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_mean_potential', '.tiff'), sep = '/'),dpi = 300)
  
  
  
  
  phi_patchiness_patch_data_all_for_plotting = phi_patchiness_patch_data_all
  phi_patchiness_patch_data_all_for_plotting[phi_patchiness_patch_data_all$ref_point == 'N_terminal','ref_point'] = 'N/C Terminal'
  phi_patchiness_patch_data_all_for_plotting[phi_patchiness_patch_data_all$ref_point == 'C_terminal','ref_point'] = 'N/C Terminal'
  
  mu_accessibility = phi_patchiness_patch_data_all_for_plotting %>% group_by(ref_point) %>% summarise(grp.mean=mean(patch_accessibility_percentage)) %>% as.data.frame()
  mu_patch_dist = phi_patchiness_patch_data_all_for_plotting %>% group_by(ref_point) %>% summarise(grp.mean=mean(patch_centroid_idr_surf_dist)) %>% as.data.frame()
  
  colnames(phi_patchiness_patch_data_all_for_plotting)[which(colnames(phi_patchiness_patch_data_all_for_plotting) == 'ref_point')] = 'FD:IDR Junction'
  
  print(ggplot(phi_patchiness_patch_data_all_for_plotting, aes(x = patch_accessibility_percentage, fill = `FD:IDR Junction`, color = `FD:IDR Junction`))+
          geom_histogram(position="identity", alpha=.5) +
          theme_classic() +
          labs(x = 'Patch Accessibility Fraction', y = 'Count') +
          geom_vline(data=mu_accessibility, aes(xintercept=grp.mean, color = ref_point), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_accessibility', '.tiff'), sep = '/'),dpi = 300)
  

  print(ggplot(phi_patchiness_patch_data_all_for_plotting, aes(x = patch_centroid_idr_surf_dist, fill = `FD:IDR Junction`, color = `FD:IDR Junction`))+
          geom_histogram(position="identity", alpha=.5) +
          theme_classic() +
          labs(x = '', y = 'Count') +
          geom_vline(data=mu_patch_dist, aes(xintercept=grp.mean, color = ref_point), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_centroid_surf_dist', '.tiff'), sep = '/'),dpi = 300)
  
  
  summary_stats_per_protein_mean_size = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(size = mean(patch_size_percentage)) %>% mutate(stat = 'Mean') %>% as.data.frame() 
  summary_stats_per_protein_max_size = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(size = max(patch_size_percentage)) %>% mutate(stat = 'Max') %>% as.data.frame() 
  summary_stats_per_protein_min_size = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(size = min(patch_size_percentage)) %>% mutate(stat = 'Min') %>% as.data.frame() 
  summary_stats_per_protein_size = bind_rows(summary_stats_per_protein_mean_size,
                                             summary_stats_per_protein_max_size,
                                             summary_stats_per_protein_min_size)
  
  mu_size = summary_stats_per_protein_size %>% group_by(stat) %>% summarise(grp.mean=mean(size)) %>% as.data.frame()
  
  print(ggplot(summary_stats_per_protein_size, aes(x = size, fill = stat, color = stat))+
          geom_density(alpha=.5) +
          theme_classic() +
          labs(x = 'Patch Size Fraction', y = 'Density') +
          geom_vline(data=mu_size, aes(xintercept=grp.mean, color = stat), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_size_per_protein', '.tiff'), sep = '/'),dpi = 300)
  
  summary_stats_per_protein_mean_dist = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(dist = mean(patch_centroid_idr_surf_dist)) %>% mutate(stat = 'Mean') %>% as.data.frame() 
  summary_stats_per_protein_max_dist = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(dist = max(patch_centroid_idr_surf_dist)) %>% mutate(stat = 'Max') %>% as.data.frame() 
  summary_stats_per_protein_min_dist = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(dist = min(patch_centroid_idr_surf_dist)) %>% mutate(stat = 'Min') %>% as.data.frame() 
  summary_stats_per_protein_dist = bind_rows(summary_stats_per_protein_mean_dist,
                                             summary_stats_per_protein_max_dist,
                                             summary_stats_per_protein_min_dist)
  
  nrow(summary_stats_per_protein_mean_dist)
  
  mu_dist = summary_stats_per_protein_dist %>% group_by(stat) %>% summarise(grp.mean=mean(dist)) %>% as.data.frame()
  
  print(ggplot(summary_stats_per_protein_dist, aes(x = dist, fill = stat, color = stat))+
          geom_density(alpha=.5) +
          theme_classic() +
          labs(x = 'Surface Distance To Patch Centroid', y = 'Density') +
          geom_vline(data=mu_dist, aes(xintercept=grp.mean, color = stat), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_dist_per_protein', '.tiff'), sep = '/'),dpi = 300)
  
  summary_stats_per_protein_mean_accessibility = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(accessibility = mean(patch_accessibility_percentage)) %>% mutate(stat = 'Mean') %>% as.data.frame() 
  summary_stats_per_protein_max_accessibility = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(accessibility = max(patch_accessibility_percentage)) %>% mutate(stat = 'Max') %>% as.data.frame() 
  summary_stats_per_protein_min_accessibility = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(accessibility = min(patch_accessibility_percentage)) %>% mutate(stat = 'Min') %>% as.data.frame() 
  summary_stats_per_protein_accessibility = bind_rows(summary_stats_per_protein_mean_accessibility,
                                                      summary_stats_per_protein_max_accessibility,
                                                      summary_stats_per_protein_min_accessibility)
  
  nrow(summary_stats_per_protein_mean_accessibility)
  
  
  mu_accessibility = summary_stats_per_protein_accessibility %>% group_by(stat) %>% summarise(grp.mean=mean(accessibility)) %>% as.data.frame()
  
  print(ggplot(summary_stats_per_protein_accessibility, aes(x = accessibility, fill = stat, color = stat))+
          geom_density(alpha=.5) +
          theme_classic() +
          labs(x = 'Patch Accessibility Fraction', y = 'Density') +
          geom_vline(data=mu_accessibility, aes(xintercept=grp.mean, color = stat), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_accessibility_per_protein', '.tiff'), sep = '/'),dpi = 300)
  
  
  summary_stats_per_protein_mean_potential = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(potential = mean(patch_mean_q)) %>% mutate(stat = 'Mean') %>% as.data.frame() 
  summary_stats_per_protein_max_potential = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(potential = max(patch_mean_q)) %>% mutate(stat = 'Max') %>% as.data.frame() 
  summary_stats_per_protein_min_potential = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(potential = min(patch_mean_q)) %>% mutate(stat = 'Min') %>% as.data.frame() 
  summary_stats_per_protein_potential = bind_rows(summary_stats_per_protein_mean_potential,
                                                  summary_stats_per_protein_max_potential,
                                                  summary_stats_per_protein_min_potential)
  
  mu_potential = summary_stats_per_protein_potential %>% group_by(stat) %>% summarise(grp.mean=mean(potential)) %>% as.data.frame()
  
  print(ggplot(summary_stats_per_protein_potential, aes(x = potential, fill = stat, color = stat))+
          geom_density(alpha=.5) +
          theme_classic() +
          labs(x = expression(hat(phi)), y = 'Density') +
          geom_vline(data=mu_potential, aes(xintercept=grp.mean, color = stat), linetype="dashed") +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_mean_potential_per_protein', '.tiff'), sep = '/'),dpi = 300)
  
  
  
 
}

gen_histogram_plots()


gen_patch_accessibility_correlation = function()
{
  tmp = phi_patchiness_patch_data_all %>% filter(ref_point != 'Random') %>% as.data.frame()
  tmp$combined_metric = tmp$idr_len*tmp$patch_size_percentage/tmp$patch_centroid_idr_surf_dist
  cor(tmp$patch_accessibility_percentage, tmp$idr_len, method = 'spearman')
  cor(tmp$patch_accessibility_percentage, tmp$patch_size_percentage, method = 'spearman')
  cor(tmp$patch_accessibility_percentage, tmp$patch_centroid_idr_surf_dist, method = 'spearman')
  cor(tmp$patch_accessibility_percentage, tmp$combined_metric, method = 'spearman')
  
  
  save_dir = paste0('./Output/Surface_Charge_Characteristic_Plots/Patch_Accessibility_Correlation')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  ##filtering just for plotting purposes due to 2 instances where IDR len > 1000 -- correlation coefficient is the same; 
  p = (ggscatter(tmp %>% filter(idr_len <= 1000) %>% as.data.frame(), x = "patch_accessibility_percentage", y = "idr_len",
                 add = "reg.line",  
                 xlab = 'Patch Accessibility Fraction', ylab = 'IDR Length',
                 ylim = c(0,1000),
                 font.label = c(20),
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman")
  print(p)
  ggsave(paste(save_dir, paste0('accessibility_idrlen_corr.png'), sep = '/'),dpi = 300)
  

  p = (ggscatter(tmp, x = "patch_accessibility_percentage", y = "patch_size_percentage",
                 add = "reg.line",  
                 xlab = 'Patch Accessibility Fraction', ylab = 'Patch Size %',
                 font.label = c(20),
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman")
  print(p)
  ggsave(paste(save_dir, paste0('accessibility_patchsize.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(tmp, x = "patch_accessibility_percentage", y = "patch_centroid_idr_surf_dist",
                 add = "reg.line",  
                 xlab = 'Patch Accessibility Fraction', ylab = 'Patch Centroid Distance Along Surface',
                 font.label = c(20),
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman")
  print(p)
  ggsave(paste(save_dir, paste0('accessibility_patchsurfdist.png'), sep = '/'),dpi = 300)
  
  
  ##filtering just for plotting purposes due to 1 instance where combined_metric > 9 -- correlation coefficient is the same; 
  p = (ggscatter(tmp %>% filter(combined_metric <= 9) %>% as.data.frame(), x = "patch_accessibility_percentage", y = "combined_metric",
                 add = "reg.line",  
                 xlab = 'Patch Accessibility Fraction', ylab = 'Length*Size/Distance',
                 font.label = c(20),
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE))
  p = p + stat_cor(method = "spearman")
  print(p)
  ggsave(paste(save_dir, paste0('accessibility_combined_lensizedist.png'), sep = '/'),dpi = 300)
  
  
}

gen_patch_accessibility_correlation()




gen_num_patches_plots = function(all_pos_or_neg, ylim_max)
{
  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Num_Patches', all_pos_or_neg, sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  phi_patchiness_patch_data_all_idrorign = phi_patchiness_patch_data_all %>% filter(ref_point != 'Random') %>% as.data.frame()
  
  phi_patchiness_patch_data_all_idrorign[phi_patchiness_patch_data_all_idrorign$patch_mean_q > 0, 'patch_type'] = 'Positive'
  phi_patchiness_patch_data_all_idrorign[phi_patchiness_patch_data_all_idrorign$patch_mean_q < 0, 'patch_type'] = 'Negative'
  
  phi_patchiness_patch_data_all_idrorign_pos = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_type == 'Positive') %>% as.data.frame()
  phi_patchiness_patch_data_all_idrorign_neg = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_type == 'Negative') %>% as.data.frame()
  
  if (all_pos_or_neg == 'All') {
    df = phi_patchiness_patch_data_all_idrorign 
  } else if (all_pos_or_neg == 'Positive') {
    df = phi_patchiness_patch_data_all_idrorign_pos
  } else if (all_pos_or_neg == 'Negative') {
    df = phi_patchiness_patch_data_all_idrorign_neg
  }
  
  num_patches_total_df = df %>% group_by(PDB_ID_w_residue_ref) %>% summarise(num_patches = n()) %>% mutate(type = 'All') %>% as.data.frame() 
  
  num_patches_sizegt10_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_size_percentage >= .1) %>% summarise(num_patches = n()) %>% mutate(type = '>= 10%') %>% as.data.frame()
  num_patches_sizegt20_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_size_percentage >= .2) %>% summarise(num_patches = n()) %>% mutate(type = '>= 20%') %>% as.data.frame()
  num_patches_sizegt30_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_size_percentage >= .3) %>% summarise(num_patches = n()) %>% mutate(type = '>= 30%') %>% as.data.frame()
  num_patches_sizegt40_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_size_percentage >= .4) %>% summarise(num_patches = n()) %>% mutate(type = '>= 40%') %>% as.data.frame()
  
  num_patches_distgt20_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_centroid_idr_surf_dist >= 20) %>% summarise(num_patches = n()) %>% mutate(type = '>= 20') %>% as.data.frame()
  num_patches_distgt40_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_centroid_idr_surf_dist >= 40) %>% summarise(num_patches = n()) %>% mutate(type = '>= 40') %>% as.data.frame()
  num_patches_distgt60_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_centroid_idr_surf_dist >= 60) %>% summarise(num_patches = n()) %>% mutate(type = '>= 60') %>% as.data.frame()

  num_patches_accgt25_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_accessibility_percentage >= .25) %>% summarise(num_patches = n()) %>% mutate(type = '>= 25%') %>%as.data.frame()
  num_patches_accgt50_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_accessibility_percentage >= .50) %>% summarise(num_patches = n()) %>% mutate(type = '>= 50%') %>% as.data.frame()
  num_patches_accgt75_df = df %>% group_by(PDB_ID_w_residue_ref) %>% filter(patch_accessibility_percentage >= .75) %>% summarise(num_patches = n()) %>% mutate(type = '>= 75%') %>% as.data.frame()
  
  num_patches_size_combined_df = bind_rows(num_patches_total_df,
                                      num_patches_sizegt10_df,
                                      num_patches_sizegt20_df,
                                      num_patches_sizegt30_df,
                                      num_patches_sizegt40_df)
  
  num_patches_dist_combined_df = bind_rows(num_patches_total_df,
                                           num_patches_distgt20_df,
                                           num_patches_distgt40_df,
                                           num_patches_distgt60_df)
  
  num_patches_acc_combined_df = bind_rows(num_patches_total_df, 
                                          num_patches_accgt25_df, 
                                          num_patches_accgt50_df, 
                                          num_patches_accgt75_df)
  
  colnames(num_patches_size_combined_df)[3] = 'Patch Size'
  num_patches_size_combined_df$num_patches = as.factor(num_patches_size_combined_df$num_patches)
  
  print(ggplot(num_patches_size_combined_df, aes(x = num_patches, fill = `Patch Size`, color = `Patch Size`))+
          geom_bar(stat="count",alpha=.5) +
          theme_classic() +
          coord_cartesian(ylim = c(0, ylim_max)) +
          labs(x = 'Number of Patches', y = 'Count') +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_size_count', '.tiff'), sep = '/'),dpi = 300)
 
  colnames(num_patches_dist_combined_df)[3] = 'Patch Distance'
  num_patches_dist_combined_df$num_patches = as.factor(num_patches_dist_combined_df$num_patches)
  
  print(ggplot(num_patches_dist_combined_df, aes(x = num_patches, fill = `Patch Distance`, color = `Patch Distance`))+
          geom_bar(stat="count",alpha=.5) +
          theme_classic() +
          coord_cartesian(ylim = c(0, ylim_max)) +
          labs(x = 'Number of Patches', y = 'Count') +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_dist_count', '.tiff'), sep = '/'),dpi = 300)
  
  colnames(num_patches_acc_combined_df)[3] = 'Patch Accessibility'
  num_patches_acc_combined_df$num_patches = as.factor(num_patches_acc_combined_df$num_patches)
  
  print(ggplot(num_patches_acc_combined_df, aes(x = num_patches, fill = `Patch Accessibility`, color = `Patch Accessibility`))+
          geom_bar(stat="count",alpha=.5) +
          theme_classic() +
          coord_cartesian(ylim = c(0, ylim_max)) +
          labs(x = 'Number of Patches', y = 'Count') +
          theme(text=element_text(size=20)))
  ggsave(paste(save_dir, paste0('patch_accesibility_count', '.tiff'), sep = '/'),dpi = 300)
  
    
}

gen_num_patches_plots('All', 300)
gen_num_patches_plots('Positive', 400)
gen_num_patches_plots('Negative', 500)



protein_size_plots = function()
{
  save_dir = paste('./Output/Surface_Charge_Characteristic_Plots/Protein_Dimensions', sep = '/')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  phi_patchiness_max_surf_dist = phi_patchiness_f_surf_dist %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(surface_dist)) %>% as.data.frame()
  phi_patchiness_max_surf_dist = phi_patchiness_max_surf_dist %>% filter(ref_point != 'Random') %>% as.data.frame()
  
  phi_patchiness_max_r = phi_patchiness_f_r %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(radius)) %>% as.data.frame()
  phi_patchiness_max_r = phi_patchiness_max_r %>% filter(ref_point != 'Random') %>% as.data.frame()
  
  print(mean(phi_patchiness_max_r$radius))
  
  plot_name = 'Radius_Max.tiff'
  print(ggplot(phi_patchiness_max_r, aes(x = radius)) +
          geom_histogram(position="identity", color = "black", fill = "grey") +
          geom_vline(aes(xintercept = mean(radius)), linetype = "dashed") +
          labs(x = 'Maximum Radius', title=paste0("N = ", nrow(phi_patchiness_max_r)), y = "Count") +
          theme_classic() + 
          theme(text=element_text(size=25)))
  ggsave(paste(save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)
  
  print(mean(phi_patchiness_max_surf_dist$surface_dist))
  
  plot_name = 'Surface_Distance_Max.tiff'
  print(ggplot(phi_patchiness_max_surf_dist, aes(x = surface_dist)) +
          geom_histogram(position="identity", color = "black", fill = "grey") +
          geom_vline(aes(xintercept = mean(surface_dist)), linetype = "dashed") +
          labs(x = 'Maximum Surface Distance', title=paste0("N = ", nrow(phi_patchiness_max_surf_dist)), y = "Count") +
          theme_classic() + 
          theme(text=element_text(size=25)))
  ggsave(paste(save_dir, plot_name, sep = '/'), width = 8, height = 6, dpi=300)
  
}

protein_size_plots()




gen_fd_idr_corr_scatterplot = function()
{
  save_dir = paste0('./Output/Surface_Charge_Characteristic_Plots/FD_IDR_Correlation')
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  phi_patchiness_max_surf_dist = phi_patchiness_f_surf_dist %>% group_by(PDB_ID_w_residue_ref) %>% slice(which.max(surface_dist)) %>% as.data.frame()
  phi_patchiness_max_surf_dist = phi_patchiness_max_surf_dist %>% filter(ref_point != 'Random') %>% as.data.frame()
  
  
  p = (ggscatter(phi_patchiness_max_surf_dist, x = "mean_phi", y = "ncpr_fd",
            add = "reg.line",  
            xlab = expression(hat(phi)), ylab = 'FD NCPR',
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('mean_phi_ncpr_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(phi_patchiness_max_surf_dist, x = "ncpr_idr", y = "ncpr_fd",
                  add = "reg.line", 
                  xlab = 'IDR NCPR', ylab = 'FD NCPR',
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('ncpr_idr_ncpr_fd_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(phi_patchiness_max_surf_dist, x = "tcpr_idr", y = "tcpr_fd",
                 add = "reg.line", 
                 xlab = 'IDR FCR', ylab = 'FD FCR',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('tcpr_idr_tcpr_fd_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(phi_patchiness_max_surf_dist, x = "ncpr_idr", y = "mean_phi",
                 add = "reg.line", 
                 xlab = 'IDR NCPR', ylab = expression(hat(phi)),
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('ncpr_idr_mean_phi_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(phi_patchiness_max_surf_dist, x = "tcpr_idr", y = "mean_phi",
                 add = "reg.line", 
                 xlab = 'IDR FCR', ylab = expression(hat(phi)),
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('tcpr_idr_mean_phi_corr.png'), sep = '/'),dpi = 300)
  
  
  ####patchiness correlation####
  
  phi_patchiness_patch_data_all_idrorign = phi_patchiness_patch_data_all %>% filter(ref_point != 'Random') %>% as.data.frame()
  phi_patchiness_patch_data_all_idrorign$size_times_potential = phi_patchiness_patch_data_all_idrorign$patch_size_percentage*phi_patchiness_patch_data_all_idrorign$patch_mean_q
  phi_patchiness_patch_data_all_idrorign$acc_times_potential = phi_patchiness_patch_data_all_idrorign$patch_accessibility_percentage*phi_patchiness_patch_data_all_idrorign$patch_mean_q
  
  summary_stats_per_protein = phi_patchiness_patch_data_all_idrorign %>% group_by(PDB_ID) %>% summarise(total_size = sum(patch_size_percentage), 
                                                                                                        total_accessibility = sum(patch_accessibility_percentage), 
                                                                                                        net_potential = sum(patch_mean_q),
                                                                                                        net_size_potential = sum(size_times_potential),
                                                                                                        net_acc_potential = sum(acc_times_potential)) %>% as.data.frame() 
  
  summary_stats_per_protein = summary_stats_per_protein %>% left_join(phi_patchiness_max_surf_dist[,c('PDB_ID', 'ncpr_idr', 'tcpr_idr')]) %>% as.data.frame()
  
  p = (ggscatter(summary_stats_per_protein, x = "ncpr_idr", y = "net_size_potential",
                 add = "reg.line", 
                 xlab = 'IDR NCPR', ylab = 'Net Size*Potential',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('ncpr_idr_sizepotential_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(summary_stats_per_protein, x = "ncpr_idr", y = "net_acc_potential",
                 add = "reg.line", 
                 xlab = 'IDR NCPR', ylab = 'Net Accessibility*Potential',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('ncpr_idr_accpotential_corr.png'), sep = '/'),dpi = 300)
  
  
  
  p = (ggscatter(summary_stats_per_protein, x = "ncpr_idr", y = "total_size",
                 add = "reg.line", 
                 xlab = 'IDR NCPR', ylab = 'Total Patch Size',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('ncpr_idr_patchsize_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(summary_stats_per_protein, x = "tcpr_idr", y = "net_size_potential",
                 add = "reg.line", 
                 xlab = 'IDR FCR', ylab = 'Net Size*Potential',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('tcpr_idr_sizepotential_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(summary_stats_per_protein, x = "tcpr_idr", y = "net_acc_potential",
                 add = "reg.line", 
                 xlab = 'IDR FCR', ylab = 'Net Accessibility*Potential',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('tcpr_idr_accpotential_corr.png'), sep = '/'),dpi = 300)
  
  
  p = (ggscatter(summary_stats_per_protein, x = "tcpr_idr", y = "total_size",
                 add = "reg.line", 
                 xlab = 'IDR FCR', ylab = 'Total Patch Size',
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 font.label = c(20),
                 conf.int = TRUE))
  p = p + stat_cor(method = "pearson")
  print(p)
  ggsave(paste(save_dir, paste0('tcpr_idr_patchsize_corr.png'), sep = '/'),dpi = 300)
  
}

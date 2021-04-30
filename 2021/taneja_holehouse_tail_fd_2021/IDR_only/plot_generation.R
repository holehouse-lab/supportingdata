library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix) #requires installation of 'coin' packages as well

fcr_data = read.csv('./Output_Data/fcr.csv')
ncpr_data = read.csv('./Output_Data/ncpr.csv')
fp_data = read.csv('./Output_Data/fp.csv')
fn_data = read.csv('./Output_Data/fn.csv')
idr_len_data = read.csv('./Output_Data/idr_len.csv')
idr_len_sequence_len_data = read.csv('./Output_Data/idr_len_sequence_len.csv')
kappa_data = read.csv('./Output_Data/kappa.csv')
phospho_data = read.csv('./Output_Data/phospho.csv')


fcr_long_df = fcr_data %>% gather(idr_type, fcr, all:C_w_pdb) %>% filter(!is.na(fcr)) %>% as.data.frame()
ncpr_long_df = ncpr_data %>% gather(idr_type, ncpr, all:C_w_pdb) %>% filter(!is.na(ncpr)) %>% as.data.frame()
fp_long_df = fp_data %>% gather(idr_type, fp, all:C_w_pdb) %>% filter(!is.na(fp)) %>% as.data.frame()
fn_long_df = fn_data %>% gather(idr_type, fn, all:C_w_pdb) %>% filter(!is.na(fn)) %>% as.data.frame()
idr_len_long_df = idr_len_data %>% gather(idr_type, idr_len, all:C_w_pdb) %>% filter(!is.na(idr_len)) %>% as.data.frame()
idr_len_sequence_len_long_df = idr_len_sequence_len_data %>% gather(idr_type, idr_len_sequence_len, all:C_w_pdb) %>% filter(!is.na(idr_len_sequence_len)) %>% as.data.frame()
kappa_long_df = kappa_data %>% gather(idr_type, kappa, all:C_w_pdb) %>% filter(!is.na(kappa)) %>% as.data.frame()
phospho_long_df = phospho_data %>% gather(idr_type, phospho, all:C_w_pdb) %>% filter(!is.na(phospho)) %>% as.data.frame()


combine_n_c_term = function(data)
{
  n_c_tail_data = data %>% filter(idr_type %in% c('N', 'C')) %>% as.data.frame()
  n_c_tail_w_pdb_data = data %>% filter(idr_type %in% c('N_w_pdb', 'C_w_pdb')) %>% as.data.frame()
  
  n_c_tail_data$idr_type = 'tail'
  n_c_tail_w_pdb_data$idr_type = 'tail_w_pdb'
  
  data = bind_rows(data, n_c_tail_data)
  data = bind_rows(data, n_c_tail_w_pdb_data)
  
  data = data %>% filter((idr_type %in% c('N', 'C', 'N_w_pdb', 'C_w_pdb')) == FALSE) %>% as.data.frame()
  
  return(data)
}

fcr_long_df_processed = combine_n_c_term(fcr_long_df)
ncpr_long_df_processed = combine_n_c_term(ncpr_long_df)
fp_long_df_processed = combine_n_c_term(fp_long_df)
fn_long_df_processed = combine_n_c_term(fn_long_df)
idr_len_long_df_processed = combine_n_c_term(idr_len_long_df)
idr_len_sequence_len_long_df_processed = combine_n_c_term(idr_len_sequence_len_long_df)
kappa_long_df_processed = combine_n_c_term(kappa_long_df)
phospho_long_df_processed = combine_n_c_term(phospho_long_df)

all_features_long_df = cbind(fcr_long_df, 
                             ncpr_long_df,
                             fp_long_df,
                             fn_long_df,
                             idr_len_long_df,
                             idr_len_sequence_len_long_df,
                             kappa_long_df,
                             phospho_long_df)



plot_single_feature = function(data, feature, feature_str, comparison_type_str)
{
  
  if (comparison_type_str == 'All') { 
    data = data %>% filter(idr_type %in% c('all', 'tail', 'tail_w_pdb')) %>% as.data.frame()
    data[data$idr_type == 'all', 'idr_type'] = 'All'
    data[data$idr_type == 'tail', 'idr_type'] = 'Tail'
    data[data$idr_type == 'tail_w_pdb', 'idr_type'] = 'Tail W/ PDB'
    data$idr_type = factor(data$idr_type, levels = c('All', 'Tail', 'Tail W/ PDB'))
    plot_save_dir = './Output_Plots/All_Comparison'
  } else {
    data = data %>% filter(idr_type %in% c('non_tail', 'tail', 'tail_w_pdb')) %>% as.data.frame()
    data[data$idr_type == 'non_tail', 'idr_type'] = 'Non-Tail'
    data[data$idr_type == 'tail', 'idr_type'] = 'Tail'
    data[data$idr_type == 'tail_w_pdb', 'idr_type'] = 'Tail W/ PDB'
    data$idr_type = factor(data$idr_type, levels = c('Non-Tail', 'Tail', 'Tail W/ PDB'))
    plot_save_dir = './Output_Plots/Non_Tail_Comparison'
  }
  
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  if (feature_str == 'Kappa')
  {
    print(ggplot(data = data, aes_string(x = 'idr_type', y = feature)) + 
            geom_violin(trim=FALSE, fill="grey80", alpha = .5) + 
            labs(y = feature_str) + 
            ylim(0,1) +
            theme_classic() + 
            geom_boxplot(width=0.1, outlier.shape = NA) +
            theme(plot.title = element_text(size=15), text=element_text(size=25), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))
    plot_name = gsub('/','_',feature_str)
    plot_name = gsub('% ', '',plot_name)
    ggsave(paste(plot_save_dir, paste0(plot_name, '.tiff'), sep = '/'), width = 4, height = 8)
  }
  
  else
  {
    print(ggplot(data = data, aes_string(x = 'idr_type', y = feature)) + 
            geom_violin(trim=FALSE, fill="grey80", alpha = .5) + 
            labs(y = feature_str) + 
            theme_classic() + 
            geom_boxplot(width=0.1, outlier.shape = NA) +
            theme(plot.title = element_text(size=15), text=element_text(size=25), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))
    plot_name = gsub('/','_',feature_str)
    plot_name = gsub('% ', '',plot_name)
    ggsave(paste(plot_save_dir, paste0(plot_name, '.tiff'), sep = '/'), width = 4, height = 8)
    
    if (feature == 'idr_len')
    {
      data = data %>% filter(idr_len <= 300) %>% as.data.frame()
      print(ggplot(data = data, aes_string(x = 'idr_type', y = feature)) + 
              geom_violin(trim=FALSE, fill="grey80", alpha = .5) + 
              labs(y = feature_str) + 
              theme_classic() + 
              geom_boxplot(width=0.1, outlier.shape = NA) + 
              theme(plot.title = element_text(size=15), text=element_text(size=25), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))
      plot_name = paste(plot_name, 'lte_300', sep = '_')
      ggsave(paste(plot_save_dir, paste0(plot_name, '.tiff'), sep = '/'), width = 4, height = 8)
    }
  }
    
}

for (comparison_type in c('All', 'Non_Tail'))
{
  plot_single_feature(fcr_long_df_processed, 'fcr', 'FCR', comparison_type)
  plot_single_feature(ncpr_long_df_processed, 'ncpr', 'NCPR', comparison_type)
  plot_single_feature(fp_long_df_processed, 'fp', 'Fraction Positive', comparison_type)
  plot_single_feature(fn_long_df_processed, 'fn', 'Fraction Negative', comparison_type)
  plot_single_feature(idr_len_long_df_processed, 'idr_len', 'IDR Length', comparison_type)
  plot_single_feature(phospho_long_df_processed, 'phospho', '% Phospho', comparison_type)
  plot_single_feature(kappa_long_df_processed, 'kappa', 'Kappa', comparison_type)
}

plot_idr_len_power_law = function()
{
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  
  plot_save_dir = './Output_Plots/Non_Tail_Comparison'
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  get_hist_info = function(data)
  {
    hist_output = hist(data$idr_len, breaks = seq(10,300,10), include.lowest = FALSE)
    print(hist_output$breaks)
    print(hist_output$counts)
    return(data.frame(idr_len = hist_output$breaks, num = c(hist_output$counts,NA), idr_type = unique(data$idr_type)))
  }
  
  data = idr_len_long_df_processed
  #data = data %>% filter(idr_len <= 300) %>% filter(idr_type != 'all') %>% as.data.frame()
  data = data %>% filter(idr_type != 'all') %>% as.data.frame()
  

  data[data$idr_type == 'non_tail', 'idr_type'] = 'Non-Tail'
  data[data$idr_type == 'tail', 'idr_type'] = 'Tail'
  data[data$idr_type == 'tail_w_pdb', 'idr_type'] = 'Tail W/ PDB'
  data$idr_type = factor(data$idr_type, levels = c('Non-Tail', 'Tail', 'Tail W/ PDB'))
  
  data_median = data %>% group_by(idr_type) %>% summarise(grp.median = median(idr_len)) %>% as.data.frame()
  
  print(ggplot(data = data, aes(x = idr_len, color = idr_type, fill = idr_type)) + 
          stat_density(geom="line",position="identity", size = 1.5) +
          labs(x = 'IDR Length', y = 'Density') + 
          scale_color_manual(values = c(ggplotColours(3)[3], ggplotColours(3)[1], ggplotColours(3)[2])) +
          geom_vline(data=data_median, aes(xintercept=grp.median, color=idr_type),
                     linetype="dashed", size = 1, show.legend = F) +
          theme_classic() + 
          scale_x_continuous(breaks=seq(0, 300, 50)) +
          coord_cartesian(xlim=c(0,300)) +
          theme(legend.title = element_blank(), legend.position = c(.5,.5), text=element_text(size=20)))
  plot_name = paste('idr_len_all_density', sep = '_')
  ggsave(paste(plot_save_dir, paste0(plot_name, '.tiff'), sep = '/'), width=6, height=4)
  
  
  
  plot_save_dir = './Output_Plots/All_Comparison'
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  data = idr_len_long_df_processed
  data = data %>% filter(idr_type != 'non_tail') %>% as.data.frame()
  
  data[data$idr_type == 'all', 'idr_type'] = 'All'
  data[data$idr_type == 'tail', 'idr_type'] = 'Tail'
  data[data$idr_type == 'tail_w_pdb', 'idr_type'] = 'Tail W/ PDB'
  data$idr_type = factor(data$idr_type, levels = c('All', 'Tail', 'Tail W/ PDB'))
  
  data_median = data %>% group_by(idr_type) %>% summarise(grp.median = median(idr_len)) %>% as.data.frame()
  
  print(ggplot(data = data, aes(x = idr_len, color = idr_type, fill = idr_type)) + 
          stat_density(geom="line",position="identity", size = 1.5) +
          labs(x = 'IDR Length', y = 'Density') + 
          scale_color_manual(values = c(ggplotColours(3)[3], ggplotColours(3)[1], ggplotColours(3)[2])) +
          geom_vline(data=data_median, aes(xintercept=grp.median, color=idr_type),
                     linetype="dashed", size = 1, show.legend = F) +
          theme_classic() + 
          scale_x_continuous(breaks=seq(0, 300, 50)) +
          coord_cartesian(xlim=c(0,300)) +
          theme(legend.title = element_blank(), legend.position = c(.5,.5), text=element_text(size=20)))
  plot_name = paste('idr_len_all_density', sep = '_')
  ggsave(paste(plot_save_dir, paste0(plot_name, '.tiff'), sep = '/'), width=6, height=4)
                    
  
}



fcr_long_df_processed %>% group_by(idr_type) %>% summarise(mean(fcr)) %>% as.data.frame()
ncpr_long_df_processed %>% group_by(idr_type) %>% summarise(mean(ncpr)) %>% as.data.frame()
idr_len_sequence_len_long_df_processed %>% group_by(idr_type) %>% summarise(mean(idr_len_sequence_len)) %>% as.data.frame()
idr_len_sequence_len_long_df_processed %>% group_by(idr_type) %>% summarise(median(idr_len_sequence_len)) %>% as.data.frame()
idr_len_long_df_processed %>% group_by(idr_type) %>% summarise(mean(idr_len)) %>% as.data.frame()

pairwise.wilcox.test(phospho_long_df_processed$phospho, phospho_long_df_processed$idr_type, p.adjust.method="none")
pairwise.wilcox.test(kappa_long_df_processed$kappa, kappa_long_df_processed$idr_type, p.adjust.method="none")
pairwise.wilcox.test(ncpr_long_df_processed$ncpr, fcr_long_df_processed$idr_type, p.adjust.method="none")
pairwise.wilcox.test(fcr_long_df_processed$fcr, fcr_long_df_processed$idr_type, p.adjust.method="none")
pairwise.wilcox.test(fp_long_df_processed$fp, fp_long_df_processed$idr_type, p.adjust.method="none")
pairwise.wilcox.test(fn_long_df_processed$fn, fn_long_df_processed$idr_type, p.adjust.method="none")
pairwise.wilcox.test(idr_len_long_df_processed$idr_len, idr_len_long_df_processed$idr_type, p.adjust.method="none")

phospho_long_df_processed %>% wilcox_effsize(phospho ~ idr_type) %>% as.data.frame()
kappa_long_df_processed %>% wilcox_effsize(kappa ~ idr_type) %>% as.data.frame()
ncpr_long_df_processed %>% wilcox_effsize(ncpr ~ idr_type) %>% as.data.frame()
fcr_long_df_processed %>% wilcox_effsize(fcr ~ idr_type) %>% as.data.frame()
fp_long_df_processed %>% wilcox_effsize(fp ~ idr_type) %>% as.data.frame()
fn_long_df_processed %>% wilcox_effsize(fn ~ idr_type) %>% as.data.frame()
idr_len_long_df_processed %>% wilcox_effsize(idr_len ~ idr_type) %>% as.data.frame()




plot_phase_diagram_individual = function(data, idr_type_str)
{
  plot_save_dir = './Output_Plots'
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  data = data %>% filter(idr_type == idr_type_str) %>% as.data.frame()
  
  data[data$idr_type == 'all', 'idr_type'] = 'All'
  data[data$idr_type == 'tail', 'idr_type'] = 'Tail'
  data[data$idr_type == 'tail_w_pdb', 'idr_type'] = 'Tail W PDB'
  
  print(ggplot(data = data, aes(fp, fn) ) + 
          geom_density_2d_filled(binwidth=2) +
          ylim(0,.4) +
          xlim(0,.4) + 
          labs(y = 'Fraction Negative', x = 'Fraction Positive') + 
          theme_classic() + 
          theme(plot.title = element_text(size=15), text=element_text(size=25)))
  
  plot_name = paste0(paste('phase_diagram', idr_type_str, sep = '_'), '.tiff')
  ggsave(paste(plot_save_dir, plot_name, sep = '/'),dpi = 300)
  
}


plot_phase_diagram_facet = function(data)
{
  data = data %>% filter(idr_type != 'non_tail') %>% as.data.frame()
  
  plot_save_dir = './Output_Plots'
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  data[data$idr_type == 'all', 'idr_type'] = 'All'
  data[data$idr_type == 'tail', 'idr_type'] = 'Tail'
  data[data$idr_type == 'tail_w_pdb', 'idr_type'] = 'Tail W PDB'
  
  data$idr_type = factor(data$idr_type, levels = c('All', 'Tail', 'Tail W PDB'))
  
  print(ggplot(data = data, aes(fp, fn) ) + 
          geom_density_2d_filled(binwidth=2) +
          facet_wrap(vars(idr_type)) + 
          ylim(0,.4) +
          xlim(0,.4) + 
          labs(y = 'Fraction Negative', x = 'Fraction Positive') + 
          theme_classic() + 
          theme(plot.title = element_text(size=15), text=element_text(size=25), axis.text.x = element_text(angle = 45, hjust = 1)))
  
  plot_name = paste0(paste('phase_diagram_facet', sep = '_'), '.tiff')
  ggsave(paste(plot_save_dir, plot_name, sep = '/'),width = 8, height = 6)
  
}

phase_diagram_df = cbind(fp_long_df_processed,fn_long_df_processed)
phase_diagram_df = phase_diagram_df[,c(1,2,4)]

plot_phase_diagram_individual(phase_diagram_df, 'all')
plot_phase_diagram_individual(phase_diagram_df, 'tail')
plot_phase_diagram_individual(phase_diagram_df, 'tail_w_pdb')
plot_phase_diagram_facet(phase_diagram_df)
  
phase_diagram_df %>% group_by(idr_type) %>% filter((fp+fn) > .35) %>% filter(abs((fp-fn)) <= .35) %>% summarise(count = n()) %>% as.data.frame()
phase_diagram_df %>% group_by(idr_type) %>% filter((fp+fn) > .35) %>% filter(abs((fp-fn)) > .35) %>% summarise(count = n()) %>% as.data.frame()

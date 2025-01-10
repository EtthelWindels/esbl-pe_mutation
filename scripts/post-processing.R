# -----------------------------------------------------------------------------
## Process BEAST log files
## 2024-04-04 Etthel Windels
## -----------------------------------------------------------------------------

## Load libraries

library(dplyr)
library(stringr)
library(ape)
library(MCMCvis)
library(ggplot2)
library(ggtree)
library(lubridate)
library(treeio)
library(cowplot)

## ---------------------------

## Read files

coli_meta <- read.csv("path_to_coli_metadata")
kleb_meta <- read.csv("path_to_kleb_metadata")

LOG_coli_M1 <- "path_to_M1_coli_combined_log"
LOG_kleb_M1 <- "path_to_M1_kleb_combined_log"
LOG_coli_M2 <- "path_to_M2_coli_combined_log"
LOG_kleb_M2 <- "path_to_M2_kleb_combined_log"
LOG_coli_M3 <- "path_to_M3_coli_combined_log"
LOG_kleb_M3 <- "path_to_M3_kleb_combined_log"
LOG_coli_S11 <- "path_to_S11_coli_combined_log"
LOG_kleb_S11 <- "path_to_S11_kleb_combined_log"
LOG_coli_S21 <- "path_to_S21_coli_combined_log"
LOG_kleb_S21 <- "path_to_S21_kleb_combined_log"
LOG_coli_S22 <- "path_to_S22_coli_combined_log"
LOG_kleb_S22 <- "path_to_S22_kleb_combined_log"
LOG_coli_S23 <- "path_to_S23_coli_combined_log"
LOG_kleb_S23 <- "path_to_S23_kleb_combined_log"
LOG_coli_S24 <- "path_to_S24_coli_combined_log"
LOG_kleb_S24 <- "path_to_S24_kleb_combined_log"
LOG_coli_S31 <- "path_to_S31_coli_combined_log"
LOG_kleb_S31 <- "path_to_S31_kleb_combined_log"

## ---------------------------

## Data loading functions

load_trajectories <- function(filename, burninFrac=0, subsample=NA) {
  df_in <- as.matrix(read.table(filename, header=T))
  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }
  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }
  return(df_in)
}

get_parameter_summary <- function(logfile, parameter){
  summary = MCMCsummary(logfile, HPD=TRUE)
  parameter_summary = summary[parameter,]
  parameter_mean = parameter_summary[,1]
  parameter_hpdL = parameter_summary[,3]
  parameter_hpdU = parameter_summary[,4]
  parameter_median = median(select(as.data.frame(logfile), parameter)[,1])
  return(list(mean=parameter_mean, hpdL=parameter_hpdL, hpdU=parameter_hpdU, median=parameter_median))
}

extract_sampling_period <- function(species, patient){
  alignment <- ape::read.nexus.data(paste0("path_to_alignments/",species,"/pat",patient,"_concat_newid_varpos.nex"))
  dates <- str_split_fixed(names(alignment),"/",3)[,3]
  decimal_dates <- decimal_date(ymd(dates))
  sampling_period <- max(decimal_dates) - min(decimal_dates) 
}

extract_mrsd <- function(species, patient){
  alignment <- ape::read.nexus.data(paste0("path_to_alignments/",species,"/pat",patient,"_concat_newid_varpos.nex"))
  dates <- str_split_fixed(names(alignment),"/",3)[,3]
  decimal_dates <- decimal_date(ymd(dates))
  mrsd <- max(decimal_dates)
}

generate_dataset <- function(logfile, species, patient_list){
  log = as.data.frame(logfile, row.names = F)
  new_data <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(new_data) <- c("patient", "treeheight", "clockrate.mean", "clockrate.hpdL", "clockrate.hpdU", "clockrate.median")
  for (i in patient_list){
    variable_tree = paste0("Tree",i,".height")
    height = mean(select(log, variable_tree)[,1])
    clockrate.mean = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$mean)
    clockrate.hpdL = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$hpdL)
    clockrate.hpdU = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$hpdU)
    clockrate.median = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$median)
    new_data[nrow(new_data) + 1,] = c(i, height, clockrate.mean, clockrate.hpdL, clockrate.hpdU, clockrate.median)
  }
  new_data$treeheight <- as.numeric(new_data$treeheight)
  new_data$clockrate.mean <- as.numeric(new_data$clockrate.mean)
  new_data$clockrate.median <- as.numeric(new_data$clockrate.median)
  new_data$sampling_period <- rep(NA, dim(new_data)[1])
  for (i in patient_list){
    if (i %in% c("27.1","27.2","30.1","30.2","48.1","48.2")){
      i=str_replace(i, "\\.", "_")
    }
    new_data$sampling_period[new_data$patient==i] = extract_sampling_period(species,i)
  }
  return(new_data)
}

generate_dataset_ST <- function(logfile, species, ST_list){
  log = as.data.frame(logfile, row.names = F)
  new_data <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(new_data) <- c("patient", "treeheight", "clockrate.mean", "clockrate.hpdL", "clockrate.hpdU", "clockrate.median")
  for (i in ST_list){
    variable_tree = paste0("Tree",i,".height")
    height = mean(select(log, variable_tree)[,1])
    clockrate.mean = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$mean)
    clockrate.hpdL = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$hpdL)
    clockrate.hpdU = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$hpdU)
    clockrate.median = suppressWarnings(get_parameter_summary(logfile, paste0("clockRate",i))$median)
    new_data[nrow(new_data) + 1,] = c(i, height, clockrate.mean, clockrate.hpdL, clockrate.hpdU, clockrate.median)
  }
  new_data$treeheight <- as.numeric(new_data$treeheight)
  new_data$clockrate.mean <- as.numeric(new_data$clockrate.mean)
  new_data$clockrate.median <- as.numeric(new_data$clockrate.median)
  new_data$sampling_period <- rep(NA, dim(new_data)[1])
  return(new_data)
}

generate_clockdata <- function(logfile, species, patient_list){
  log = as.data.frame(logfile, row.names = F)    
  new_data <- as.data.frame(cbind(rep("prior",10000),rlnorm(10000,-13.815,1.25)))
  for (i in patient_list){
    variable = paste0("clockRate",i)
    clock = select(log, variable)[,1]
    new_data <- rbind(new_data, cbind(rep(i,length(clock)), as.numeric(clock)))
  }
  colnames(new_data) <- c("rep", "clock")
  new_data$clock <- as.numeric(new_data$clock)
  new_data$rep <- factor(new_data$rep, levels=c("prior", patient_list))
  return(new_data)
}

generate_clockdata_sensitivity <- function(logfile, species, patient_list){
  log = as.data.frame(logfile, row.names = F)    
  new_data <- as.data.frame(cbind(rep("prior",10000),rlnorm(10000,-14.59676,1.25)))
  for (i in patient_list){
    variable = paste0("clockRate",i)
    clock = select(log, variable)[,1]
    new_data <- rbind(new_data, cbind(rep(i,length(clock)), as.numeric(clock)))
  }
  colnames(new_data) <- c("rep", "clock")
  new_data$clock <- as.numeric(new_data$clock)
  new_data$rep <- factor(new_data$rep, levels=c("prior", patient_list))
  return(new_data)
}

generate_multiplierdata <- function(logfile, species, patient_list){
  new_data <- as.data.frame(cbind(rep("prior",10000),rlnorm(10000,0,0.5)))
  log = as.data.frame(logfile, row.names = F)
  for (i in patient_list){
    variable = paste0("multiplier",i)
    multiplier = select(log, variable)[,1]
    new_data <- rbind(new_data, cbind(rep(i,length(multiplier)), as.numeric(multiplier)))
  }
  colnames(new_data) <- c("rep", "multiplier")
  new_data$multiplier <- as.numeric(new_data$multiplier)
  new_data$rep <- factor(new_data$rep, levels=c("prior", patient_list))
  return(new_data)
}

## ---------------------------

## Plotting functions

plot_clockrate <- function(data, species, patient_list){
  if (species=="E. coli"){
    colour = "steelblue4"
  } else {
    colour = "darkorange2"
  }
  plot <- ggplot(data, aes(x=rep, y=clock, fill=rep, col=rep)) +
    geom_violin(aes(alpha=0.5)) +
    labs(title=species, y="Clock rate")+ 
    theme_classic() +
    theme(axis.text.x = element_text(hjust=0.5, size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.position = "none",
          plot.title = element_text(face="italic", size=14, hjust=0.5)) +
    ylim(0,1E-5) +
    scale_fill_manual(values=c('darkgrey', rep(colour,length(patient_list))), labels=c("prior",patient_list)) +
    scale_colour_manual(values=c('darkgrey', rep(colour,length(patient_list))), labels=c("prior",patient_list))
  return(plot)
}

plot_trees <- function(species, analysis, patient_list){
  plot_list <- lapply(patient_list, function(patient_nr){
    tree <- treeio::read.beast(paste0("path_to_results/",analysis,"/",species,"/",analysis,"_",species,".p",patient_nr,".mcc.txt"))  # MCC tree created with TreeAnnotator
    tree_data <- fortify(tree)
    MRSD <- extract_mrsd(species,str_replace(patient_nr, "\\.", "_"))
    mrsd <- str_split_fixed(date_decimal(MRSD),pattern=" ",3)[1]
    treeplot <- ggtree(tree,layout='rectangular', size=0.8, mrsd=mrsd, ignore.negative.edge=TRUE) # +
    treeplot %<+% tree_data + 
      geom_range(range='height_0.95_HPD', color='grey', alpha=.6, size=2) +
      theme_tree2(legend.position = "none",
                  plot.title = element_text(size=14, hjust=0.5)) +
      guides(fill = guide_legend(byrow = TRUE))+
      labs(title=paste0("Patient ", patient_nr), x="Year")
  })
  if(length(patient_list)==length(coli_nr_all)){
    plot <- plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
                      plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],
                      plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],
                      plot_list[[25]],plot_list[[26]],plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],
                      plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
                      plot_list[[41]],plot_list[[42]],plot_list[[43]],plot_list[[44]],plot_list[[45]],plot_list[[46]],plot_list[[47]],plot_list[[48]], ncol=6)
  } else if(length(patient_list)==length(kleb_nr_all)){
    plot <- plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
                      plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],
                      plot_list[[17]],plot_list[[18]],plot_list[[19]], ncol=6)
  } else if (length(patient_list)==length(coli_nr_sel)){
    plot <- plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
                      plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],
                      plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]], ncol=5)
  } else if (length(patient_list)==length(kleb_nr_sel)){
    plot <- plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]], ncol=5)    
  }
    return(plot)
}

plot_multiplier <- function(data, species, patient_list){
  if (species=="E. coli"){
    colour = "steelblue4"
  } else {
    colour = "darkorange2"
  }
  plot <- ggplot(data, aes(x=rep, y=multiplier, fill=rep, col=rep)) +
    geom_violin(aes(alpha=0.5)) +
    geom_hline(yintercept=1, linetype=2, size=0.75)+
    labs(title=species, y="Clock rate multiplier")+ 
    theme_classic() +
    theme(axis.text.x = element_text(hjust=0.5, size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.position = "none",
          plot.title = element_text(face="italic", size=14, hjust=0.5)) +
    ylim(0,5) +
    scale_fill_manual(values=c('darkgrey', rep(colour,length(patient_list))), labels=c("prior",patient_list)) +
    scale_colour_manual(values=c('darkgrey', rep(colour,length(patient_list))), labels=c("prior",patient_list))
  return(plot)
}

## ---------------------------

## Patient numbers

# patients with at least 3 samples
coli_nr_sel <- c("17","22","24","27.1","27.2","28","30.1","31","34","35","38","40","42","44","45","46","48.1","55","56","63","69","73")
kleb_nr_sel <- c("01","03","10","11","14","15","18")
coli_ST_sel <- c("2","3","4","6","11")
kleb_ST_sel <- c("8")

# all patients 
coli_nr_all <- c("17","18","21","22","24","26","27.1","27.2","28","29","30.1","30.2","31","32","33","34","35","37","38","39","40","41","42","43","44","45","46","47","48.1","48.2","50","51","52","53","55","56","60","61","62","63","64","65","66","68","69","71","72","73")
kleb_nr_all <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19")
coli_ST_all <- 1:19
kleb_ST_all <- 1:16

## ---------------------------

## Load log files

log_coli_M1 <- load_trajectories(LOG_coli_M1)
log_kleb_M1 <- load_trajectories(LOG_kleb_M1)
log_coli_M2 <- load_trajectories(LOG_coli_M2)
log_kleb_M2 <- load_trajectories(LOG_kleb_M2)
log_coli_M3 <- load_trajectories(LOG_coli_M3)
log_kleb_M3 <- load_trajectories(LOG_kleb_M3)
log_coli_S11 <- load_trajectories(LOG_coli_S11)
log_kleb_S11 <- load_trajectories(LOG_kleb_S11)
log_coli_S21 <- load_trajectories(LOG_coli_S21)
log_kleb_S21 <- load_trajectories(LOG_kleb_S21)
log_coli_S22 <- load_trajectories(LOG_coli_S22)
log_kleb_S22 <- load_trajectories(LOG_kleb_S22)
log_coli_S23 <- load_trajectories(LOG_coli_S23)
log_kleb_S23 <- load_trajectories(LOG_kleb_S23)
log_coli_S24 <- load_trajectories(LOG_coli_S24)
log_kleb_S24 <- load_trajectories(LOG_kleb_S24)
log_coli_S31 <- load_trajectories(LOG_coli_S31)
log_kleb_S31 <- load_trajectories(LOG_kleb_S31)


## ---------------------------

## M1 - SHARED CLOCK RATES

# Posterior means datasets

data_coli_M1 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(data_coli_M1) <- c("patient", "treeheight")
for (i in coli_nr_sel){
  variable = paste0("Tree",i,".height")
  log = as.data.frame(log_coli_M1, row.names = F)
  height = mean(select(log, variable)[,1])
  data_coli_M1[nrow(data_coli_M1) + 1,] = c(i, height)
}
data_coli_M1$treeheight <- as.numeric(data_coli_M1$treeheight)
data_coli_M1$sampling_period <- rep(NA, dim(data_coli_M1)[1])
for (i in coli_nr_sel){
  if (i %in% c("27.1","27.2","30.1","30.2","48.1","48.2")){
    i=str_replace(i, "\\.", "_")
  }
  data_coli_M1$sampling_period[data_coli_M1$patient==i] = extract_sampling_period("coli",i)
}

data_kleb_M1 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(data_kleb_M1) <- c("patient", "treeheight")
for (i in kleb_nr_sel){
  variable = paste0("Tree",i,".height")
  log = as.data.frame(log_kleb_M1, row.names = F)
  height = mean(select(log, variable)[,1])
  data_kleb_M1[nrow(data_kleb_M1) + 1,] = c(i, height)
}
data_kleb_M1$treeheight <- as.numeric(data_kleb_M1$treeheight)
data_kleb_M1$sampling_period <- rep(NA, dim(data_kleb_M1)[1])
for (i in kleb_nr_sel){
  data_kleb_M1$sampling_period[data_kleb_M1$patient==i] = extract_sampling_period("kleb",i)
}

# Clock rate datasets

set.seed(1)
priordata <- cbind(rep("prior",10000),rlnorm(10000,-13.815,1.25))

clockdata_coli_M1 <- as.data.frame(rbind(priordata,
                                         cbind(rep("posterior",length(as.data.frame(log_coli_M1)$clockRate)), as.data.frame(log_coli_M1)$clockRate)))
colnames(clockdata_coli_M1) <- c("rep", "clock")
clockdata_coli_M1$clock <- as.numeric(clockdata_coli_M1$clock)
clockdata_coli_M1$rep <- factor(clockdata_coli_M1$rep, levels = c("prior", "posterior"))

clockdata_kleb_M1 <- as.data.frame(rbind(priordata,
                                         cbind(rep("posterior",length(as.data.frame(log_kleb_M1)$clockRate)), as.data.frame(log_kleb_M1)$clockRate)))
colnames(clockdata_kleb_M1) <- c("rep", "clock")
clockdata_kleb_M1$clock <- as.numeric(clockdata_kleb_M1$clock)
clockdata_kleb_M1$rep <- factor(clockdata_kleb_M1$rep, levels = c("prior", "posterior"))

# Clock rate estimates

get_parameter_summary(log_coli_M1, "clockRate")
get_parameter_summary(log_kleb_M1, "clockRate")

clockplot_coli_M1 <- ggplot(clockdata_coli_M1, aes(x=rep, y=clock, fill=rep, col=rep)) +
  geom_violin(aes(alpha=0.5)) +
  labs(title=expression(atop(bold(italic("E. coli")))), y="Clock rate")+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5)) +
  ylim(0,1E-5) +
  scale_fill_manual(values=c('darkgrey', 'steelblue4'), labels=c("prior","posterior")) +
  scale_colour_manual(values=c('darkgrey', 'steelblue4'), labels=c("prior","posterior")) 

clockplot_kleb_M1 <- ggplot(clockdata_kleb_M1, aes(x=rep, y=clock, fill=rep, col=rep)) +
  geom_violin(aes(alpha=0.5)) +
  labs(title=expression(atop(bold(italic("K. pneumoniae")))), y="Clock rate") + 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5)) +
  ylim(0,1E-5) +
  scale_fill_manual(values=c('darkgrey', 'darkorange2'), labels=c("prior","posterior")) +
  scale_colour_manual(values=c('darkgrey', 'darkorange2'), labels=c("prior","posterior")) 

clockplot_M1 <- plot_grid(clockplot_coli_M1,clockplot_kleb_M1,ncol=2)

ggsave("M1_clock.png", plot=clockplot_M1, path="path_to_figures", width=250, height=150, units="mm", bg='white')
ggsave("M1_clock.svg", plot=clockplot_M1, path="path_to_figures", width=250, height=150, units="mm", bg='white')

## ---------------------------

## M2 - PATIENT-SPECIFIC CLOCK RATES (MULTIPLIER)

# Shared clock rate estimates

get_parameter_summary(log_coli_M2, "clock_shared")
get_parameter_summary(log_kleb_M2, "clock_shared")

# Posterior means datasets

coli_data_M2 <- generate_dataset(log_coli_M2, "coli", coli_nr_sel)
kleb_data_M2 <- generate_dataset(log_kleb_M2, "kleb", kleb_nr_sel)

# Clock rate datasets

clockdata_coli_M2 <- generate_clockdata(log_coli_M2, "coli", coli_nr_sel)
clockdata_kleb_M2 <- generate_clockdata(log_kleb_M2, "kleb", kleb_nr_sel)

# Multiplier datasets

multiplierdata_coli_M2 <- generate_multiplierdata(log_coli_M2, "coli", coli_nr_sel)
multiplierdata_kleb_M2 <- generate_multiplierdata(log_kleb_M2, "kleb", kleb_nr_sel)

# Clock rate estimates

clockplot_coli_M2 <- plot_clockrate(clockdata_coli_M2, "E. coli", coli_nr_sel)
clockplot_kleb_M2 <- plot_clockrate(clockdata_kleb_M2, "K. pneumoniae", kleb_nr_sel)

ggsave("M2_coli_clock.png", plot=clockplot_coli_M2, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("M2_coli_clock.svg", plot=clockplot_coli_M2, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("M2_kleb_clock.png", plot=clockplot_kleb_M2, path="/path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("M2_kleb_clock.svg", plot=clockplot_kleb_M2, path="path_to_figures", width=150, height=150, units="mm", bg='white')

# Multiplier estimates

multiplierplot_coli_M2 <- plot_multiplier(multiplierdata_coli_M2, "E. coli", coli_nr_sel)
multiplierplot_kleb_M2 <- plot_multiplier(multiplierdata_kleb_M2, "K. pneumoniae", kleb_nr_sel)

ggsave("M2_coli_multiplier.png", plot=multiplierplot_coli_M2, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("M2_coli_multiplier.svg", plot=multiplierplot_coli_M2, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("M2_kleb_multiplier.png", plot=multiplierplot_kleb_M2, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("M2_kleb_multiplier.svg", plot=multiplierplot_kleb_M2, path="path_to_figures", width=150, height=150, units="mm", bg='white')

# Tree plots

treeplot_coli_M2 <- plot_trees("coli", "M2-patspecifclock_multiplier_3sampl", coli_nr_sel)
treeplot_kleb_M2 <- plot_trees("kleb", "M2-patspecifclock_multiplier_3sampl", kleb_nr_sel)

ggsave("M2_coli_trees.png", plot=treeplot_coli_M2, path="path_to_figures", width=400, height=350, units="mm", bg='white')
ggsave("M2_coli_trees.svg", plot=treeplot_coli_M2, path="path_to_figures", width=400, height=350, units="mm", bg='white')
ggsave("M2_kleb_trees.png", plot=treeplot_kleb_M2, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("M2_kleb_trees.svg", plot=treeplot_kleb_M2, path="path_to_figures", width=400, height=150, units="mm", bg='white')

## ---------------------------

## M3 - SEQUENCE-TYPE-SPECIFIC CLOCK RATES (MULTIPLIER)

# Shared clock rate estimates

get_parameter_summary(log_coli_M3, "clock_shared")
get_parameter_summary(log_kleb_M3, "clock_shared")

# Posterior means datasets

coli_data_M3 <- generate_dataset_ST(log_coli_M3, "coli", coli_ST_sel)
kleb_data_M3 <- generate_dataset_ST(log_kleb_M3, "kleb", kleb_ST_sel)

# Clock rate datasets

clockdata_coli_M3 <- generate_clockdata(log_coli_M3, "coli", coli_ST_sel)
clockdata_kleb_M3 <- generate_clockdata(log_kleb_M3, "kleb", kleb_ST_sel)

# Multiplier datasets

multiplierdata_coli_M3 <- generate_multiplierdata(log_coli_M3, "coli", coli_ST_sel)
multiplierdata_kleb_M3 <- generate_multiplierdata(log_kleb_M3, "kleb", kleb_ST_sel)

# Clock rate estimates

clockplot_coli_M3 <- plot_clockrate(clockdata_coli_M3, "E. coli", coli_ST_sel)
clockplot_kleb_M3 <- plot_clockrate(clockdata_kleb_M3, "K. pneumoniae", kleb_ST_sel)

ggsave("M3_coli_clock.png", plot=clockplot_coli_M3, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("M3_coli_clock.svg", plot=clockplot_coli_M3, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("M3_kleb_clock.png", plot=clockplot_kleb_M3, path="path_to_figures", width=130, height=150, units="mm", bg='white')
ggsave("M3_kleb_clock.svg", plot=clockplot_kleb_M3, path="path_to_figures", width=130, height=150, units="mm", bg='white')

# Multiplier estimates

multiplierplot_coli_M3 <- plot_multiplier(multiplierdata_coli_M3, "E. coli", coli_ST_sel)
multiplierplot_kleb_M3 <- plot_multiplier(multiplierdata_kleb_M3, "K. pneumoniae", kleb_ST_sel)

ggsave("M3_coli_multiplier.png", plot=multiplierplot_coli_M3, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("M3_coli_multiplier.svg", plot=multiplierplot_coli_M3, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("M3_kleb_multiplier.png", plot=multiplierplot_kleb_M3, path="path_to_figures", width=130, height=150, units="mm", bg='white')
ggsave("M3_kleb_multiplier.svg", plot=multiplierplot_kleb_M3, path="path_to_figures", width=130, height=150, units="mm", bg='white')

## ---------------------------

## S1-1 - SHARED CLOCK RATES - ALL SAMPLES

# Posterior means datasets

data_coli_S11 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(data_coli_S11) <- c("patient", "treeheight")
for (i in coli_nr_all){
  variable = paste0("Tree",i,".height")
  log = as.data.frame(log_coli_S11, row.names = F)
  height = mean(select(log, variable)[,1])
  data_coli_S11[nrow(data_coli_S11) + 1,] = c(i, height)
}
data_coli_S11$treeheight <- as.numeric(data_coli_S11$treeheight)
data_coli_S11$sampling_period <- rep(NA, dim(data_coli_S11)[1])
for (i in coli_nr_all){
  if (i %in% c("27.1","27.2","30.1","30.2","48.1","48.2")){
    i=str_replace(i, "\\.", "_")
  }
  data_coli_S11$sampling_period[data_coli_S11$patient==i] = extract_sampling_period("coli",i)
}

data_kleb_S11 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(data_kleb_S11) <- c("patient", "treeheight")
for (i in kleb_nr_all){
  variable = paste0("Tree",i,".height")
  log = as.data.frame(log_kleb_S11, row.names = F)
  height = mean(select(log, variable)[,1])
  data_kleb_S11[nrow(data_kleb_S11) + 1,] = c(i, height)
}
data_kleb_S11$treeheight <- as.numeric(data_kleb_S11$treeheight)
data_kleb_S11$sampling_period <- rep(NA, dim(data_kleb_S11)[1])
for (i in kleb_nr_all){
  data_kleb_S11$sampling_period[data_kleb_S11$patient==i] = extract_sampling_period("kleb",i)
}

# Clock rate datasets

set.seed(1)
priordata <- cbind(rep("prior",10000),rlnorm(10000,-13.815,1.25))

clockdata_coli_S11 <- as.data.frame(rbind(priordata,
                                        cbind(rep("posterior",length(as.data.frame(log_coli_S11)$clockRate)), as.data.frame(log_coli_S11)$clockRate)))
colnames(clockdata_coli_S11) <- c("rep", "clock")
clockdata_coli_S11$clock <- as.numeric(clockdata_coli_S11$clock)
clockdata_coli_S11$rep <- factor(clockdata_coli_S11$rep, levels = c("prior", "posterior"))

clockdata_kleb_S11 <- as.data.frame(rbind(priordata,
                                        cbind(rep("posterior",length(as.data.frame(log_kleb_S11)$clockRate)), as.data.frame(log_kleb_S11)$clockRate)))
colnames(clockdata_kleb_S11) <- c("rep", "clock")
clockdata_kleb_S11$clock <- as.numeric(clockdata_kleb_S11$clock)
clockdata_kleb_S11$rep <- factor(clockdata_kleb_S11$rep, levels = c("prior", "posterior"))

# Clock rate estimates

get_parameter_summary(log_coli_S11, "clockRate")
get_parameter_summary(log_kleb_S11, "clockRate")

clockplot_coli_S11 <- ggplot(clockdata_coli_S11, aes(x=rep, y=clock, fill=rep, col=rep)) +
  geom_violin(aes(alpha=0.5)) +
  labs(title=expression(atop(bold(italic("E. coli")))), y="Clock rate")+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5)) +
  ylim(0,1E-5) +
  scale_fill_manual(values=c('darkgrey', 'steelblue4'), labels=c("prior","posterior")) +
  scale_colour_manual(values=c('darkgrey', 'steelblue4'), labels=c("prior","posterior")) 

clockplot_kleb_S11 <- ggplot(clockdata_kleb_S11, aes(x=rep, y=clock, fill=rep, col=rep)) +
  geom_violin(aes(alpha=0.5)) +
  labs(title=expression(atop(bold(italic("K. pneumoniae")))), y="Clock rate") + 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.title = element_text(size=18, hjust=0.5)) +
  ylim(0,1E-5) +
  scale_fill_manual(values=c('darkgrey', 'darkorange2'), labels=c("prior","posterior")) +
  scale_colour_manual(values=c('darkgrey', 'darkorange2'), labels=c("prior","posterior")) 

clockplot_S11 <- plot_grid(clockplot_coli_S11,clockplot_kleb_S11,ncol=2)

ggsave("S11_clock.png", plot=clockplot_S11, path="path_to_figures", width=250, height=150, units="mm", bg='white')
ggsave("S11_clock.svg", plot=clockplot_S11, path="path_to_figures", width=250, height=150, units="mm", bg='white')

## ---------------------------

## S2-1 - PATIENT-SPECIFIC CLOCK RATES (MULTIPLIER) - ALL SAMPLES

# Shared clock rate estimates

get_parameter_summary(log_coli_S21, "clock_shared")
get_parameter_summary(log_kleb_S21, "clock_shared")

# Posterior means datasets

coli_data_S21 <- generate_dataset(log_coli_S21, "coli", coli_nr_all)
kleb_data_S21 <- generate_dataset(log_kleb_S21, "kleb", kleb_nr_all)

# Clock rate datasets

clockdata_coli_S21 <- generate_clockdata(log_coli_S21, "coli", coli_nr_all)
clockdata_kleb_S21 <- generate_clockdata(log_kleb_S21, "kleb", kleb_nr_all)

# Multiplier datasets

multiplierdata_coli_S21 <- generate_multiplierdata(log_coli_S21, "coli", coli_nr_all)
multiplierdata_kleb_S21 <- generate_multiplierdata(log_kleb_S21, "kleb", kleb_nr_all)

# Clock rate estimates

clockplot_coli_S21 <- plot_clockrate(clockdata_coli_S21, "E. coli", coli_nr_all)
clockplot_kleb_S21 <- plot_clockrate(clockdata_kleb_S21, "K. pneumoniae", kleb_nr_all)

ggsave("S21_coli_clock.png", plot=clockplot_coli_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S21_coli_clock.svg", plot=clockplot_coli_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S21_kleb_clock.png", plot=clockplot_kleb_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S21_kleb_clock.svg", plot=clockplot_kleb_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')

# Multiplier estimates

multiplierplot_coli_S21 <- plot_multiplier(multiplierdata_coli_S21, "E. coli", coli_nr_all)
multiplierplot_kleb_S21 <- plot_multiplier(multiplierdata_kleb_S21, "K. pneumoniae", kleb_nr_all)

ggsave("S21_coli_multiplier.png", plot=multiplierplot_coli_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S21_coli_multiplier.svg", plot=multiplierplot_coli_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S21_kleb_multiplier.png", plot=multiplierplot_kleb_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S21_kleb_multiplier.svg", plot=multiplierplot_kleb_S21, path="path_to_figures", width=400, height=150, units="mm", bg='white')

## ---------------------------

## S2-2 - PATIENT-SPECIFIC CLOCK RATES - CLOCK PRIOR

# Shared clock rate estimates

get_parameter_summary(log_coli_S22, "clock_shared")
get_parameter_summary(log_kleb_S22, "clock_shared")

# Posterior means datasets

coli_data_S22 <- generate_dataset(log_coli_S22, "coli", coli_nr_sel)
kleb_data_S22 <- generate_dataset(log_kleb_S22, "kleb", kleb_nr_sel)

# Clock rate datasets

clockdata_coli_S22 <- generate_clockdata_sensitivity(log_coli_S22, "coli", coli_nr_sel)
clockdata_kleb_S22 <- generate_clockdata_sensitivity(log_kleb_S22, "kleb", kleb_nr_sel)

# Multiplier datasets

multiplierdata_coli_S22 <- generate_multiplierdata(log_coli_S22, "coli", coli_nr_sel)
multiplierdata_kleb_S22 <- generate_multiplierdata(log_kleb_S22, "kleb", kleb_nr_sel)

# Clock rate estimates

clockplot_coli_S22 <- plot_clockrate(clockdata_coli_S22, "E. coli", coli_nr_sel)
clockplot_kleb_S22 <- plot_clockrate(clockdata_kleb_S22, "K. pneumoniae", kleb_nr_sel)

ggsave("S22_coli_clock.png", plot=clockplot_coli_S22, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S22_coli_clock.svg", plot=clockplot_coli_S22, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S22_kleb_clock.png", plot=clockplot_kleb_S22, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("S22_kleb_clock.svg", plot=clockplot_kleb_S22, path="path_to_figures", width=150, height=150, units="mm", bg='white')

# Multiplier estimates

multiplierplot_coli_S22 <- plot_multiplier(multiplierdata_coli_S22, "E. coli", coli_nr_sel)
multiplierplot_kleb_S22 <- plot_multiplier(multiplierdata_kleb_S22, "K. pneumoniae", kleb_nr_sel)

ggsave("S22_coli_multiplier.png", plot=multiplierplot_coli_S22, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S22_coli_multiplier.svg", plot=multiplierplot_coli_S22, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S22_kleb_multiplier.png", plot=multiplierplot_kleb_S22, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("S22_kleb_multiplier.svg", plot=multiplierplot_kleb_S22, path="path_to_figures", width=150, height=150, units="mm", bg='white')

## ---------------------------

## S2-3 - PATIENT-SPECIFIC CLOCK RATES - EXPONENTIAL COALESCENT

# Shared clock rate estimates

get_parameter_summary(log_coli_S23, "clock_shared")
get_parameter_summary(log_kleb_S23, "clock_shared")

# Posterior means datasets

coli_data_S23 <- generate_dataset(log_coli_S23, "coli", coli_nr_sel)
kleb_data_S23 <- generate_dataset(log_kleb_S23, "kleb", kleb_nr_sel)

# Clock rate datasets

clockdata_coli_S23 <- generate_clockdata(log_coli_S23, "coli", coli_nr_sel)
clockdata_kleb_S23 <- generate_clockdata(log_kleb_S23, "kleb", kleb_nr_sel)

# Multiplier datasets

multiplierdata_coli_S23 <- generate_multiplierdata(log_coli_S23, "coli", coli_nr_sel)
multiplierdata_kleb_S23 <- generate_multiplierdata(log_kleb_S23, "kleb", kleb_nr_sel)

# Clock rate estimates

clockplot_coli_S23 <- plot_clockrate(clockdata_coli_S23, "E. coli", coli_nr_sel)
clockplot_kleb_S23 <- plot_clockrate(clockdata_kleb_S23, "K. pneumoniae", kleb_nr_sel)

ggsave("S23_coli_clock.png", plot=clockplot_coli_S23, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S23_coli_clock.svg", plot=clockplot_coli_S23, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S23_kleb_clock.png", plot=clockplot_kleb_S23, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("S23_kleb_clock.svg", plot=clockplot_kleb_S23, path="path_to_figures", width=150, height=150, units="mm", bg='white')

# Multiplier estimates

multiplierplot_coli_S23 <- plot_multiplier(multiplierdata_coli_S23, "E. coli", coli_nr_sel)
multiplierplot_kleb_S23 <- plot_multiplier(multiplierdata_kleb_S23, "K. pneumoniae", kleb_nr_sel)

ggsave("S23_coli_multiplier.png", plot=multiplierplot_coli_S23, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S23_coli_multiplier.svg", plot=multiplierplot_coli_S23, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S23_kleb_multiplier.png", plot=multiplierplot_kleb_S23, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("S23_kleb_multiplier.svg", plot=multiplierplot_kleb_S23, path="path_to_figures", width=150, height=150, units="mm", bg='white')

## ---------------------------

## S2-4 - PATIENT-SPECIFIC CLOCK RATES (INDEPENDENT)

# Posterior means datasets

coli_data_S24 <- generate_dataset(log_coli_S24, "coli", coli_nr_sel)
kleb_data_S24 <- generate_dataset(log_kleb_S24, "kleb", kleb_nr_sel)

# Clock rate datasets

clockdata_coli_S24 <- generate_clockdata(log_coli_S24, "coli", coli_nr_sel)
clockdata_kleb_S24 <- generate_clockdata(log_kleb_S24, "kleb", kleb_nr_sel)

ggsave("S24_clockhist.png", plot=clockhist_S24, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S24_clockhist.svg", plot=clockhist_S24, path="path_to_figures", width=300, height=150, units="mm", bg='white')

# Clock rate estimates

clockplot_coli_S24 <- plot_clockrate(clockdata_coli_S24, "E. coli", coli_nr_sel)
clockplot_kleb_S24 <- plot_clockrate(clockdata_kleb_S24, "K. pneumoniae", kleb_nr_sel)

ggsave("S24_coli_clock.png", plot=clockplot_coli_S24, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S24_coli_clock.svg", plot=clockplot_coli_S24, path="path_to_figures", width=300, height=150, units="mm", bg='white')
ggsave("S24_kleb_clock.png", plot=clockplot_kleb_S24, path="path_to_figures", width=150, height=150, units="mm", bg='white')
ggsave("S24_kleb_clock.svg", plot=clockplot_kleb_S24, path="path_to_figures", width=150, height=150, units="mm", bg='white')

## ---------------------------

## S3-1 - ST-SPECIFIC CLOCK RATES (MULTIPLIER) - ALL SAMPLES

# Shared clock rate estimates

get_parameter_summary(log_coli_S31, "clock_shared")
get_parameter_summary(log_kleb_S31, "clock_shared")

# Posterior means datasets

coli_data_S31 <- generate_dataset_ST(log_coli_S31, "coli", coli_ST_all)
kleb_data_S31 <- generate_dataset_ST(log_kleb_S31, "kleb", kleb_ST_all)

# Clock rate datasets

clockdata_coli_S31 <- generate_clockdata(log_coli_S31, "coli", coli_ST_all)
clockdata_kleb_S31 <- generate_clockdata(log_kleb_S31, "kleb", kleb_ST_all)

# Multiplier datasets

multiplierdata_coli_S31 <- generate_multiplierdata(log_coli_S31, "coli", coli_ST_all)
multiplierdata_kleb_S31 <- generate_multiplierdata(log_kleb_S31, "kleb", kleb_ST_all)

# Clock rate estimates

clockplot_coli_S31 <- plot_clockrate(clockdata_coli_S31, "E. coli", coli_ST_all)
clockplot_kleb_S31 <- plot_clockrate(clockdata_kleb_S31, "K. pneumoniae", kleb_ST_all)

ggsave("S31_coli_clock.png", plot=clockplot_coli_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S31_coli_clock.svg", plot=clockplot_coli_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S31_kleb_clock.png", plot=clockplot_kleb_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S31_kleb_clock.svg", plot=clockplot_kleb_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')

# Multiplier estimates

multiplierplot_coli_S31 <- plot_multiplier(multiplierdata_coli_S31, "E. coli", coli_ST_all)
multiplierplot_kleb_S31 <- plot_multiplier(multiplierdata_kleb_S31, "K. pneumoniae", kleb_ST_all)

ggsave("S31_coli_multiplier.png", plot=multiplierplot_coli_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S31_coli_multiplier.svg", plot=multiplierplot_coli_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S31_kleb_multiplier.png", plot=multiplierplot_kleb_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
ggsave("S31_kleb_multiplier.svg", plot=multiplierplot_kleb_S31, path="path_to_figures", width=400, height=150, units="mm", bg='white')
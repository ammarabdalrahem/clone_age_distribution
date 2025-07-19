---
title: "Clone age analysis"
author: "Ammar Abdalrahem"
output: html_document
date: "2024-10-21"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

------------------------------------------------------------------------

## 1. Dependencies

Remember to re-run this code every time you re-open this R Notebook.

```{r, echo=FALSE}
#Code to install packages if necessary, and read them with library function

required_packages <- c("knitr","rmarkdown","ggridges","data.table","ggplot2","gganimate","dplyr","tidyr","cowplot", "reshape2","ggcorrplot","corrplot","viridis")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

```

## 2. Retrieve data



```{r, echo=FALSE}
# Get the path of the current R script
path <- dirname(rstudioapi::getSourceEditorContext()$path)

# Set the working directory to the path of the current R script
setwd(path)

# Check the current working directory
#getwd()

#new data (different mutation rates)
file_path_dif_mutation <- "2024-12-12-17h 55minsimulationS3_without_burnin.txt"
data_dif_mutation <- fread(file_path_dif_mutation, header = TRUE, sep = "\t", fill = TRUE)

colnames(data_dif_mutation) <- c("Replicate","Generation","Mutation_Rate", "clonalrate", "Nb_alleles_tot", "Number_alleles_Tot" , "Number_fixed_loci_Tot", "Mean_He_Tot" ,  "Mean_Ho_Tot" , "Mean_FIS_Tot", "Var_FIS_Tot", "Number_Genotypes_Tot","R_Tot", "Pareto_beta_Tot","List_distribution_gen_clonal_genotypes")


#Remove rows contain NA 
#data_dif_mutation <- data_dif_mutation %>% drop_na() 

# Filter the data, alleles 4 and mutation rate 1e-06 and generation every 10
data <- data_dif_mutation %>% filter(Nb_alleles_tot == "4") %>% filter(Mutation_Rate == 0.001) %>% filter (Generation %% 10 == 0) 

#data_dif_mutation_m4_N4$Mean_FIS_Tot <- as.numeric(data_dif_mutation_m4_N4$Mean_FIS_Tot)
```



## The relation between population genetic indices and clonerate 
# As already shown in Stoeckel et al. (2021), this relationship, with an initial population size of 10^4 Figure S1a.
# Reproduce the analysis again 


```{r, echo=FALSE}

data_filt_6000 <- data %>% filter (Generation == 6000)


# Create the boxplot for R
plot_R <- ggplot(data_filt_6000, aes(x = as.factor(clonalrate), y = R_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "R"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Create the boxplot for Pareto_beta
plot_Pareto_beta <- ggplot(data_filt_6000, aes(x = as.factor(clonalrate), y = Pareto_beta_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression(italic(β) * " Pareto")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


#Creat the boxplot for Var_FIS
plot_Var_FIS <- ggplot(data_filt_6000, aes(x = as.factor(clonalrate), y = Var_FIS_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = "Clonal rate",
    y = expression("Var("*italic(F)[IS]*")")
) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Creat the boxplot for Mean_FIS
plot_Mean_FIS <- ggplot(data_filt_6000, aes(x = as.factor(clonalrate), y = Mean_FIS_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = "Clonal rate",
    y = expression("Mean("*italic(F)[IS]*")")
) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Creat the boxplot for Mean_He
plot_Mean_He <- ggplot(data_filt_6000, aes(x = as.factor(clonalrate), y = Mean_He_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression(italic(H)[e])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

#Creat the boxplot for Mean_Ho
plot_Mean_Ho <- ggplot(data_filt_6000, aes(x = as.factor(clonalrate), y = Mean_Ho_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression(italic(H)[o])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )




# Combine the plots 
# cowplot package
combined_plot <- plot_grid( plot_Mean_He, plot_Mean_Ho, plot_R , plot_Pareto_beta,  plot_Var_FIS, plot_Mean_FIS,   label_size = 12, labels = c("a", "b", "c", "d", "e","f"), ncol = 2)

# Save plot
ggsave(
  filename = "combined_pop_indices_plot.png",  
  plot = combined_plot,           
  width = 18,  # Double-column width (cm)
  height = 15,   # Adjust height to maintain aspect ratio
  units = "cm",                   
  dpi = 1200                       
)

combined_plot



```



# prepare data for joyplot

```{r, echo=FALSE}
# Create a list to store clone ages for each group
clone_ages_by_group <- list()

# Process (loop) the data by row 
for (i in seq_len(nrow(data))) {
  group_value <- as.character(data[i, 4])  # Convert clonerate to character (ex. "0.7")
  age_distribution <- gsub("[{}]", "", data[i, 15])  # Remove curly braces (ex. {0: 300, 1: 213} → "0: 300, 1: 213")
  split_distribution <- strsplit(age_distribution, ", ")[[1]]  # Split by comma (ex. "0: 300" "1: 213")
  
  # If the group not included, add it
  if (!group_value %in% names(clone_ages_by_group)) {
    clone_ages_by_group[[group_value]] <- c()
  }
  
  # Process each age:count pair
  for (pair in split_distribution) {
    # Split each pair by colon and directly append the replicated ages
    split_pair <- strsplit(pair, ":")[[1]] #(ex. "0: 300" )
    age <- as.numeric(split_pair[1])       #(ex. 0)
    count <- as.numeric(split_pair[2])     #(ex. 300)
    
    # Append the replicated ages to the group
    clone_ages_by_group[[group_value]] <- c(clone_ages_by_group[[group_value]], rep(age, count))
  }
}

# Combine the data for all groups into one df
all_ages <- c()
all_groups <- c()

for (group_key in names(clone_ages_by_group)) {
  clone_ages <- clone_ages_by_group[[group_key]]
  if (length(clone_ages) > 0) {
    all_ages <- c(all_ages, clone_ages)
    all_groups <- c(all_groups, rep(group_key, length(clone_ages)))
  }
}



```

##joyplot

```{r, echo=FALSE}
# Create a data frame for plotting
df <- data.frame(Clone_Age = all_ages, Group = factor(all_groups, levels = rev(sort(unique(all_groups)))))


# Joyplot 
Joyplot_plot<- ggplot(df, aes(x = Clone_Age, y = Group, fill = Group)) +
  geom_density_ridges(scale = 3, alpha = 0.7) +
  labs(
    x = "Generations",
    y = "Clonal rate"
  ) +
  theme_ridges() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.title.x = element_text(size = 14, hjust = 0.5),
    axis.title.y = element_text(size = 14, angle = 90,  hjust = 0.5),
    legend.position = "none",
    
  )

# Save plot
ggsave(
  filename = "Joyplot_clone_ages_by_clonerateplot.png",  
  plot = Joyplot_plot,           
  width = 15,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
# Check the lowest value in clone age 
#min(df$Clone_Age)
#  curve height indicating the number of individuals at each age 

```


## Joyplot limit to 40

```{r, echo=FALSE}
# Joyplot limit to 40
joyplot_40 <- ggplot(df, aes(x = Clone_Age, y = Group, fill = Group)) +
  geom_density_ridges(scale = 3, alpha = 0.7) +
  #scale_fill_viridis_d(option = "C") +
  labs(
    #title = "Clone Ages by Clonerate (Limited to 40)",
    x = "Generations",
    y = "Clonal rate"
  ) +
  xlim(NA ,40) +  # Set x-axis limits to end at 40
  theme_ridges() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    #plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14, hjust = 0.5),
    axis.title.y = element_text(size = 14, angle = 90,  hjust = 0.5),
    legend.position = "none",
    
  )



# Save plot
ggsave(
  filename = "Joyplot_clone_ages_40_by_clonerate_plot.png",  
  plot = joyplot_40,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)

```


## Log scale plot

```{r, echo=FALSE}

# Prepare x-values (age keys) and absolute frequencies for each group
age_keys <- sort(unique(unlist(clone_ages_by_group)))  # Get unique ages

# Create a df directly from the clone ages
df_log <- data.frame(
  Age = rep(age_keys, times = length(clone_ages_by_group)),
  Count = unlist(lapply(clone_ages_by_group, function(ages) {
    sapply(age_keys, function(age) sum(ages == age))
  })),
  Group = rep(names(clone_ages_by_group), each = length(age_keys))
)



# Create the plot
log_age<- ggplot(df_log, aes(x = Age, y = Count, color = Group)) +
  geom_line(alpha = 0.6) +  # transparency
  scale_y_log10() +  # Set y-axis to log scale
  scale_x_log10() +  # Set x-axis to log scale
  labs(
    title = "Absolute frequencies of clone ages by group (Log scale)",
    x = "Age of clones (Log scale)",
    y = "Number of individuals",
    color = "Clone group"
  ) +
  theme_minimal() +  
  theme(
    legend.position = "bottom",  
    legend.direction = "horizontal",  # Arrange legend horizontally
    legend.box = "horizontal" 
  )

# Save plot
ggsave(
  filename = "log_ages_by_clonerate_plot.png",  
  plot = log_age,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 300                       
)

```


## Scatter plot for all replicates

```{r, echo=FALSE}

parsed_data <- data %>%
  rowwise() %>%
  mutate(parsed_list = strsplit(gsub("[{}]", "", List_distribution_gen_clonal_genotypes), ", ")) %>%
  unnest(parsed_list) %>%
  separate(parsed_list, into = c("Age", "Count"), sep = ":") %>%
  mutate(
    Age = as.numeric(Age),
    Count = as.numeric(Count)
  ) %>%
  ungroup()


# Only specific clonerate 
# Create movie 0 to 1 clonerate
#parsed_data_0.98 <- parsed_data %>% filter(clonalrate == "0.1", Replicate== "140")


# Plot the data 
scatter_plot_all_ages <- ggplot(parsed_data, aes(x = Age, y = Count, color = as.factor(clonalrate))) +
  geom_point(alpha = 0.7, size = 2) +  # Increase point size slightly
  scale_y_log10() +
  scale_x_log10() +
  #scale_color_viridis_d(option = "C", end = 0.9) +  # Use perceptually-uniform colors
  labs(
    x = "Clone age (log scale)",
    y = "Number of individuals (log scale)",
    color = "Clonal rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(nrow = 1))

print(scatter_plot_all_ages)

# Save plot
ggsave(
  filename = "scatter_ages_by_clonerate_plot.png",  
  plot = scatter_plot_all_ages,           
  width = 40,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)


# Animated HTML
p <- ggplot(parsed_data, aes(x = Age, y = Count, color = as.factor(clonalrate))) +
  geom_point(alpha = 0.7) +  # transparency
  scale_y_log10() +  # Set y-axis to log scale
  scale_x_log10() +  # Set x-axis to log scale
  labs(
    title = "Number of individuals vs. Clone age by clonal rate",
    x = "Age of clones ",
    y = "Number of individuals (log scale)",
    color = "Clonal rate"
  ) +
  theme_minimal() +  
  theme(
    legend.position = "top",  
    legend.direction = "horizontal"  
  ) +
  guides(color = guide_legend(nrow = 1)) +  # Arrange legend in one row
  transition_states(clonalrate, transition_length = 2, state_length = 1) +  # Transition based on clonal rate
  labs(title = "Clone age distribution for clonal rate: {closest_state}")  # Dynamic title showing clonal rate

# Animate the plot
animate(p, nframes = 100, fps = 10, width = 1500, height = 800, res = 150, renderer = gifski_renderer())
#anim <- animate(p, nframes = 100, fps = 10, width = 1500, height = 800, res = 150, renderer = gifski_renderer())

# Save animation
#anim_save("clone_age_animation_wide_high_res.gif", animation = anim)



```



## The relation between population genetic indices, generations and clonerate

```{r, echo=FALSE}
# filter the data to minimize the clonal rates for clear visualization
data_min_c_10_6 <- data %>% filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999", "0.9999"))

# group by generation and calculate the median per each indices
data_grouped_g_c_10_6 <- data_min_c_10_6 %>%
  group_by(Generation, clonalrate) %>%
  summarise(
    Mean_FIS_Tot = median(Mean_FIS_Tot, na.rm = TRUE),
    Var_FIS_Tot = median(Var_FIS_Tot, na.rm = TRUE),
    Mean_He_Tot = median(Mean_He_Tot, na.rm = TRUE),
    Mean_Ho_Tot = median(Mean_Ho_Tot, na.rm = TRUE),
    R_Tot = median(R_Tot, na.rm = TRUE),
    Pareto_beta_Tot = median(Pareto_beta_Tot, na.rm = TRUE)
  ) %>%
  ungroup() # 


# Plot mean FIS by clonal rate and generation
plot_Mean_FIS_10_6 <-  ggplot(data_grouped_g_c_10_6, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Mean("*italic(F)[IS]*")")) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot variance of FIS by clonal rate and generation
plot_Var_FIS_10_6 <- ggplot(data_grouped_g_c_10_6, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Var_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Var("*italic(F)[IS]*")")) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


 # Plot mean He by clonal rate and generation
plot_Mean_He_10_6 <- ggplot(data_grouped_g_c_10_6, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_He_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[e])) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot mean Ho by clonal rate and generation
plot_Mean_Ho_10_6 <- ggplot(data_grouped_g_c_10_6, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_Ho_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[o])) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Plot R by clonal rate and generation
plot_R_10_6 <-  ggplot(data_grouped_g_c_10_6, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = R_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(R))) +
  labs(
    x = "Generation",
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(margin = margin(r = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot Pareto beta by clonal rate and generation
plot_Pareto_beta_10_6 <- ggplot(data_grouped_g_c_10_6, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Pareto_beta_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(β) * " Pareto")) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )

```


## data for mutation rate 10-4

```{r, echo=FALSE}
# Filter the data for mutation rate 10^-4
data_mutation_10_4 <- data_dif_mutation %>% filter(Nb_alleles_tot == "4") %>% filter(Mutation_Rate == 0.0001) %>% filter (Generation %% 10 == 0) 



```


## The relation between population genetic indices, generations and clonerate for mutation rate 10^-4

```{r, echo=FALSE}
# filter the data to minimize the clonal rates for clear visualization
data_min_c_10_4 <- data_mutation_10_4 %>% filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999", "0.9999"))

# group by generation and calculate the median per each indices
data_grouped_g_c_10_4 <- data_min_c_10_4 %>%
  group_by(Generation, clonalrate) %>%
  summarise(
    Mean_FIS_Tot = median(Mean_FIS_Tot, na.rm = TRUE),
    Var_FIS_Tot = median(Var_FIS_Tot, na.rm = TRUE),
    Mean_He_Tot = median(Mean_He_Tot, na.rm = TRUE),
    Mean_Ho_Tot = median(Mean_Ho_Tot, na.rm = TRUE),
    R_Tot = median(R_Tot, na.rm = TRUE),
    Pareto_beta_Tot = median(Pareto_beta_Tot, na.rm = TRUE)
  ) %>%
  ungroup() # 


# Plot mean FIS by clonal rate and generation
plot_Mean_FIS_10_4 <-  ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Mean("*italic(F)[IS]*")")) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 30), # 30 pts on the left
    axis.title.y = element_blank(),
  )


# Plot variance of FIS by clonal rate and generation
plot_Var_FIS_10_4 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Var_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Var("*italic(F)[IS]*")")) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 30), # 30 pts on the left
    axis.title.y = element_blank()
  )


 # Plot mean He by clonal rate and generation
plot_Mean_He_10_4 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_He_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[e])) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 30), # 30 pts on the left
    axis.title.y = element_blank()
  )


# Plot mean Ho by clonal rate and generation
plot_Mean_Ho_2 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_Ho_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[o])) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 30), # 30 pts on the left
    axis.title.y = element_blank()
  )

# Plot R by clonal rate and generation
plot_R_10_4 <-  ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = R_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(R))) +
  labs(
    x = "Generation",
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(margin = margin(r = 10)),
    axis.title.y = element_blank(),
    plot.margin = margin(5, 5, 5, 30)  # 30 pts on the left

  )


# Plot Pareto beta by clonal rate and generation
plot_Pareto_beta_10_4 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Pareto_beta_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(β) * " Pareto")) +
  labs(
    y = "Clonal rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(5, 5, 5, 30) # 30 pts on the left
  )

```


```{r, echo=FALSE}
# Combine the plots 
combined_plot_2 <- plot_grid (plot_Mean_FIS_10_6, plot_Mean_FIS_10_4, plot_Var_FIS_10_6, plot_Var_FIS_10_4, plot_Mean_He_10_6,plot_Mean_He_10_4, plot_R_10_6, plot_R_10_4
                              ,  label_size = 12, labels = c("a", "b", "c", "d","e","f","g","h"), ncol = 2)


combined_plot_2


# Save plot
ggsave(
  filename = "combined_pop_indices_plot_2.png",
  plot = combined_plot_2,
  width = 30,  
  height = 25,  
  units = "cm",
  dpi = 1200
)



```


## R^2 per each replicate
# we calculate R^2 for each replicate and clonal rate at generation 6000
```{r, echo=FALSE}
# data for generation 6000
parsed_data_g <- parsed_data %>%
  filter(Generation == 6000) 

# Initialize a data frame to store the results
r2_results <- data.frame(clonalrate = character(), Replicate = integer(), R2 = numeric(),slope=numeric(), stringsAsFactors = FALSE)

# Get unique clonal rates
unique_clonal_rates <- unique(parsed_data_g$clonalrate)

# Loop over each clonal rate to calculate R^2 and slopes
for (clonal_rate in unique_clonal_rates) {
  # Filter data for the current clonal rate
  parsed_data_clonal <- parsed_data_g %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log-transform Count

  # Loop over unique replicates for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate
    replicate_data <- parsed_data_clonal %>%
      filter(Replicate == replicate)

    # Calculate R^2  and slopes if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data)  # Fit linear model
      r2 <- summary(lm_fit)$r.squared  # Extract R^2 from the model summary
      # Calculate slope
      slope <- coef(lm_fit)["Age"]  # Extract slope (coefficient for Age)
    } else {
      r2 <- NA  # If not enough data, assign NA
      slope <- NA  # If not enough data, assign NA
    }

    # Store the results
    r2_results <- rbind(r2_results, data.frame(clonalrate = clonal_rate, Replicate = replicate, R2 = r2, slope = slope))
  }
}


# Remove rows where NA
r2_results <- r2_results %>% filter(!is.na(R2))
r2_results <- r2_results %>% filter(!is.na(slope))

# Convert 'clonalrate' to a factor
r2_results$clonalrate <- as.factor(r2_results$clonalrate)

# Display the results
head(r2_results)


# max,mod, var age
# Calculate the maximum, median, and variance of clone ages for each clonal rate and replicate

age_stats <- parsed_data_g %>%
  group_by(clonalrate, Replicate) %>%
  summarise(
    Max_Age = max(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    Var_Age = var(Age, na.rm = TRUE),
    .groups = "drop"
  )


# Plot 
R_2_plot<- ggplot(r2_results, aes(x = clonalrate, y = R2)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +  labs(
    y = expression(R^2),
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


R_2_plot

# Plot slope
slope_plot <- ggplot(r2_results, aes(x = clonalrate, y = slope)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) + labs(
    x = "Clonal rate",
    y = "Slope",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

slope_plot

# combine the plots
combined_r2_plot <- plot_grid(slope_plot, R_2_plot, labels = c("a", "b"), ncol = 1)
# Save plot
ggsave(
  filename = "combined_r2_plot.png",  
  plot = combined_r2_plot,           
  width = 18,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)

combined_r2_plot
```

## clonal age summary as function of genetic indices 6000


```{r}
# Merge slope results with population genetics data on clonal rate and replicate
data_6000 <- data %>% filter(Generation == 6000) 
# Convert factors to character before merging
r2_results$clonalrate <- as.character(r2_results$clonalrate)
data_6000$clonalrate <- as.character(data_6000$clonalrate)
age_stats$clonalrate <- as.character(age_stats$clonalrate)

#merge by replicate 

merged_slope_parsed_df <- merge(r2_results, data_6000 , age_stats %>% select(clonalrate, Replicate, Mean_FIS_Tot, Var_FIS_Tot, Mean_Ho_Tot, Mean_He_Tot, R_Tot, Pareto_beta_Tot), by.x = c("clonalrate", "Replicate"), by.y = c("clonalrate", "Replicate"))

# Minimize the  clonal rates for better visualization
merged_slope_parsed_data <- merged_slope_parsed_df %>%
  filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999","0.9999"))

```




## spearman correlation using merged_slope_parsed_data

```{r, echo=FALSE}
# Prepare numeric-only data (exclude categorical variables)
cor_data <- merged_slope_parsed_data %>%
  dplyr::select(-clonalrate, -Replicate)

# Compute Spearman correlation matrix
cor_matrix <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")

# Rename columns and rows for better readability
colnames(cor_matrix) <- c("R²", "Slope", "Mean (FIS)", "Var (FIS)", "Mean (Ho)", "Mean (He)", "R", "Pareto β")
rownames(cor_matrix) <- c("R²", "Slope", "Mean (FIS)", "Var (FIS)", "Mean (Ho)", "Mean (He)", "R", "Pareto β")

# Set up file path for saving
file_path <- "Spearman_correlation_plot_merged_data.png"
png(height = 20, width = 20, file = file_path, units = "cm", res = 1200, bg = "white")

# Create correlation plot
Spearman_correlation_plot <- corrplot(cor_matrix,
                                     method = "color",
                                     type = "upper",
                                     tl.cex = 0.9,                    # Slightly smaller text
                                     addCoef.col = "black",
                                     tl.col = "black",
                                     mar = c(0,0,1,0),
                                     bg = "white",
                                     number.cex = 0.7,                # Reduced coefficient text size
                                     col = viridis(200, option = "viridis"),  # Colorblind-friendly palette
                                     tl.srt = 45,                     # Rotate labels for better readability
                                     cl.cex = 0.8,                    # Color legend text size
                                     cl.ratio = 0.2)                 # Color legend width



# Close the PNG device
dev.off()



```




# Plotting the relationship between slope and various population genetic indices

```{r, echo=FALSE}

# Plot Slope vs Var_FIS for each clonal rate
slope_var_fis <- ggplot(merged_slope_parsed_data, aes(x = slope , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope" ,
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Mean_FIS for each clonal rate
slope_mean_fis <- ggplot(merged_slope_parsed_data, aes(x = slope, y =  Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope", 
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Mean_Ho_Tot for each clonal rate
slope_Ho <- ggplot(merged_slope_parsed_data, aes(x = slope , y = Mean_Ho_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope", 
    y = "Mean_Ho",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Mean_He_Tot for each clonal rate
slope_He <- ggplot(merged_slope_parsed_data, aes(x = slope , y = Mean_He_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y =  "Mean_He",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )



# Plot Slope vs R for each clonal rate
slope_R <- ggplot(merged_slope_parsed_data, aes(x = slope , y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope", 
    y = "R",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Pareto_beta for each clonal rate
slope_Pareto_beta <- ggplot(merged_slope_parsed_data, aes(x = slope , y =  Pareto_beta_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y =  "Pareto_beta",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )


# combine all plots into one figure two cols
combined_plot_3 <- plot_grid (
  slope_var_fis, slope_mean_fis, slope_Ho, slope_He, slope_R, slope_Pareto_beta,
  labels = c("a", "b", "c", "d", "e", "f"), ncol = 2
)

combined_plot_3

ggsave(
  filename = "combined_slope_population_genetic_indices.png",  # Save the combined plot
  plot = combined_plot_3,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )
```


# Plotting the relationship between exp(slope) and various population genetic indices

```{r, echo=FALSE}
# Plot exp(Slope) vs Var_FIS for each clonal rate
exp_slope_var_fis <- ggplot(merged_slope_parsed_data, aes(x = exp(slope) , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "exp(Slope)",
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

#Plot exp(Slope) vs Mean_FIS for each clonal rate
exp_slope_mean_fis <- ggplot(merged_slope_parsed_data, aes(x = exp(slope), y = Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "exp(Slope)",
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot exp(Slope) vs Mean_Ho_Tot for each clonal rate
exp_slope_Ho <- ggplot(merged_slope_parsed_data, aes(x = exp(slope) , y = Mean_Ho_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "exp(Slope)",
    y = "Mean_Ho",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot exp(Slope) vs Mean_He_Tot for each clonal rate
exp_slope_He <- ggplot(merged_slope_parsed_data, aes(x = exp(slope) , y = Mean_He_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "exp(Slope)",
    y = "Mean_He",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot exp(Slope) vs R for each clonal rate
exp_slope_R <- ggplot(merged_slope_parsed_data, aes(x = exp(slope) , y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "exp(Slope)",
    y = "R",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot exp(Slope) vs Pareto_beta for each clonal rate
exp_slope_Pareto_beta <- ggplot(merged_slope_parsed_data, aes(x = exp(slope) , y = Pareto_beta_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "exp(Slope)",
    y = "Pareto_beta",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# combine all plots into one figure two cols
combined_plot_5 <- plot_grid (
  exp_slope_var_fis, exp_slope_mean_fis, exp_slope_Ho, exp_slope_He, exp_slope_R, exp_slope_Pareto_beta,
  labels = c("a", "b", "c", "d", "e", "f"), ncol = 2
)

#save the combined plot
ggsave(
  filename = "combined_exp_slope_population_genetic_indices.png",  # Save the combined plot
  plot = combined_plot_5,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color


```

# Plotting the relationship between slope  and various population genetic indices without clone rate = 0.999

```{r, echo=FALSE}
# Filter out clonal rate 0.999 for clearer visualization
merged_slope_parsed_data_filtered <- merged_slope_parsed_data %>%
  filter(clonalrate != "0.999")

# Plot Slope vs Var_FIS for each clonal rate (excluding 0.999)
slope_var_fis_filtered <- ggplot(merged_slope_parsed_data_filtered, aes(x = slope , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot Slope vs Mean_FIS for each clonal rate (excluding 0.999)
slope_mean_fis_filtered <- ggplot(merged_slope_parsed_data_filtered, aes(x = slope, y = Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Mean_Ho_Tot for each clonal rate (excluding 0.999)
slope_Ho_filtered <- ggplot(merged_slope_parsed_data_filtered, aes(x = slope , y = Mean_Ho_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Mean_Ho",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Mean_He_Tot for each clonal rate (excluding 0.999)
slope_He_filtered <- ggplot(merged_slope_parsed_data_filtered, aes(x = slope , y = Mean_He_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Mean_He",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs R for each clonal rate (excluding 0.999)
slope_R_filtered <- ggplot(merged_slope_parsed_data_filtered, aes(x = slope , y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "R",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )


# Plot Slope vs Pareto_beta for each clonal rate (excluding 0.999)
slope_Pareto_beta_filtered <- ggplot(merged_slope_parsed_data_filtered, aes(x = slope , y = Pareto_beta_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Pareto_beta",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# combine all plots into one figure two cols
combined_plot_6 <- plot_grid (
  slope_var_fis_filtered, slope_mean_fis_filtered, slope_Ho_filtered, slope_He_filtered, slope_R_filtered, slope_Pareto_beta_filtered,
  labels = c("a", "b", "c", "d", "e", "f"), ncol = 2
)

#save the combined plot
ggsave(
  filename = "combined_slope_population_genetic_indices_filtered.png",  # Save the combined plot
  plot = combined_plot_6,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color

```


# Plotting the relationship between slope  and various population genetic indices till vlone rate ≤ 0.9
```{r, echo=FALSE}
# Filter out clonal rates greater than 0.9 for clearer visualization
merged_slope_parsed_data_filtered_0_9 <- merged_slope_parsed_data %>%
  filter(clonalrate <= "0.9")
# Plot Slope vs Var_FIS for each clonal rate (excluding clonal rates > 0.9)
slope_var_fis_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_0_9, aes(x = slope , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot Slope vs Mean_FIS for each clonal rate (excluding clonal rates > 0.9)
slope_mean_fis_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_0_9, aes(x = slope, y = Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot Slope vs Mean_Ho_Tot for each clonal rate (excluding clonal rates > 0.9)
slope_Ho_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_0_9, aes(x = slope , y = Mean_Ho_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Mean_Ho",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot Slope vs Mean_He_Tot for each clonal rate (excluding clonal rates > 0.9)
slope_He_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_0_9, aes(x = slope , y = Mean_He_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Mean_He",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot Slope vs R for each clonal rate (excluding clonal rates > 0.9)
slope_R_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_0_9, aes(x = slope , y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "R",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot Slope vs Pareto_beta for each clonal rate (excluding clonal rates > 0.9)
slope_Pareto_beta_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_0_9, aes(x = slope , y = Pareto_beta_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "Slope",
    y = "Pareto_beta",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# combine all plots into one figure two cols
combined_plot_7 <- plot_grid (
  slope_var_fis_filtered_0_9, slope_mean_fis_filtered_0_9, slope_Ho_filtered_0_9, slope_He_filtered_0_9, slope_R_filtered_0_9, slope_Pareto_beta_filtered_0_9,
  labels = c("a", "b", "c", "d", "e", "f"), ncol = 2
)

#save the combined plot
ggsave(
  filename = "combined_slope_population_genetic_indices_filtered_0_9.png",  # Save the combined plot
  plot = combined_plot_7,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color
```



# Plotting the relationship between R^2 and various population genetic indices

```{r, echo=FALSE}
# Plot R^2 vs Var_FIS for each clonal rate
R2_var_fis <- ggplot(merged_slope_parsed_data, aes(x = R2 , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_FIS for each clonal rate
R2_mean_fis <- ggplot(merged_slope_parsed_data, aes(x = R2 , y = Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_Ho_Tot for each clonal rate
R2_Ho <- ggplot(merged_slope_parsed_data, aes(x = R2 , y = Mean_Ho_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_Ho",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_He_Tot for each clonal rate
R2_He <- ggplot(merged_slope_parsed_data, aes(x = R2 , y = Mean_He_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_He",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs R for each clonal rate
R2_R <- ggplot(merged_slope_parsed_data, aes(x = R2 , y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "R",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Pareto_beta for each clonal rate
R2_Pareto_beta <- ggplot(merged_slope_parsed_data, aes(x = R2 , y = Pareto_beta_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Pareto_beta",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# combine all plots into one figure two cols
combined_plot_4 <- plot_grid (
  R2_var_fis, R2_mean_fis, R2_Ho, R2_He, R2_R, R2_Pareto_beta,
  labels = c("a", "b", "c", "d", "e", "f"), ncol = 2
)

# save the combined plot
ggsave(
  filename = "combined_R2_population_genetic_indices.png",  # Save the combined plot
  plot = combined_plot_4,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )

combined_plot_4
```


## scatter plot for R^2 vs Var_FIS and R^2 vs Mean_FIS for c= 0.999 , 0.9999
```{r, echo=FALSE}
# Filter the merged_slope_parsed_data for clonal rates 0.999 and 0.9999

merged_slope_parsed_df_extreme <- merged_slope_parsed_data %>%
  filter(clonalrate %in% c("0.99","0.999", "0.9999"))

R2_var_fis_extreme  <- ggplot(merged_slope_parsed_df_extreme, aes(x = R2 , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_FIS for each clonal rate
R2_mean_fis_extreme  <- ggplot(merged_slope_parsed_df_extreme, aes(x = R2 , y = Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )




# combine the plots into one figure two cols
combined_plot_8 <- plot_grid (
  R2_var_fis_0.999, R2_mean_fis_0.999, R2_var_fis_0.9999, R2_mean_fis_0.9999,
  labels = c("a", "b", "c", "d"), ncol = 2
)

# save
ggsave(
  filename = "combined_R2_var_fis_mean_fis_0.999_0.9999.png",  # Save the combined plot
  plot = combined_plot_8,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 20,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color
```

# Scatter plot for 4 replicates within specific clonerate
#0.1
```{r, echo=FALSE}
# Take 4 randome replicates from the clonal rate 0.1

set.seed(123)  
random_reps_0.1 <- parsed_data %>%
  filter(clonalrate == "0.1") %>%
  pull(Replicate) %>%
  unique() %>%
  sample(4)

# View the selected replicates
#random_reps_0.1

df_0.1 <- parsed_data %>% filter(clonalrate == "0.1") %>% filter (Generation == 6000)  %>% filter(Replicate %in% random_reps_0.1)

# create scatter plot for each replicate
# create list to store plots

scatter_plots_0.1 <- list()

for (replicate in random_reps_0.1) {
  p <- ggplot(df_0.1 %>% filter(Replicate == replicate), aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Adds regression line only
    labs(
      x = "Clone age",
      y = "Number of individuals (log scale)",
      title = paste("Replicate", replicate)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  scatter_plots_0.1[[paste0("scatter_plot_", replicate)]] <- p
}

# combine the plots into one
combined_scatter_plots_0.1 <- plot_grid(plotlist = scatter_plots_0.1, ncol = 2)

combined_scatter_plots_0.1

# Save the combined scatter plots
ggsave(
  filename = "scatter_plots_0.1.png",  
  plot = combined_scatter_plots_0.1,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
```

#0.7
```{r, echo=FALSE}
# Take 4 randome replicates from the clonal rate 

set.seed(123)  
random_reps_0.7 <- parsed_data %>%
  filter(clonalrate == "0.7") %>%
  pull(Replicate) %>%
  unique() %>%
  sample(4)

# View the selected replicates
df_0.7 <- parsed_data %>% filter(clonalrate == "0.7") %>% filter (Generation == 6000) %>% filter(Replicate %in% random_reps_0.7)

# create scatter plot for each replicate
# create list to store plots

scatter_plots_0.7 <- list()

for (replicate in random_reps_0.7) {
  p <- ggplot(df_0.7 %>% filter(Replicate == replicate), aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Adds regression line only
    labs(
      x = "Clone age",
      y = "Number of individuals (log scale)",
      title = paste("Replicate", replicate)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  scatter_plots_0.7[[paste0("scatter_plot_", replicate)]] <- p
}

# combine the plots into one
combined_scatter_plots_0.7 <- plot_grid(plotlist = scatter_plots_0.7, ncol = 2)

combined_scatter_plots_0.7

# Save the combined scatter plots
ggsave(
  filename = "scatter_plots_0.7.png",  
  plot = combined_scatter_plots_0.7,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
```

#0.9
```{r, echo=FALSE}
# Take 4 randome replicates from the clonal rate 

set.seed(123)  
random_reps_0.9 <- parsed_data %>%
  filter(clonalrate == "0.9") %>%
  pull(Replicate) %>%
  unique() %>%
  sample(4)

# View the selected replicates
df_0.9 <- parsed_data %>% filter(clonalrate == "0.9") %>% filter (Generation == 6000) %>% filter(Replicate %in% random_reps_0.9)

# create scatter plot for each replicate
# create list to store plots

scatter_plots_0.9 <- list()

for (replicate in random_reps_0.9) {
  p <- ggplot(df_0.9 %>% filter(Replicate == replicate), aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Adds regression line only
    labs(
      x = "Clone age",
      y = "Number of individuals (log scale)",
      title = paste("Replicate", replicate)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  scatter_plots_0.9[[paste0("scatter_plot_", replicate)]] <- p
}

# combine the plots into one
combined_scatter_plots_0.9 <- plot_grid(plotlist = scatter_plots_0.9, ncol = 2)

combined_scatter_plots_0.9

# Save the combined scatter plots
ggsave(
  filename = "scatter_plots_0.9.png",  
  plot = combined_scatter_plots_0.9,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
```

#0.999
```{r, echo=FALSE}
# Take 4 randome replicates from the clonal rate XX
set.seed(123)
random_reps_0.999 <- parsed_data %>%
  filter(clonalrate == "0.999") %>%
  pull(Replicate) %>%
  unique() %>%
  sample(4)
# View the selected replicates
df_0.999 <- parsed_data %>% filter(clonalrate == "0.999") %>% filter (Generation == 6000) %>% filter(Replicate %in% random_reps_0.999)
# create scatter plot for each replicate
# create list to store plots
scatter_plots_0.999 <- list()

for (replicate in random_reps_0.999) {
  p <- ggplot(df_0.999 %>% filter(Replicate == replicate), aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Adds regression line only
    labs(
      x = "Clone age",
      y = "Number of individuals (log scale)",
      title = paste("Replicate", replicate)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  scatter_plots_0.999[[paste0("scatter_plot_", replicate)]] <- p
}
# combine the plots into one
combined_scatter_plots_0.999 <- plot_grid(plotlist = scatter_plots_0.999, ncol = 2)
combined_scatter_plots_0.999
# Save the combined scatter plots
ggsave(
  filename = "scatter_plots_0.999.png",  
  plot = combined_scatter_plots_0.999,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
```

# Plot the relationship between slope and R^2 and R (cross-shaped scatter)

```{r, echo=FALSE}
# Calculate mean and standard deviation for R^2 and R_Tot grouped by clonalrate
error_data <- merged_slope_parsed_data %>%
  group_by(clonalrate) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_R_Tot = mean(R_Tot),
    sd_R_Tot = sd(R_Tot),
    Pareto_beta_Tot = mean(Pareto_beta_Tot),
    sd_R_Tot = sd(R_Tot)
  )

r2_r <- ggplot(error_data, aes(x = mean_R2, y = mean_R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 4, shape = 3) +  # cross shape
  geom_errorbar(aes(ymin = mean_R_Tot - sd_R_Tot, ymax = mean_R_Tot + sd_R_Tot), width = 0.05) +  # vertical bars
  geom_errorbarh(aes(xmin = mean_R2 - sd_R2, xmax = mean_R2 + sd_R2), height = 0.05) +  # horizontal bars
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(x = "R^2", y = "R", color = "Clonal Rate") +
  theme_minimal() +
  theme(legend.position = "right")

r2_beta <- ggplot(error_data, aes(x = mean_R2, y = Pareto_beta_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 4, shape = 3) +  # cross shape
  geom_errorbar(aes(ymin = Pareto_beta_Tot - sd_R_Tot, ymax = Pareto_beta_Tot + sd_R_Tot), width = 0.05) +  # vertical bars
  geom_errorbarh(aes(xmin = mean_R2 - sd_R2, xmax = mean_R2 + sd_R2), height = 0.05) +  # horizontal bars
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(x = "R^2", y = "Pareto_beta", color = "Clonal Rate") +
  theme_minimal() +
  theme(legend.position = "right")


#save
ggsave(
  filename = "R2_R_cross_scatter.png",  # Save the cross-shaped scatter plot
  plot = r2_r,  # The cross-shaped scatter plot object
  width = 20,  # Width of the saved image
  height = 15,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color

#save
ggsave(
  filename = "R2_Pareto_beta_cross_scatter.png",  # Save the cross-shaped scatter plot
  plot = r2_beta,  # The cross-shaped scatter plot object
  width = 20,  # Width of the saved image
  height = 15,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color

```









# Plotting the relationship between R^2 and various population genetic indices ≤0.9
```{r, echo=FALSE}
# Filter out clonal rates greater than 0.9 for clearer visualization
merged_slope_parsed_data_filtered_R2_0_9 <- merged_slope_parsed_data %>%
  filter(clonalrate <= "0.9")
# Plot R^2 vs Var_FIS for each clonal rate (excluding clonal rates > 0.9)
R2_var_fis_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_R2_0_9, aes(x = R2 , y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Var(FIS)",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_FIS for each clonal rate (excluding clonal rates > 0.9)
R2_mean_fis_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_R2_0_9, aes(x = R2 , y = Mean_FIS_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_FIS",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_Ho_Tot for each clonal rate (excluding clonal rates > 0.9)
R2_Ho_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_R2_0_9, aes(x = R2 , y = Mean_Ho_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_Ho",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Mean_He_Tot for each clonal rate (excluding clonal rates > 0.9)
R2_He_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_R2_0_9, aes(x = R2 , y = Mean_He_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Mean_He",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs R for each clonal rate (excluding clonal rates > 0.9)
R2_R_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_R2_0_9, aes(x = R2 , y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "R",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# Plot R^2 vs Pareto_beta for each clonal rate (excluding clonal rates > 0.9)
R2_Pareto_beta_filtered_0_9 <- ggplot(merged_slope_parsed_data_filtered_R2_0_9, aes(x = R2 , y = Pareto_beta_Tot , color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    x = "R^2",
    y = "Pareto_beta",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )
# combine all plots into one figure two cols
combined_plot_8 <- plot_grid (
  R2_var_fis_filtered_0_9, R2_mean_fis_filtered_0_9, R2_Ho_filtered_0_9, R2_He_filtered_0_9, R2_R_filtered_0_9, R2_Pareto_beta_filtered_0_9,
  labels = c("a", "b", "c", "d", "e", "f"), ncol = 2
)
# save the combined plot
ggsave(
  filename = "combined_R2_population_genetic_indices_filtered_0_9.png",  # Save the combined plot
  plot = combined_plot_8,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" ) # Background color

```


## Ages per each clone rate 

```{r, echo=FALSE}

# Summarize data to count distinct ages per clonal rate and replicate
age_counts <- parsed_data %>%
  group_by(clonalrate, Replicate) %>%
  summarize(Age_Count = n_distinct(Age), groups = "drop")  # Count unique ages per replicate within each clonal rate

# Convert clonal rate to a factor, so it appears in order on the x-axis
age_counts$clonalrate <- factor(age_counts$clonalrate, levels = sort(unique(age_counts$clonalrate)))

# Violin-Box Plot
ggplot(age_counts, aes(x = clonalrate, y = Age_Count, fill = clonalrate)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "#104E8B", size = 0.8) +  # Violin plot with transparency
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA, color = "black") +  # Overlay box plot
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +  # Add median point
  scale_fill_brewer(palette = "Blues") +  # Color palette for fill
  labs(
    title = "Distribution of Age Group Counts by Clonal Rate",
    x = "Clonal Rate (C value)",
    y = "Number of Age Groups in Replicate"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  theme(
    legend.position = "none",  # Remove legend for clarity
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  )

```


## oldest age per clonalrate

```{r, echo=FALSE}

# Summarize data to find the oldest age per clonal rate and replicate
oldest_ages <- parsed_data %>%
  group_by(clonalrate, Replicate) %>%
  summarize(Max_Age = max(Age), .groups = "drop")  # Get the max age per replicate within each clonal rate

# Convert clonal rate to a factor for ordered plotting
oldest_ages$clonalrate <- factor(oldest_ages$clonalrate, levels = sort(unique(oldest_ages$clonalrate)))

#Plot the oldest age per clonal rate using a box plot
ggplot(oldest_ages, aes(x = clonalrate, y = log(Max_Age), fill = clonalrate)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.color = "red", alpha = 0.7) +  # Box plot for oldest ages
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "darkblue") +  # Add jittered points for individual data
  labs(
    title = "Oldest Age Observed per Clonal Rate",
    x = "Clonal Rate (C value)",
    y = "Oldest Age in Replicate"
  ) +
  scale_fill_brewer(palette = "Blues") +  # Color palette
  theme_minimal(base_size = 15) +  # Minimal theme with larger font size
  theme(
    legend.position = "none",  # Remove legend for clarity
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  )

#max age at clone rate = 0.9
# Find the maximum age at clonal rate 0.9
max_age_09 <- parsed_data %>%
  filter(clonalrate == "0.9") %>%
  summarize(Max_Age = max(Age))

# Add mean age of clones on same plot diff color or another plot
```



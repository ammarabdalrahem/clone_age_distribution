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

required_packages <- c("knitr","rmarkdown","ggridges","data.table","ggplot2","gganimate","dplyr","tidyr","cowplot", "reshape2")
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


# Import the data 
file_path_S1 <- "2024-08-30-18h_04minsimulationS1.txt"
data_S1 <- fread(file_path_S1, header = TRUE, sep = "\t", fill = TRUE)

# Import data for another simulation which clonerate =  0.98 and 0.9999
file_path_0.98 <- "2024-10-28-14h_23minsimulationS1.txt"
data_0.98 <- fread(file_path_0.98, header = TRUE, sep = "\t", fill = TRUE)

# Import data for another simulation which clonerate =  0.95, 0.96, and 0.97
file_path_0.95 <- "2024-11-04-17h_21minsimulationS1.txt"
data_0.95 <- fread(file_path_0.95, header = TRUE, sep = "\t", fill = TRUE)

# Without and low mutation u=0


file_path_no_mutation<- "2024-11-13-10h_14minsimulationS1.txt"
data_no_mutation <- fread(file_path_no_mutation, header = TRUE, sep = "\t", fill = TRUE)

#u=10^-6 

file_path_mutation <- "2024-11-13-09h_40minsimulationS1.txt"
data_mutation <- fread(file_path_mutation, header = TRUE, sep = "\t", fill = TRUE)




#data <- data_mutation 

# Combine two df
data <- rbind(data_S1,data_0.98,data_0.95)



## Data polish
# Remove last three columns 
data <- data[,-16:-18]

# Change header names
colnames(data) <- c("Replicate", "Mutation_Rate", "clonalrate", "Number_alleles","Number_fixed_loci", "Mean_He","Mean_Ho","Mean_FIS","Var_FIS","Mean_r_bar_D", "SD_r_bar_D", "Number_Genotypes", "R", "Pareto_beta","List_distribution_gen_clonal_genotypes")


head(data)
```




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
data <- data_dif_mutation %>% filter(Nb_alleles_tot == "4") %>% filter(Mutation_Rate == 0.000001) %>% filter (Generation %% 10 == 0)
head (data)

#data_dif_mutation_m4_N4$Mean_FIS_Tot <- as.numeric(data_dif_mutation_m4_N4$Mean_FIS_Tot)
```


##??
```{r, echo=FALSE}

# Prepare data
data_prepared <- data_dif_mutation_m4_N4 %>%
  # Convert to numeric
  mutate(Mean_FIS_Tot = as.numeric(Mean_He_Tot)) %>%
  # Remove exact duplicates
  distinct(Generation, clonalrate, Mean_He_Tot, .keep_all = TRUE) # duplicated rows didn't removed??


# Create the contour plot
ggplot(data_prepared, aes(x = Generation, y = clonalrate, z = Mean_He_Tot)) +
  geom_contour_filled() +  # Filled contour plot
  geom_contour(color = "black", size = 0.5, alpha = 0.5) +  # Add contour lines
  scale_fill_viridis_d(name = "Mean FIS ") +  # Color palette
  labs(
    title = " Mean FIS ",
    x = "Generation",
    y = "Clonal Rate"
  ) +
  theme_minimal() +
  guides(fill = guide_colorsteps(barwidth = 1, barheight = 10))


ggplot(data_prepared, aes(x = log(Generation), y = clonalrate)) +
  geom_point(aes(size = 10, color = Mean_He_Tot), alpha = 0.2) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(1, 10)) +
  labs(
    title = "Mean FIS Visualization",
    x = "Generation",
    y = "Clonal Rate",
    size = "Number fixed loci",
    color = "Mean_He_Tot"
  ) +
  theme_minimal()

```


## The relation between population genetic indices and clonerate 
# As already shown in Stoeckel et al. (2021), this relationship, with an initial population size of 10^4 Figure S1a.
# Reproduce the analysis again 


```{r, echo=FALSE}

# Create the boxplot for R
plot_R <- ggplot(data, aes(x = as.factor(clonalrate), y = R_Tot)) +
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
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Create the boxplot for Pareto_beta
plot_Pareto_beta <- ggplot(data, aes(x = as.factor(clonalrate), y = Pareto_beta_Tot)) +
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
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


#Creat the boxplot for Var_FIS
plot_Var_FIS <- ggplot(data, aes(x = as.factor(clonalrate), y = Var_FIS_Tot)) +
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
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Creat the boxplot for Mean_FIS
plot_Mean_FIS <- ggplot(data, aes(x = as.factor(clonalrate), y = Mean_FIS_Tot)) +
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
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Creat the boxplot for Mean_He
plot_Mean_He <- ggplot(data, aes(x = as.factor(clonalrate), y = Mean_He_Tot)) +
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
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

#Creat the boxplot for Mean_Ho
plot_Mean_Ho <- ggplot(data, aes(x = as.factor(clonalrate), y = Mean_Ho_Tot)) +
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
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )




# Combine the plots 
# cowplot package
combined_plot <- plot_grid( plot_Mean_He, plot_Mean_Ho, plot_R , plot_Pareto_beta,  plot_Var_FIS, plot_Mean_FIS, 
                            labels = c("A", "B", "C", "D", "E","F"), ncol = 2)

# Save plot
ggsave(
  filename = "combined_pop_indices_plot.png",  
  plot = combined_plot,           
  width = 25,                  
  height = 20,    
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
    title = "Clone ages by clonerate",
    x = "Clone age (Generations)",
    y = "Clonerates"
  ) +
  theme_ridges() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    
  )

# Save plot
ggsave(
  filename = "Joyplot_clone_ages_by_clonerateplot.png",  
  plot = Joyplot_plot,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
# Check the lowest value in clone age 
#min(df$Clone_Age)


```


## Joyplot limit to 40

```{r, echo=FALSE}
# Joyplot limit to 40
joyplot_40 <- ggplot(df, aes(x = Clone_Age, y = Group, fill = Group)) +
  geom_density_ridges(scale = 3, alpha = 0.7) +
  #scale_fill_viridis_d(option = "C") +
  labs(
    #title = "Clone Ages by Clonerate (Limited to 40)",
    x = "Clone Age (Generations)",
    y = "Clonerates"
  ) +
  xlim(NA ,40) +  # Set x-axis limits to end at 40
  theme_ridges() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    #plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
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

# Plot mean FIS by clonal rate and generation
plot_Mean_FIS_2 <-  ggplot(data, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Mean("*italic(F)[IS]*")")) +
  labs(
    x = "Generation",
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot variance of FIS by clonal rate and generation
plot_Var_FIS_2 <- ggplot(data, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Var_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Var("*italic(F)[IS]*")")) +
  labs(
    x = "Generation",
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


 # Plot mean He by clonal rate and generation
plot_Mean_He_2 <- ggplot(data, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_He_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[e])) +
  labs(
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot mean Ho by clonal rate and generation
plot_Mean_Ho_2 <- ggplot(data, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_Ho_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[o])) +
  labs(
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Plot R by clonal rate and generation
plot_R_2 <-  ggplot(data, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = R_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "R") +
  labs(
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot Pareto beta by clonal rate and generation
plot_Pareto_beta_2 <- ggplot(data, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Pareto_beta_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(β) * " Pareto")) +
  labs(
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )



# Combine the plots 
combined_plot_2 <- plot_grid(plot_Mean_He_2, plot_Mean_Ho_2, plot_R_2 , plot_Pareto_beta_2,  plot_Var_FIS_2, plot_Mean_FIS_2, 
                            labels = c("A", "B", "C", "D", "E","F"), ncol = 2)

# Save plot
ggsave(
  filename = "combined_pop_indices_plot_2.png",  
  plot = combined_plot_2,           
  width = 25,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)

combined_plot_2

```


## R^2 per each replicate

```{r, echo=FALSE}
# Initialize a data frame to store the results
r2_results <- data.frame(clonalrate = character(), Replicate = integer(), R2 = numeric(), stringsAsFactors = FALSE)

# Get unique clonal rates
unique_clonal_rates <- unique(parsed_data$clonalrate)

# Loop over each clonal rate
for (clonal_rate in unique_clonal_rates) {
  # Filter data for the current clonal rate
  parsed_data_clonal <- parsed_data %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log-transform Count

  # Loop over unique replicates for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate
    replicate_data <- parsed_data_clonal %>%
      filter(Replicate == replicate)

    # Calculate R^2 if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data)  # Fit linear model
      r2 <- summary(lm_fit)$r.squared  # Extract R^2 from the model summary
    } else {
      r2 <- NA  # If not enough data, assign NA
    }

    # Store the results
    r2_results <- rbind(r2_results, data.frame(clonalrate = clonal_rate, Replicate = replicate, R2 = r2))
  }
}

# Remove rows where R2 is NA
r2_results <- r2_results %>% filter(!is.na(R2))

# Convert 'clonalrate' to a factor
r2_results$clonalrate <- as.factor(r2_results$clonalrate)

# Display the results
print(r2_results)


# Plot 
R_2_plot<- ggplot(r2_results, aes(x = clonalrate, y = R2)) +
  geom_violin(fill = "#A4D3EE", color = "#104E8B", alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.color = "red", outlier.shape = 16) +
  labs(
    title = "Distribution of R^2 Values by Clonal Rate",
    x = "Clonal Rate",
    y = "R^2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save plot
ggsave(
  filename = "R_2_plot.png",  
  plot = R_2_plot,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 300                       
)

# Without Age zero 

# Initialize a data frame to store the results
r2_results_Zero <- data.frame(clonalrate = character(), Replicate = integer(), R2 = numeric(), stringsAsFactors = FALSE)

# Get unique clonal rates
unique_clonal_rates <- unique(parsed_data$clonalrate)

# Loop over each clonal rate
for (clonal_rate in unique_clonal_rates) {
  # Filter data for the current clonal rate
  parsed_data_clonal <- parsed_data %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log-transform Count

  # Loop over unique replicates for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate and exclude Age == 0
    replicate_data <- parsed_data_clonal %>%
      filter(Replicate == replicate, Age > 0)  # Exclude Age == 0

    # Calculate R^2 if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data)  # Fit linear model
      r2 <- summary(lm_fit)$r.squared  # Extract R^2 from the model summary
    } else {
      r2 <- NA  # If not enough data, assign NA
    }

    # Store the results
    r2_results_Zero <- rbind(r2_results_Zero, data.frame(clonalrate = clonal_rate, Replicate = replicate, R2 = r2))
  }
}

# Remove rows where R2 is NA
r2_results_Zero <- r2_results_Zero %>% filter(!is.na(R2))

# Convert 'clonalrate' to a factor
r2_results_Zero$clonalrate <- as.factor(r2_results_Zero$clonalrate)

# Display the results
print(r2_results)


# Plot 
ggplot(r2_results_Zero, aes(x = clonalrate, y = R2)) +
  geom_violin(fill = "#A4D3EE", color = "#104E8B", alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.color = "red", outlier.shape = 16) +
  labs(
    title = "Distribution of R^2 Values by Clonal Rate (Excluding Age 0)",
    x = "Clonal Rate",
    y = "R^2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

```

```{r, echo=FALSE}


```


















## Age zero R and C value
```{r, echo=FALSE}
# Boxplot
ggplot(data, aes(x = as.factor(clonalrate), y = R_Tot)) +
  geom_boxplot(fill = "#FFDDC1", color = "#FF5733", outlier.color = "red", outlier.shape = 16) +
  labs(
    title = "Boxplot of R ratio by clonal rate",
    x = "Clonal rate (c value)",
    y = "R"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Filter parsed_data at Age 0
age_zero_data <- parsed_data %>%
  filter(Age == 0)

# C value and Age zero
# boxplot
ggplot(age_zero_data, aes(x = as.factor(clonalrate), y = (Count/1000))) +
  geom_boxplot(fill = "#A4D3EE", color = "#104E8B", outlier.color = "red", outlier.shape = 16) +
  labs(
    title = "Boxplot of individuals at age zero by clonal rate",
    x = "Clonal rate (c value)",
    y = "Number of individuals at age zero"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Calculate the summary table (C,mean R, mean age zero)
summary_table <- parsed_data %>%
  group_by(clonalrate) %>%
  summarise(
    Mean_R = mean(R_Tot, na.rm = TRUE),  # Mean of R for each clonalrate
    Mean_Age_Zero = mean(Count[Age == 0], na.rm = TRUE)  # Mean count for Age zero for each clonalrate
  ) %>%
  ungroup()

# View the result
print(summary_table)

plot (summary_table$Mean_R~summary_table$Mean_Age_Zero)

# Scatter plot 
ggplot(summary_table, aes(x = Mean_Age_Zero, y = Mean_R)) +
  geom_point(size = 3, color = "blue", alpha = 0.7) +  # Points with transparency
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Optional trend line with confidence interval
  labs(
    title = "Relationship between mean count at age zero and mean R",
    x = "Mean count at age zero",
    y = "Mean R"
  ) +
  theme_minimal()



```

## Extract specific clonerate

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
parsed_data_0.6 <- parsed_data %>% filter(clonalrate == "0.6")


ggplot(parsed_data_0.6, aes(x = Age, y = Count, color = as.factor(clonalrate))) +
  geom_point(alpha = 0.7) +  # transparency
  scale_y_log10() +  # Set y-axis to log scale
  scale_x_log10() +  # Set x-axis to log scale
  labs(
    title = "Number of individuals vs. Clone age by clonal rate",
    x = "Age of clones (log scale)",
    y = "Number of individuals (log scale)",
    color = "Clonal rate"
  ) +
  theme_minimal() +  # Minimal theme for clarity
  theme(
    legend.position = "top",  # Legend at the top
    legend.direction = "horizontal"  # Horizontal legend layout
  ) +
  guides(color = guide_legend(nrow = 1))  # Arrange legend in one row



df_06 <- df %>% filter(Group == "0.6")

# Calculate mean and confidence intervals
summary_data_06 <- df_06 %>%
  group_by(Clone_Age) %>%
  summarise(
    mean_count = n(),  # Count the number of individuals for each age
    se = sd(Clone_Age) / sqrt(n()),  # Standard Error
    lower_ci = mean_count - 1.96 * se,  # 95% CI lower bound
    upper_ci = mean_count + 1.96 * se   # 95% CI upper bound
  )

# Plot
ggplot(df_06, aes(x = Clone_Age)) +
  geom_histogram(aes(y = after_stat(count)), bins = 30, fill = "lightgray", color = "black", alpha = 0.5) +  # Histogram for individual counts
  geom_line(data = summary_data_06, aes(x = Clone_Age, y = mean_count), color = "blue", linewidth = 1) +  # Mean line
  geom_ribbon(data = summary_data_06, aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "blue") +  # CI
  labs(
    title = "Clone age distribution for clonal rate 0.6",
    x = "Clone age",
    y = "Number of individuals"
  ) +
  theme_minimal()

```



## Slope per each replicate

```{r, echo=FALSE}

# Initialize a data frame to store the slopes
slope_results <- data.frame(clonalrate = character(), Replicate = integer(), Slope = numeric(), stringsAsFactors = FALSE)

# Get unique clonal rates
unique_clonal_rates <- unique(parsed_data$clonalrate)

# Loop over for each clonal rates
for (clonal_rate in unique_clonal_rates) {
  # Filter data for the current clonal rate
  parsed_data_clonal <- parsed_data %>% # To take each clone rate
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log transform Count

  # Loop over unique replicates for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate
    replicate_data <- parsed_data_clonal %>% # To take each replicate within the smae clonerate 
      filter(Replicate == replicate)

    # Fit a linear model if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data) # Run linear regression model to see how log_Count changes with Age
      slope <- coef(lm_fit)["Age"]  # Extract the slope # how much log_Count changes with each unit of Age
    } else {
      slope <- NA  # If not enough data, assign NA
    }

    # Store the results
    slope_results <- rbind(slope_results, data.frame(clonalrate = clonal_rate, Replicate = replicate, Slope = slope))
  }
}

# Remove rows where Slope is NA
slope_results <- slope_results %>% filter(!is.na(Slope))

# Slope_results df 'Slope' and 'Clonal_Rate'
slope_results$clonalrate <- as.factor(slope_results$clonalrate) 

ggplot(slope_results, aes(x = clonalrate, y = Slope)) +
  geom_violin(fill = "#A4D3EE", color = "#104E8B", alpha = 0.7) +  # Create the violin plot
  geom_boxplot(width = 0.1, color = "black", outlier.color = "red", outlier.shape = 16) +  # Add a box plot on top
  labs(
    title = "Violin Plot of Slopes by Clonal Rate",
    x = "Clonal Rate",
    y = "Slope"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# Create a joyplot (ridgeline plot)
ggplot(slope_results, aes(x = exp(Slope), y = clonalrate, fill = clonalrate)) +
  geom_density_ridges(scale = 0.9, alpha = 0.6, rel_min_height = 0.01) +  # Draw ridgelines
  labs(
    title = "Joyplot of Slopes by Clonal Rate",
    x = "Slope",
    y = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove the legend
    axis.text.y = element_text(angle = 0)  # Keep y-axis text horizontal
  )


# Print the slope results
#print(slope_results)

# Without age zero

# Initialize a data frame to store the slopes
slope_results_Zero <- data.frame(clonalrate = character(), Replicate = integer(), Slope = numeric(), stringsAsFactors = FALSE)

# Loop over each clonal rate
for (clonal_rate in unique_clonal_rates) {
  # Filter data for the current clonal rate
  parsed_data_clonal <- parsed_data %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log transform Count
  
  # Loop over unique replicates for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate and exclude Age == 0
    replicate_data <- parsed_data_clonal %>%
      filter(Replicate == replicate, Age > 0)  # Exclude Age == 0
    
    # Fit a linear model if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data)  # Run linear regression model
      slope <- coef(lm_fit)["Age"]  # Extract the slope
    } else {
      slope <- NA  # If not enough data, assign NA
    }

    # Store the results
    slope_results_Zero <- rbind(slope_results_Zero, data.frame(clonalrate = clonal_rate, Replicate = replicate, Slope = slope))
  }
}

# Remove rows where Slope is NA
slope_results_Zero <- slope_results_Zero %>% filter(!is.na(Slope))

# Convert 'clonalrate' column to factor
slope_results_Zero$clonalrate <- as.factor(slope_results_Zero$clonalrate) 

# Plot the slopes using a violin plot with a boxplot overlay
ggplot(slope_results_Zero, aes(x = clonalrate, y = Slope)) +
  geom_violin(fill = "#A4D3EE", color = "#104E8B", alpha = 0.7) +  # Create the violin plot
  geom_boxplot(width = 0.1, color = "black", outlier.color = "red", outlier.shape = 16) +  # Add a box plot on top
  labs(
    title = "Violin Plot of Slopes by Clonal Rate (Excluding Age 0)",
    x = "Clonal Rate",
    y = "Slope"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )





```








```{r}
# Initialize an empty data frame to store residuals for all clone rates
residuals_long <- data.frame(Clonal_Rate = character(), Replicate = integer(), Age = numeric(), Residual = numeric(), stringsAsFactors = FALSE)

# Loop through each unique clonal rate
for (clonal_rate in unique(parsed_data$clonalrate)) {
  # Filter the data for the current clonal rate
  parsed_data_clonal <- parsed_data %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log-transform Count for the linear model
  
  # Loop through each replicate for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate
    replicate_data <- parsed_data_clonal %>%
      filter(Replicate == replicate)
    
    # Fit a linear model if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data)
      residuals <- residuals(lm_fit)  # Calculate residuals
      
      # Append the residuals to the data frame
      residuals_long <- rbind(
        residuals_long,
        data.frame(
          Clonal_Rate = clonal_rate,
          Replicate = replicate,
          Age = replicate_data$Age,
          Residual = residuals
        )
      )
    }
  }
}

head(residuals_long)

residuals_long_c01 <- filter(residuals_long , Clonal_Rate == 0.2)

# Plot the residuals
ggplot(residuals_long_c01, aes(x = Age, y = Residual, color = as.factor(Clonal_Rate))) +
  geom_point(alpha = 0.6, size = 3) +
  labs(
    title = "Residuals by clone age and clonal rate",
    x = "Clone Age",
    y = "Residuals"
  ) +
  xlim(NA ,10) +  # Set x-axis limits to end at 50
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Without Zero

# Initialize an empty data frame to store residuals for all clone rates
residuals_long_zero <- data.frame(Clonal_Rate = character(), Replicate = integer(), Age = numeric(), Residual = numeric(), stringsAsFactors = FALSE)

# Loop through each unique clonal rate
for (clonal_rate in unique(parsed_data$clonalrate)) {
  # Filter the data for the current clonal rate
  parsed_data_clonal <- parsed_data %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log-transform Count for the linear model
  
  # Loop over unique replicates for the current clonal rate
  for (replicate in unique(parsed_data_clonal$Replicate)) {
    # Filter data for the current replicate and exclude Age == 0
    replicate_data <- parsed_data_clonal %>%
      filter(Replicate == replicate, Age > 0)  # Exclude Age == 0
    
    # Fit a linear model if there are enough data points
    if (nrow(replicate_data) > 1) {  # Ensure there are at least 2 points
      lm_fit <- lm(log_Count ~ Age, data = replicate_data)
      residuals <- residuals(lm_fit)  # Calculate residuals
      
      # Append the residuals to the data frame
      residuals_long_zero <- rbind(
        residuals_long_zero,
        data.frame(
          Clonal_Rate = clonal_rate,
          Replicate = replicate,
          Age = replicate_data$Age,
          Residual = residuals
        )
      )
    }
  }
}

head(residuals_long_zero)

residuals_long_c01 <- filter(residuals_long_zero , Clonal_Rate == 0.2)

# Plot the residuals
ggplot(residuals_long_c01, aes(x = Age, y = Residual, color = as.factor(Clonal_Rate))) +
  geom_point(alpha = 0.6, size= 3) +
  labs(
    title = "Residuals by clone age and clonal rate (Excluding Age 0)",
    x = "Clone Age",
    y = "Residuals"
  ) +
  xlim(NA ,10) +  # Set x-axis limits to end at 50
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )




```





```{r}
# Merge slope results with population genetics data on clonal rate and replicate
merged_slope_parsed_data <- merge(slope_results, parsed_data %>% select(clonalrate, Replicate,Mean_FIS_Tot, Var_FIS_Tot, R_Tot, Pareto_beta_Tot),
                     by.x = c("clonalrate", "Replicate"),
                     by.y = c("clonalrate", "Replicate"))


# Plot Slope vs Var_FIS for each clonal rate
ggplot(merged_slope_parsed_data, aes(x = Var_FIS_Tot, y = exp(Slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    title = "Relationship between Slope and Var(FIS) across Clone Rates",
    x = "Var(FIS)",
    y = "Slope",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Mean_FIS for each clonal rate
ggplot(merged_slope_parsed_data, aes(x = Mean_FIS_Tot, y = exp(Slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    title = "Relationship between Slope and Mean_FIS across Clone Rates",
    x = "Mean_FIS",
    y = "Slope",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )



# Plot Slope vs R for each clonal rate
ggplot(merged_slope_parsed_data, aes(x = R_Tot, y = exp(Slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    title = "Relationship between Slope and R across Clone Rates",
    x = "R",
    y = "Slope",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )

# Plot Slope vs Pareto_beta for each clonal rate
ggplot(merged_slope_parsed_data, aes(x = Pareto_beta_Tot, y = exp(Slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
    title = "Relationship between Slope and Pareto_beta across Clone Rates",
    x = "Pareto_beta",
    y = "Slope",
    color = "Clonal Rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 10)
  )





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
```{r}
# Merge slope results with population genetics data on clonal rate and replicate
merged_slope_age_count <- merge(slope_results, parsed_data %>% select(clonalrate, Replicate,Mean_FIS, Var_FIS,Mean_r_bar_D, SD_r_bar_D, R, Pareto_beta),
                     by.x = c("Clonal_Rate", "Replicate"),
                     by.y = c("clonalrate", "Replicate"))


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


# Add mean age of clones on same plot diff color or another plot
```



```{r}
# Assuming parsed_data is your data frame containing Age, Count, and clonalrate
# Filter for a specific clonal rate and create log_Count
parsed_data_new <- parsed_data %>%
  mutate(log_Count = log(Count)) %>%  # Create log_Count variable
  filter(clonalrate == "0.1") %>%      # Adjust clonal rate as needed
  arrange(Age)                          # Sort by Age

window_size <- 5  # Size of the sliding window
step_size <- 1    # Step size for the sliding window


# Initialize a data frame to store results
sliding_window_results <- data.frame(Age = numeric(), Mean = numeric(), Slope = numeric(), Volatility = numeric())

# Sliding window loop
for (i in seq(1, nrow(parsed_data_new) - window_size + 1, by = step_size)) {
  window_data <- parsed_data_new[i:(i + window_size - 1), ]  # Extract the current window
  
  # Calculate mean of the Count
  mean_age <- mean(window_data$Age)
  
  # Fit a linear model if there are enough data points
  if (nrow(window_data) > 1) {  # Ensure there are at least 2 points
    lm_fit <- lm(log_Count ~ Age, data = window_data)  # Fit linear model
    slope <- coef(lm_fit)["Age"]  # Extract the slope
  } else {
    slope <- NA  # If not enough data, assign NA
  }
  
  # Calculate volatility (standard deviation)
  volatility <- sd(window_data$Count)  # Standard deviation of Count
  
  # Store the results
  sliding_window_results <- rbind(sliding_window_results, data.frame(Age = mean_age, Mean = mean(window_data$Count), Slope = slope, Volatility = volatility))
}

# Plot Mean and Volatility over Sliding Window
ggplot(sliding_window_results, aes(x = Age)) +
  geom_line(aes(y = Mean, color = "Mean")) +
  geom_line(aes(y = Volatility, color = "Volatility")) +
  labs(title = "Mean and Volatility Over Sliding Windows",
       x = "Mean Age (Sliding Window)",
       y = "Values") +
  scale_color_manual(values = c("Mean" = "blue", "Volatility" = "red")) +
  theme_minimal()


```

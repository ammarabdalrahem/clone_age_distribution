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

required_packages <- c("knitr","rmarkdown","ggridges","data.table","ggplot2","gganimate","dplyr","tidyr","cowplot", "reshape2","ggcorrplot","corrplot")
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
data <- data_dif_mutation %>% filter(Nb_alleles_tot == "4") %>% filter(Mutation_Rate == 0.000001) %>% filter (Generation %% 10 == 0) 

#data_dif_mutation_m4_N4$Mean_FIS_Tot <- as.numeric(data_dif_mutation_m4_N4$Mean_FIS_Tot)
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
  width = 16,  # Double-column width (cm)
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
# filter the data to minimize the clonal rates for clear visualization
data_min_c <- data %>% filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999"))

# group by generation and calculate the median per each indices
data_grouped_g_c <- data_min_c %>%
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
plot_Mean_FIS_2 <-  ggplot(data_grouped_g_c, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Mean("*italic(F)[IS]*")")) +
  labs(
    x = "Generation",
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot variance of FIS by clonal rate and generation
plot_Var_FIS_2 <- ggplot(data_grouped_g_c, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Var_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Var("*italic(F)[IS]*")")) +
  labs(
    x = "Generation",
    y = "Clonal Rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


 # Plot mean He by clonal rate and generation
plot_Mean_He_2 <- ggplot(data_grouped_g_c, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_He_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[e])) +
  labs(
    y = "Clonal Rate",
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
plot_Mean_Ho_2 <- ggplot(data_grouped_g_c, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_Ho_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[o])) +
  labs(
    y = "Clonal Rate",
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
plot_R_2 <-  ggplot(data_grouped_g_c, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = R_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "R") +
  labs(
    y = "Clonal Rate",
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


# Plot Pareto beta by clonal rate and generation
plot_Pareto_beta_2 <- ggplot(data_grouped_g_c, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Pareto_beta_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(β) * " Pareto")) +
  labs(
    y = "Clonal Rate",
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



# Combine the plots 
combined_plot_2 <- plot_grid(plot_Mean_He_2, plot_Mean_Ho_2, plot_R_2 , plot_Pareto_beta_2,  plot_Var_FIS_2, plot_Mean_FIS_2, label_size = 12, labels = c("a", "b", "c", "d", "e","f"), ncol = 2)

# Save plot
ggsave(
  filename = "combined_pop_indices_plot_2.png",  
  plot = combined_plot_2,           
  width = 16,  
  height = 15,  
  units = "cm",                   
  dpi = 1200                       
)

combined_plot_2


```



## R^2 per each replicate

```{r, echo=FALSE}
# data for generation 100
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


# Plot 
R_2_plot<- ggplot(r2_results, aes(x = clonalrate, y = R2)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +  labs(
    x = "Clonal rate",
    y = expression(R^2),
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
combined_r2_plot <- plot_grid(R_2_plot, slope_plot, labels = c("A", "B"), ncol = 1)
# Save plot
ggsave(
  filename = "combined_r2_plot.png",  
  plot = combined_r2_plot,           
  width = 25,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)

combined_r2_plot
```


# Summary 


```{r, echo=FALSE}
# Create a summary table for the population genetic indices per generation 100

# Filter and prepare data for generation 6000
gen6000_data <- parsed_data %>% filter(Generation == 6000)

# Summarise genetic indices
summary_indices <- data_grouped_g_c %>%
  filter(Generation == 6000) %>%
  group_by(clonalrate) %>%
  summarise(
    median_R_Tot = median(R_Tot, na.rm = TRUE),
    median_Var_FIS = median(Var_FIS_Tot, na.rm = TRUE),
    median_FIS = median(Mean_FIS_Tot, na.rm = TRUE),
    median_He = median(Mean_He_Tot, na.rm = TRUE),
    median_Ho = median(Mean_Ho_Tot, na.rm = TRUE),
    median_Pareto_beta = median(Pareto_beta_Tot, na.rm = TRUE),
    .groups = 'drop'
  )

# Summarise R2 and slope
r2_summary <- r2_results %>%
  group_by(clonalrate) %>%
  summarise(
    mean_R2 = median(R2, na.rm = TRUE),
    mean_slope = mean(slope, na.rm = TRUE),
    .groups = 'drop'
  )

# Summarise age statistics from parsed_data
age_summary <- gen6000_data %>%
  group_by(clonalrate) %>%
  summarise(
    max_age = max(Age, na.rm = TRUE),
    median_age = median(rep(Age, Count), na.rm = TRUE),
    var_age = var(Age, na.rm = TRUE),
    .groups = 'drop'
  )


# Ensure clonalrate is numeric in all summary tables
summary_indices$clonalrate <- as.numeric(as.character(summary_indices$clonalrate))
r2_summary$clonalrate <- as.numeric(as.character(r2_summary$clonalrate))
age_summary$clonalrate <- as.numeric(as.character(age_summary$clonalrate))


# Find common clonalrate values across all summary tables
common_clonalrates <- unique(r2_summary$clonalrate)
# Filter each table to keep only common clonalrates
summary_indices_filtered <- summary_indices %>%
  filter(clonalrate %in% common_clonalrates)
age_summary_filtered <- age_summary %>%
  filter(clonalrate %in% common_clonalrates)



# Now join only the filtered tables
summary_table <- summary_indices_filtered %>%
  left_join(r2_summary, by = "clonalrate") %>%
  left_join(age_summary_filtered, by = "clonalrate")


```

## spearman correlation

```{r, echo=FALSE}

# Prepare numeric-only data
cor_data <- summary_table %>%
  dplyr::select(-clonalrate)

# Compute Spearman correlation
cor_matrix <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")

# Rename the columns and rows
colnames(cor_matrix) <- c("R", "Var_FIS", "Mean_FIS", "He", "Ho", "Pareto_beta", "R^2", "Slope",
                           "Max_age", "Median_age","Var_age")
rownames(cor_matrix) <- c("R", "Var_FIS", "Mean_FIS", "He", "Ho", "Pareto_beta", "R^2", "Slope",
                           "Max_age",  "Median_age","Var_age")

file_path= "Spearman_correlation_plot_6000_matrix.png"
png(height=20, width=20, file=file_path, units = "cm", res=1200,bg="white")

Spearman_correlation_plot_6000 <-  corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         tl.cex= 1,
         addCoef.col = "black",
         tl.col = "black",
         mar = c(0,0,1,0),
         bg = "white")

dev.off()




```

```{r, echo=FALSE}

#plotting the median R_tot and slope of age distribution (each dot represent clonerate)
plot_R_slope <- ggplot(summary_table, aes(x = mean_slope, y = median_R_Tot)) +
  geom_point(aes(color =as.factor(clonalrate)), size = 3) +  # Points with transparency
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Optional trend line with confidence interval
  labs(
    x = "Slope of age distribution",
    y = "Median R"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",  
    legend.direction = "horizontal"  
  ) +
  guides(color = guide_legend(nrow = 1)) +  # Arrange legend in one row
  scale_color_viridis_d(option = "C", end = 0.9)  # Use perceptually-uniform colors

plot_R_slope


# plot R_tot with max age 
plot_R_max_age <- ggplot(summary_table, aes(x = max_age, y = median_R_Tot)) +
  geom_point(aes(color =as.factor(clonalrate)), size = 3) +  # Points with transparency
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Optional trend line with confidence interval
  labs(
    x = "Max age",
    y = "Median R"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",  
    legend.direction = "horizontal"  
  ) +
  guides(color = guide_legend(nrow = 1)) +  # Arrange legend in one row
  scale_color_viridis_d(option = "C", end = 0.9)  # Use perceptually-uniform colors
plot_R_max_age


```













## Slope vs Var_FIS, Mean_FIS, R_Tot and Pareto_beta



```{r}
# Merge slope results with population genetics data on clonal rate and replicate
merged_slope_parsed_data <- merge(r2_results, parsed_data %>% select(clonalrate, Replicate,Mean_FIS_Tot, Var_FIS_Tot, R_Tot, Pareto_beta_Tot),
                     by.x = c("clonalrate", "Replicate"),
                     by.y = c("clonalrate", "Replicate"))


# Plot Slope vs Var_FIS for each clonal rate
ggplot(merged_slope_parsed_data, aes(x = Var_FIS_Tot, y = exp(slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
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
ggplot(merged_slope_parsed_data, aes(x = Mean_FIS_Tot, y = exp(slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
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
ggplot(merged_slope_parsed_data, aes(x = R_Tot, y = exp(slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
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
ggplot(merged_slope_parsed_data, aes(x = Pareto_beta_Tot, y = exp(slope), color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line for each clonal rate
  labs(
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



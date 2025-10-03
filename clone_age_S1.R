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

```{r}
#Code to install packages if necessary, and read them with library function

required_packages <- c("knitr","rmarkdown","ggridges","data.table","ggplot2","gganimate","dplyr","tidyr","cowplot", "reshape2","ggcorrplot","corrplot","viridis","ggforce","ggbeeswarm")
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



```{r}
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

# Filter the data, alleles 4 and mutation rate 1e-03 and generation every 10
data <- data_dif_mutation %>% filter(Nb_alleles_tot == "4") %>% filter(Mutation_Rate == 0.001) %>% filter (Generation %% 10 == 0) 

head (data)

```



## The relation between population genetic indices and clonerate 
# As already shown in Stoeckel et al. (2021), this relationship, with an initial population size of 10^4 Figure S1a.

```{r}
# Filter the data for generation 6000
data_6000 <- data %>% filter (Generation == 6000)


# Create the boxplot for R
plot_R <- ggplot(data_6000, aes(x = as.factor(clonalrate), y = R_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = "Clonality rate",
    y = expression(italic(R))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10),size = 12),
    axis.title.y = element_text(margin = margin(r = 10), size = 12)
  )

# Create the boxplot for Pareto_beta
plot_Pareto_beta <- ggplot(data_6000, aes(x = as.factor(clonalrate), y = Pareto_beta_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = "Clonality rate",
    y = expression(italic(β) * " Pareto")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10),size = 12),
    axis.title.y = element_text(margin = margin(r = 10), size = 12)
  )


#Creat the boxplot for Var_FIS
plot_Var_FIS <- ggplot(data_6000, aes(x = as.factor(clonalrate), y = Var_FIS_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression("Var("*italic(F)[IS]*")")
) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10),size = 12),
    axis.title.y = element_text(margin = margin(r = 10), size = 12)
  )


# Creat the boxplot for Mean_FIS
plot_Mean_FIS <- ggplot(data_6000, aes(x = as.factor(clonalrate), y = Mean_FIS_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression("Mean("*italic(F)[IS]*")")
) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10),size = 12),
    axis.title.y = element_text(margin = margin(r = 10), size = 12)
  )

# Creat the boxplot for Mean_He
plot_Mean_He <- ggplot(data_6000, aes(x = as.factor(clonalrate), y = Mean_He_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression(italic(H)[E])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10),size = 12),
    axis.title.y = element_text(margin = margin(r = 10), size = 12)
  )

#Creat the boxplot for Mean_Ho
plot_Mean_Ho <- ggplot(data_6000, aes(x = as.factor(clonalrate), y = Mean_Ho_Tot)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +
  labs(
    x = NULL,
    y = expression(italic(H)[O])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10),size = 12),
    axis.title.y = element_text(margin = margin(r = 10), size = 12)
  )




# Combine the plots 
combined_pop_indices_plot <- plot_grid(  plot_Mean_FIS, plot_Var_FIS, 
                                        plot_Mean_He, plot_Mean_Ho,
                                        plot_R, plot_Pareto_beta,
                                        label_size = 12, 
                                        labels = c("a", "b", "c", "d", "e","f"),
                                        ncol = 2)

# Save plot
ggsave(
  filename = "combined_pop_indices_plot.png",  
  plot = combined_pop_indices_plot,           
  width = 18,  # Double-column width (cm)
  height = 15,   # Adjust height to maintain aspect ratio
  units = "cm",                   
  dpi = 1200                       
)

combined_pop_indices_plot



```



# prepare data for joyplot

```{r, echo=FALSE}
# Create a list to store clone ages for each group
clone_ages_by_group <- list()
# max generation  data 9000
data_9000 <- data %>% filter(Generation == 9000) # Filter the data for generation 9000
# Process (loop) the data by row 
for (i in seq_len(nrow(data_9000))) {
  group_value <- as.character(data_9000$clonalrate[i]) # Convert clonerate to character (ex. "0.7")
  age_distribution <- gsub("[{}]", "", data_9000$List_distribution_gen_clonal_genotypes[i]) # Remove curly braces (ex. {0: 300, 1: 213} → "0: 300, 1: 213")
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



# bound function

bounded_kde <- function(x, 
                        lower = 0, 
                        upper = 6000, 
                        method = c("reflection", "logit"),
                        bw = "nrd0", 
                        n = 512) {
  
  method <- match.arg(method)
  
  if (method == "reflection") {
    # --- Méthode par réflexion ---
    x_reflect <- c(x, 2*lower - x, 2*upper - x)
    d <- density(x_reflect, from = lower, to = upper, bw = bw, n = n)
    
    # Normalisation
    dx <- d$x[2] - d$x[1]
    d$y <- d$y / sum(d$y * dx)
    return(list(x = d$x, y = d$y, method = "reflection"))
    
  } else if (method == "logit") {
    # --- Méthode par logit ---
    if (!requireNamespace("ks", quietly = TRUE)) {
      stop("Le package 'ks' est requis pour la méthode logit. Installez-le avec install.packages('ks').")
    }
    
    logit_transform <- function(x) log((x - lower)/(upper - x))
    inv_logit <- function(z) (upper*exp(z) + lower)/(1 + exp(z))
    
    # Transformation
    z <- logit_transform(x)
    fit <- ks::kde(z)
    
    # Grille sur l'échelle originale
    grid_x <- seq(lower, upper, length.out = n)
    grid_z <- logit_transform(grid_x)
    
    # Dérivée de la transformation pour le changement de variable
    dzdx <- 1/(grid_x - lower) + 1/(upper - grid_x)
    
    # Densité estimée
    dens <- ks::predict(fit, x = grid_z) * dzdx
    
    return(list(x = grid_x, y = dens, method = "logit"))
  }
}





```

##joyplot

```{r, echo=FALSE}
# Create a data frame for plotting
df <- data.frame(Clone_Age = all_ages, Group = factor(all_groups, levels = rev(sort(unique(all_groups)))))


# Define bounds
lower_bound <- 0
upper_bound <- 9000


# Define bounds (min and max of Clone_Age)
lower_bound <- 0
upper_bound <- max(df$Clone_Age, na.rm = TRUE)  # e.g. 9000

# Compute bounded densities for each Clonality rate group
densities <- df %>%
  group_by(Group) %>%
  do({
    ages <- .$Clone_Age
    
    # Case 1: all values identical (e.g., all 0)
    if (length(unique(ages)) == 1) {
      tibble(Clone_Age = unique(ages),
             Density   = 1,   # spike at 0
             Group     = unique(.$Group),
             special   = TRUE)
      
    # Case 2: variable ages → use bounded KDE
    } else {
      res <- bounded_kde(ages,
                         lower = lower_bound,
                         upper = upper_bound,
                         method = "reflection",
                         bw = "nrd0",
                         n = 512)
      tibble(Clone_Age = res$x,
             Density   = res$y / max(res$y),   # normalize per group
             Group     = unique(.$Group),
             special   = FALSE)
    }
  }) %>%
  ungroup()



Joyplot_plot <- ggplot(densities, aes(x = Clone_Age, y = Group,
                                      height = Density, fill = Group)) +
  geom_ridgeline(scale = 3, alpha = 0.7, color = NA) +
  labs(
    x = "Clone ages (Generations)",
    y = "Clonality rate"
  ) +
  theme_ridges() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14, hjust = 0.5),
    axis.title.y = element_text(size = 14, angle = 90, hjust = 0.5),
    legend.position = "none"
  ) + # zoom in
  facet_zoom(x = Clone_Age >= 0 & Clone_Age <= 2000, zoom.size = 1)


# summary stats
df_summary <- df %>%
    group_by(Group) %>%
    summarise(
        mean_val   = mean(Clone_Age, na.rm = TRUE),
        median_val = median(Clone_Age, na.rm = TRUE),
        .groups = "drop"
    )


# jittered raw data (grey)
jittered_plot <- ggplot(df, aes(x = Clone_Age, y = Group)) +
  geom_quasirandom(width = 0.3, alpha = 0.4, size = 0.4, aes(color = Group)) +
  geom_point(data = df_summary,
             aes(x = mean_val, y = Group, shape = "Mean"),
             color = "red", size = 2, inherit.aes = FALSE) +
  geom_point(data = df_summary,
             aes(x = median_val, y = Group, shape = "Median"),
             color = "blue", size = 2, inherit.aes = FALSE) +
  scale_shape_manual(values = c(Mean = 16, Median = 17), breaks = c("Mean", "Median"), name = "") +
  guides(color = "none",  # This removes the clonality rate legend #color = guide_legend(reverse = TRUE)
         shape = guide_legend(override.aes = list(color = c("red", "blue"), size = 3))) +
  labs(x = "Clone ages (Generations)", y = "Clonality rate", shape = "") +
  theme_classic() +
  facet_zoom(x = Clone_Age >= 0 & Clone_Age <= 150, zoom.size = 1) +
  theme(legend.position = "right")

#save new plot
ggsave(
  filename = "jittered_ages_by_clonerate_plot_points.png",  
  plot = jittered_plot,           
  width = 15,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)

# Joyplot 
# Save plot
ggsave(
  filename = "Joyplot_clone_ages_by_clonerateplot.png",  
  plot = Joyplot_plot,           
  width = 15,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)

Joyplot_plot

```



## Scatter plot for all replicates

```{r, echo=FALSE}
# choose the max generation 9000
# Filter the data for generation 9000 and specific Clonality rates
data_9000 <- data %>% filter(Generation == 9000) %>%
  filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", 
                           "0.95", "0.98", "0.99", "0.999", "0.9999"))

parsed_data <- data_9000 %>%
  rowwise() %>%
  mutate(parsed_list = strsplit(gsub("[{}]", "", List_distribution_gen_clonal_genotypes), ", ")) %>%
  unnest(parsed_list) %>%
  separate(parsed_list, into = c("Age", "Count"), sep = ":") %>%
  mutate(
    Age = as.numeric(Age),
    Count = as.numeric(Count)
  ) %>%
  ungroup()

# filter parsed data
#parsed_data_new <- parsed_data %>% filter(clonalrate== 0.1)


# Plot the data 
scatter_plot_all_ages <- ggplot(parsed_data, aes(x = Age, y = Count, color = as.factor(clonalrate))) +
  geom_point(alpha = 0.7, size = 2) +  # Increase point size slightly
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = "Clone age (log scale)",
    y = "Number of individuals (log scale)",
    color = "Clonality rate"
  ) +
  theme_minimal(base_size = 14) +
  #scale_color_viridis_d(option = "viridis") +
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

scatter_plot_all_ages


loess_plot_all_ages <- ggplot(parsed_data,
       aes(x = Age, y = Count,
           color = as.factor(clonalrate), fill = as.factor(clonalrate))) +
    geom_point(alpha = 0.3, size = 0.6, show.legend = FALSE) +
    geom_smooth(method = "loess", 
                span = 0.3,  # Adjust smoothness
                se = TRUE, linewidth = 1.3) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Clone age (log scale)", y = "Number of individuals (log scale)", 
         color = "Clonality rate", fill = "Clonality rate") +
    theme_minimal(base_size = 14) +
  #scale_color_viridis_d(option = "viridis") +
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


glm_plot_all_ages <- ggplot(parsed_data,
       aes(x = Age + 1, y = Count,
           color = as.factor(clonalrate), fill = as.factor(clonalrate))) +
    geom_point(alpha = 0.2, size = 0.6, show.legend = FALSE) +
    geom_smooth(method = "lm",
                formula = y ~ x,  
                se = TRUE, linewidth = 1) +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = "Clone age (log scale)", y = "Number of individuals (log scale)", 
         color = "Clonality rate", fill = "Clonality rate") +
    theme_minimal(base_size = 14) +
  #scale_color_viridis_d(option = "viridis") +
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

# Save plot
ggsave(
  filename = "loess_ages_by_clonerate_plot.png",  
  plot = loess_plot_all_ages,           
  width = 25,                  
  height = 15,    
  units = "cm",                   
  dpi = 1200                       
)

# Save plot
ggsave(
  filename = "glm_ages_by_clonerate_plot.png",  
  plot = glm_plot_all_ages,           
  width = 25,                  
  height = 15,    
  units = "cm",                   
  dpi = 1200                       
)

# Animated HTML

Animated_scatter_plot <- ggplot(parsed_data, aes(x = Age, y = Count, color = as.factor(clonalrate))) +
  geom_point(alpha = 0.7) +  # transparency
  scale_y_log10() +  # Set y-axis to log scale
  scale_x_log10() +  # Set x-axis to log scale
  labs(
    # Remove the original title so it doesn't show on top
    x = "Clone age (log scale)",
    y = "Number of individuals (log scale)",
    color = "Clonality rate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 1)) +  # Arrange legend in one row
  transition_states(clonalrate, transition_length = 2, state_length = 1) +
  # Add dynamic text label inside the plot at custom coordinates:
  geom_text(
    aes(x = 0.1, y = max(Count), label = paste("Clonality rate:", clonalrate)), # Adjust x and y for position
    color = "black",
    size = 5,
    hjust = 0,
    vjust = 1,
    inherit.aes = FALSE
  )

# Animate the plot
animate(Animated_scatter_plot, nframes = 100, fps = 10, width = 1500, height = 800, res = 150, renderer = gifski_renderer())
#anim <- animate(Animated_scatter_plot, nframes = 100, fps = 10, width = 1500, height = 800, res = 150, renderer = gifski_renderer())

# Save animation
#anim_save("clone_age_animation_wide_high_res.gif", animation = anim)



```



## The relation between population genetic indices, generations and clonerate
#heatmap

```{r, echo=FALSE}
# filter the data to minimize the Clonality rates for clear visualization
data_min_c_10_3 <- data %>% filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999", "0.9999"))


# group by generation and calculate the median per each indices
data_grouped_g_c_10_3 <- data_min_c_10_3 %>%
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


# Plot mean FIS by Clonality rate and generation
plot_Mean_FIS_10_3 <-  ggplot(data_grouped_g_c_10_3, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Mean("*italic(F)[IS]*")")) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot variance of FIS by Clonality rate and generation
plot_Var_FIS_10_3 <- ggplot(data_grouped_g_c_10_3, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Var_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Var("*italic(F)[IS]*")")) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


 # Plot mean He by Clonality rate and generation
plot_Mean_He_10_3 <- ggplot(data_grouped_g_c_10_3, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_He_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[e])) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot mean Ho by Clonality rate and generation
plot_Mean_Ho_10_3 <- ggplot(data_grouped_g_c_10_3, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_Ho_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[o])) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Plot R by Clonality rate and generation
plot_R_10_3 <-  ggplot(data_grouped_g_c_10_3, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = R_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(R))) +
  labs(
    x = "Generation",
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(margin = margin(r = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# Plot Pareto beta by Clonality rate and generation
plot_Pareto_beta_10_3 <- ggplot(data_grouped_g_c_10_3, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Pareto_beta_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(β) * " Pareto")) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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
data_mutation_10_4 <- data_dif_mutation %>% filter(Nb_alleles_tot == "4") %>% filter(Mutation_Rate == 0.000001) %>% filter (Generation %% 10 == 0) 



```


## The relation between population genetic indices, generations and clonerate for mutation rate 10^-4

```{r, echo=FALSE}
# filter the data to minimize the Clonality rates for clear visualization
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


# Plot mean FIS by Clonality rate and generation
plot_Mean_FIS_10_4 <-  ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Mean("*italic(F)[IS]*")")) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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


# Plot variance of FIS by Clonality rate and generation
plot_Var_FIS_10_4 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Var_FIS_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression("Var("*italic(F)[IS]*")")) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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


 # Plot mean He by Clonality rate and generation
plot_Mean_He_10_4 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_He_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[e])) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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


# Plot mean Ho by Clonality rate and generation
plot_Mean_Ho_10_4  <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Mean_Ho_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(H)[o])) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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

# Plot R by Clonality rate and generation
plot_R_10_4 <-  ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = R_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(R))) +
  labs(
    x = "Generation",
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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


# Plot Pareto beta by Clonality rate and generation
plot_Pareto_beta_10_4 <- ggplot(data_grouped_g_c_10_4, aes(x = as.factor(Generation), y = as.factor(clonalrate), fill = Pareto_beta_Tot)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = expression(italic(β) * " Pareto")) +
  labs(
    y = "Clonality rate",
  ) +
  theme_minimal(base_size = 14) +
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
combined_pop_indices_generations_heatmap <- plot_grid (plot_Mean_FIS_10_3, plot_Mean_FIS_10_4, plot_Var_FIS_10_3, 
                                                       plot_Var_FIS_10_4, plot_Mean_He_10_3,plot_Mean_He_10_4, plot_R_10_3, plot_R_10_4,
                                                        label_size = 12, labels = c("a", "b", "c", "d","e","f","g","h"), ncol = 2)


combined_pop_indices_generations_heatmap


# Save plot
ggsave(
  filename = "combined_pop_indices_plot_heatmap.png",
  plot = combined_pop_indices_generations_heatmap,
  width = 30,  
  height = 20,  
  units = "cm",
  dpi = 1200
)



```


## R^2 per each replicate
# we calculate R^2 for each replicate and Clonality rate at generation 6000
```{r, echo=FALSE}
# data for generation 6000

parsed_data_g <- data %>% 
  filter(Generation == 6000) %>%
  filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999", "0.9999")) %>%
  rowwise() %>% # Process each row individually
  mutate(parsed_list = strsplit(gsub("[{}]", "", List_distribution_gen_clonal_genotypes), ", ")) %>% # Remove curly braces and split by comma
  unnest(parsed_list) %>% # Unnest the list into separate rows. example "0: 300, 1: 213" becomes "0: 300" and "1: 213"
  separate(parsed_list, into = c("Age", "Count"), sep = ":") %>% # Separate age and count into two columns. example "0: 300" becomes "0" and "300"
  mutate(
    Age = as.numeric(Age),
    Count = as.numeric(Count)
  ) %>%
  ungroup()



# Initialize a data frame to store the results
r2_results <- data.frame(clonalrate = character(), Replicate = integer(), R2 = numeric(),slope=numeric(), stringsAsFactors = FALSE)

# Get unique Clonality rates
unique_clonal_rates <- unique(parsed_data_g$clonalrate)

# Loop over each Clonality rate to calculate R^2 and slopes
for (clonal_rate in unique_clonal_rates) {
  # Filter data for the current Clonality rate
  parsed_data_clonal <- parsed_data_g %>%
    filter(clonalrate == clonal_rate) %>%
    mutate(log_Count = log(Count))  # Log-transform Count

  # Loop over unique replicates for the current Clonality rate
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
# Calculate the maximum, median, and variance of clone ages for each Clonality rate and replicate

age_stats <- parsed_data_g %>%
  group_by(clonalrate, Replicate) %>%
  summarise(
    Max_Age = max(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    Var_Age = var(Age, na.rm = TRUE),
    .groups = "drop"
  )

head (age_stats)

```



## clonal age summary as function of genetic indices 6000


```{r}
# Filter data for generation 6000
data_6000 <- data %>%
  filter(Generation == 6000) %>%
  filter(clonalrate %in% c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999", "0.9999"))




# Ensure clonalrate is of the same type (character) across datasets for merging
r2_results$clonalrate   <- as.character(r2_results$clonalrate)
data_6000$clonalrate    <- as.character(data_6000$clonalrate)
age_stats$clonalrate    <- as.character(age_stats$clonalrate)

# Merge slope/R² results with raw data by Clonality rate and replicate
merged_1 <- merge(
  r2_results,
  data_6000,
  by = c("clonalrate", "Replicate")
)

# Merge with age summary statistics
merged_slope_parsed_df <- merge(
  merged_1,
  age_stats,
  by = c("clonalrate", "Replicate")
)

# Filter to keep only selected Clonality rates for focused analysis/visualization
selected_clonalrates <- c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.98", "0.99", "0.999","0.9999")



merged_slope_parsed_data <- merged_slope_parsed_df %>%
  filter(clonalrate %in% selected_clonalrates)

# Remove irrelevant columns for clarity in final dataset
merged_slope_parsed_data <- dplyr::select(merged_slope_parsed_data,
   -List_distribution_gen_clonal_genotypes,
   -Generation,
   -Mutation_Rate,
   -Nb_alleles_tot,
   -Number_alleles_Tot,
   -Number_fixed_loci_Tot,
   -Number_Genotypes_Tot
)
```

# Plot R^2 and slope for each Clonality rate
```{r, echo=FALSE}

# Plot 
R_2_plot<- ggplot(merged_slope_parsed_data, aes(x = clonalrate, y = R2)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +  labs(
    y = expression(R^2),
    x = "Clonality rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


#R_2_plot

# Plot slope
slope_plot <- ggplot(merged_slope_parsed_data, aes(x = clonalrate, y = slope)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) +  labs(
    y = "Slope",
    x = "Clonality rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


var_age_plot <- ggplot(merged_slope_parsed_data, aes(x = clonalrate, y = Var_Age)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) + labs(
    y = "Var(Age)",
    x = "Clonality rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )


max_age_plot <- ggplot(merged_slope_parsed_data, aes(x = clonalrate, y = Max_Age)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) + labs(
    y = "Max(Age)",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )

median_age_plot <- ggplot(merged_slope_parsed_data, aes(x = clonalrate, y = Median_Age)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) + labs(
    y = "Median(Age)",
    x = "Clonality rate",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Summarize data to count distinct ages per Clonality rate and replicate
age_counts <- parsed_data %>%
  group_by(clonalrate, Replicate) %>%
  summarize(Age_Count = n_distinct(Age), groups = "drop")  # Count unique ages per replicate within each Clonality rate

# Convert Clonality rate to a factor, so it appears in order on the x-axis
age_counts$clonalrate <- factor(age_counts$clonalrate, levels = sort(unique(age_counts$clonalrate)))



age_class_plot <- ggplot(age_counts, aes(x = clonalrate, y = Age_Count)) +
  geom_boxplot(alpha = 0.8,
    outlier.color = "grey50",
    outlier.shape = 16,
    fill = "#4B9CD3", 
    color = "black"
  ) + labs(
    y = "Number of age classes",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# combine the plots
combined_r2_slope_age_plot <- plot_grid(
  age_class_plot, max_age_plot, R_2_plot, slope_plot,
  labels = c("a","b", "c", "d"),
  label_size = 12,
  ncol = 2,
  align = 'v'
  )


# Save plot
ggsave(
  filename = "combined_r2_slope_age_plot.png",  
  plot = combined_r2_slope_age_plot,           
  width = 18,                  
  height = 15,    
  units = "cm",                   
  dpi = 1200                       
)


# combine var and median age plots
combined_var_median_age_plot <- plot_grid(
  var_age_plot, median_age_plot,
  labels = c("a", "b"),
  label_size = 12,
  ncol = 2,
  align = 'v'
)

# save the plot
ggsave(
  filename = "combined_var_median_age_plot.png",  
  plot = combined_var_median_age_plot,           
  width = 18,                  
  height = 8,    
  units = "cm",                   
  dpi = 1200                       
)

combined_r2_slope_age_plot
```


## spearman correlation using merged_slope_parsed_data

```{r, echo=FALSE}
# Prepare numeric-only data (exclude categorical variables)
cor_data <- merged_slope_parsed_data %>% 
  dplyr::select(-clonalrate, -Replicate)

# Compute Spearman correlation matrix
cor_matrix <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")

# Rename columns and rows for better readability
colnames(cor_matrix) <- c("R²", "Slope", "Mean (He)", "Mean (Ho)", "Mean (FIS)", "Var (FIS)",  "R", "Pareto β", "Max (Age)", "Median (Age)", "Var (Age)")
                          
                          
rownames(cor_matrix) <- c("R²", "Slope", "Mean (He)", "Mean (Ho)", "Mean (FIS)", "Var (FIS)",  "R", "Pareto β", "Max (Age)", "Median (Age)", "Var (Age)")

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
                                     col = viridis(200, option = "D"),  # Colorblind-friendly palette
                                     tl.srt = 45,                     # Rotate labels for better readability
                                     cl.cex = 0.8,                    # Color legend text size
                                     cl.ratio = 0.2)                 # Color legend width



# Close the PNG device
dev.off()



```




# Plotting the relationship between slope and various population genetic indices

```{r, echo=FALSE}


common_theme <- theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(face = "italic", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "italic", size = 13, margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )


create_slope_plot <- function(data, y_var, y_label, show_x_label = FALSE) {
  p <- ggplot(data, aes(x = slope, y = !!sym(y_var), color = as.factor(clonalrate))) +
    geom_point(size = 3, alpha = 0.7) +                    
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
    labs(
      y = y_label,
      color = "Clonality rate"
    ) +
    common_theme
  
  if (show_x_label) {
    p <- p + labs(x = "Slope")
  } else {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  return(p)
}


# Create plots, showing x-axis label only on last two
slope_var_fis <- create_slope_plot(merged_slope_parsed_data, "Var_FIS_Tot",  expression("Var("*italic(F)[IS]*")"), show_x_label = FALSE)
slope_mean_fis <- create_slope_plot(merged_slope_parsed_data, "Mean_FIS_Tot",  expression("Mean("*italic(F)[IS]*")"), show_x_label = FALSE)
slope_Ho      <- create_slope_plot(merged_slope_parsed_data, "Mean_Ho_Tot", expression(italic("Mean(H"[O]*")")), show_x_label = FALSE)
slope_He      <- create_slope_plot(merged_slope_parsed_data, "Mean_He_Tot", expression(italic("Mean(H"[E]*")")), show_x_label = FALSE)
slope_R       <- create_slope_plot(merged_slope_parsed_data, "R_Tot", expression(italic("R")), show_x_label = TRUE)
slope_Pareto_beta <- create_slope_plot(merged_slope_parsed_data, "Pareto_beta_Tot", expression(italic("Pareto " * beta)), show_x_label = TRUE)


# Extract legend from the first plot
legend_shared <- get_legend(
  slope_var_fis +
    theme(
      legend.position = "top",
      legend.justification = "center"
    )
)




# Combine all plots with 2 columns, labels a-f
plot_pop_indices_slope <- plot_grid(
  slope_var_fis, slope_mean_fis,
  slope_Ho, slope_He,
  slope_R, slope_Pareto_beta,
  labels = c("a", "b", "c", "d", "e", "f"),
  ncol = 2,
  align = "v"
)

final_plot_pop_indices_slope <- plot_grid(
  legend_shared,             
  plot_pop_indices_slope,    
  ncol = 1,
  rel_heights = c(0.1, 1)    
)


final_plot_pop_indices_slope

ggsave(
  filename = "combined_slope_population_genetic_indices.png",  # Save the combined plot
  plot = final_plot_pop_indices_slope,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )
```




# Plotting the relationship between R^2 and various population genetic indices

```{r, echo=FALSE}
# Plot R^2 vs Var_FIS for each Clonality rate

create_slope_plot <- function(data, y_var, y_label, show_x_label = FALSE) {
  p <- ggplot(data, aes(x = R2, y = !!sym(y_var), color = as.factor(clonalrate))) +
    geom_point(size = 3, alpha = 0.7) +                    
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
    labs(
      y = y_label,
      color = "Clonality rate"
    ) +
    common_theme
  
  if (show_x_label) {
    p <- p + labs(x = expression(italic("R")^2))
  } else {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  return(p)
}


# Create plots, showing x-axis label only on last two
R2_var_fis <- create_slope_plot(merged_slope_parsed_data, "Var_FIS_Tot",  expression("Var("*italic(F)[IS]*")"), show_x_label = FALSE)
R2_mean_fis <- create_slope_plot(merged_slope_parsed_data, "Mean_FIS_Tot",  expression("Mean("*italic(F)[IS]*")"), show_x_label = FALSE)
R2_Ho      <- create_slope_plot(merged_slope_parsed_data, "Mean_Ho_Tot", expression(italic("Mean(H"[O]*")")), show_x_label = FALSE)
R2_He      <- create_slope_plot(merged_slope_parsed_data, "Mean_He_Tot", expression(italic("Mean(H"[E]*")")), show_x_label = FALSE)
R2_R       <- create_slope_plot(merged_slope_parsed_data, "R_Tot", expression(italic("R")), show_x_label = TRUE)
R2_Pareto_beta <- create_slope_plot(merged_slope_parsed_data, "Pareto_beta_Tot", expression(italic("Pareto " * beta)), show_x_label = TRUE)


# Combine all plots with 2 columns, labels a-f
plot_pop_indices_R2 <- plot_grid(
  R2_var_fis, R2_mean_fis,
  R2_Ho, R2_He,
  R2_R, R2_Pareto_beta,
  labels = c("a", "b", "c", "d", "e", "f"),
  ncol = 2,
  align = "v"
)

final_plot_pop_indices_R2 <- plot_grid(
  legend_shared,             
  plot_pop_indices_R2,    
  ncol = 1,
  rel_heights = c(0.1, 1)    
)



final_plot_pop_indices_R2


# save the combined plot
ggsave(
  filename = "combined_R2_population_genetic_indices.png",  # Save the combined plot
  plot = final_plot_pop_indices_R2,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )


```


# Plotting the relationship between Max_age and various population genetic indices
```{r, echo=FALSE}
# Create plots for Max_age vs various population genetic indices
create_max_age_plot <- function(data, y_var, y_label, show_x_label = FALSE) {
  p <- ggplot(data, aes(x = Max_Age, y = !!sym(y_var), color = as.factor(clonalrate))) +
    geom_point(size = 3, alpha = 0.7) +                    
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
    labs(
      y = y_label,
      color = "Clonality rate"
    ) +
    common_theme
  
  if (show_x_label) {
    p <- p + labs(x = "Max Age")
  } else {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  return(p)
}

# Create plots, showing x-axis label only on last two
Max_age_var_fis <- create_max_age_plot(merged_slope_parsed_data, "Var_FIS_Tot",  expression("Var("*italic(F)[IS]*")"), show_x_label = FALSE)
Max_age_mean_fis <- create_max_age_plot(merged_slope_parsed_data, "Mean_FIS_Tot",  expression("Mean("*italic(F)[IS]*")"), show_x_label = FALSE)
Max_age_Ho      <- create_max_age_plot(merged_slope_parsed_data, "Mean_Ho_Tot", expression(italic("Mean(H"[O]*")")), show_x_label = FALSE)
Max_age_He      <- create_max_age_plot(merged_slope_parsed_data, "Mean_He_Tot", expression(italic("Mean(H"[E]*")")), show_x_label = FALSE)
Max_age_R       <- create_max_age_plot(merged_slope_parsed_data, "R_Tot", expression(italic("R")), show_x_label = TRUE)
Max_age_Pareto_beta <- create_max_age_plot(merged_slope_parsed_data, "Pareto_beta_Tot", expression(italic("Pareto " * beta)), show_x_label = TRUE)

# Combine all plots with 2 columns, labels a-f
plot_pop_indices_Max_age <- plot_grid(
  Max_age_var_fis, Max_age_mean_fis,
  Max_age_Ho, Max_age_He,
  Max_age_R, Max_age_Pareto_beta,
  labels = c("a", "b", "c", "d", "e", "f"),
  ncol = 2,
  align = "v"
)

# Add legend at the top centered
legend_shared <- get_legend(
  Max_age_var_fis +
    theme(
      legend.position = "top",
      legend.justification = "center"
    )
)

final_plot_pop_indices_Max_age <- plot_grid(
  legend_shared,             
  plot_pop_indices_Max_age,    
  ncol = 1,
  rel_heights = c(0.1, 1)    
)

final_plot_pop_indices_Max_age

# Save the combined plot
ggsave(
  filename = "combined_Max_age_population_genetic_indices.png",  # Save the combined plot
  plot = final_plot_pop_indices_Max_age,  # The combined plot object
  width = 30,  # Width of the saved image
  height = 25,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )



```


## scatter plot for Max_age  vs Var_FIS and Max_age vs Mean_FIS for c= 0.99, 0.999 , 0.9999
```{r, echo=FALSE}
# Filter the merged_slope_parsed_data for Clonality rates 0.999 and 0.9999

merged_slope_parsed_df_extreme <- merged_slope_parsed_data %>%
  filter(clonalrate %in% c("0.99","0.999", "0.9999"))



# Plot max_age vs Mean_FIS for each Clonality rate
max_age_vs_mean_fis_plot <- ggplot(merged_slope_parsed_df_extreme, aes(x = Max_Age, y = Mean_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +                    
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
  labs(
    x = "Max Age",
    y = expression("Mean("*italic(F)[IS]*")"),
    color = "Clonality rate"
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "italic", size = 13, margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )


# Plot max_age vs Var_FIS for each Clonality rate
max_age_vs_var_fis_plot <- ggplot(merged_slope_parsed_df_extreme, aes(x = Max_Age, y = Var_FIS_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +                    
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
  labs(
    x = "Max Age",
    y = expression("Var("*italic(F)[IS]*")"),
    color = "Clonality rate"
  ) +
    theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(face = "italic", size = 13, margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks = element_line(size = 0.3),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )


# Plot max_age vs R for each Clonality rate
max_age_vs_R_plot <- ggplot(merged_slope_parsed_df_extreme, aes(x = Max_Age, y = R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +                    
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
  labs(
    x = "Max Age",
    y = expression(italic("R")),
    color = "Clonality rate"
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(face = "italic", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "italic", size = 13, margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Plot max_age vs Pareto_beta for each Clonality rate
max_age_vs_Pareto_beta_plot <- ggplot(merged_slope_parsed_df_extreme, aes(x = Max_Age, y = Pareto_beta_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 3, alpha = 0.7) +                    
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
  labs(
    x = "Max Age",
    y = expression(italic("Pareto " * beta)),
    color = "Clonality rate"
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(face = "italic", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "italic", size = 13, margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Extract legend from the first plot
legend_shared <- get_legend(
  max_age_vs_mean_fis_plot +
    theme(
      legend.position = "top",
      legend.justification = "center"
    )
)


# Combine the two plots into one figure with 2 columns

combined_max_age_plots <- plot_grid(
  max_age_vs_var_fis_plot, max_age_vs_mean_fis_plot,
  max_age_vs_R_plot, max_age_vs_Pareto_beta_plot,
  labels = c("a", "b", "c", "d"),
  ncol = 2,
  align = "v"
)


# Add the legend on top center
final_combined_max_age_plots <- plot_grid(
  legend_shared, 
  combined_max_age_plots, 
  ncol = 1, 
  rel_heights = c(0.1, 1)
)

final_combined_max_age_plots
#save the combined plot

ggsave(
  filename = "combined_max_age_plots.png",  # Save the combined plot
  plot = final_combined_max_age_plots,  # The combined plot object
  width = 25,  # Width of the saved image
  height = 18,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )

```


# Scatter plot for 4 replicates within specific clonerate
#0.1
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


# Select 4 random replicates from Clonality rate 0.1
set.seed(123)

random_reps_0.1 <- parsed_data %>%
  filter(clonalrate == "0.1") %>%
  distinct(Replicate) %>%
  pull(Replicate) %>%
  sample(4)

# Subset data for selected replicates at generation 6000
df_0.1 <- parsed_data %>%
  filter(clonalrate == "0.1", Generation == 6000, Replicate %in% random_reps_0.1)

# Generate scatter plots of Age vs Count (log scale) for each replicate

scatter_plots_0.1 <- lapply(random_reps_0.1, function(replicate) {
  df_sub <- df_0.1 %>% filter(Replicate == replicate)
  # Fit the linear model on log-transformed Count
  model <- lm(log10(Count) ~ Age, data = df_sub)
  slope <- coef(model)[2]
  r2 <- summary(model)$r.squared

  ggplot(df_sub, aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
             label = paste0("Slope = ", round(slope, 4), 
                            "\nR² = ", round(r2, 3)), 
             size = 5, color = "black", fontface = "bold") +
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
})



# Combine the plots into a single figure with 2 columns
combined_scatter_plots_0.1 <- plot_grid(plotlist = scatter_plots_0.1, 
                                        labels= c("a", "b", "c", "d"),
                                        ncol = 2)

# Display combined plot
combined_scatter_plots_0.1

# Save the figure as high-resolution PNG
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
# Take 4 randome replicates from the Clonality rate 

set.seed(123)  
random_reps_0.7 <- parsed_data %>%
  filter(clonalrate == "0.7") %>%
  pull(Replicate) %>%
  unique() %>%
  sample(4)

# View the selected replicates
df_0.7 <- parsed_data %>% filter(clonalrate == "0.7") %>% filter (Generation == 6000) %>% filter(Replicate %in% random_reps_0.7)

# create scatter plot for each replicate


scatter_plots_0.7 <- lapply(random_reps_0.7, function(replicate) {
  df_sub <- df_0.7 %>% filter(Replicate == replicate)
  # Fit the linear model on log-transformed Count
  model <- lm(log10(Count) ~ Age, data = df_sub)
  slope <- coef(model)[2]
  r2 <- summary(model)$r.squared

  ggplot(df_sub, aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
             label = paste0("Slope = ", round(slope, 4), 
                            "\nR² = ", round(r2, 3)), 
             size = 5, color = "black", fontface = "bold") +
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
})

# combine the plots into one
combined_scatter_plots_0.7 <- plot_grid(plotlist = scatter_plots_0.7, 
                                        labels= c("a", "b", "c", "d"),
                                        ncol = 2)

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

#0.99
```{r, echo=FALSE}
# Take 4 randome replicates from the Clonality rate 

set.seed(123)  
random_reps_0.9 <- parsed_data %>%
  filter(clonalrate == "0.99") %>%
  pull(Replicate) %>%
  unique() %>%
  sample(4)

# View the selected replicates
df_0.9 <- parsed_data %>% filter(clonalrate == "0.99") %>% filter (Generation == 6000) %>% filter(Replicate %in% random_reps_0.9)

# create scatter plot for each replicate
# create list to store plots


scatter_plots_0.9 <- lapply(random_reps_0.9, function(replicate) {
  df_sub <- df_0.9 %>% filter(Replicate == replicate)
  # Fit the linear model on log-transformed Count
  model <- lm(log10(Count) ~ Age, data = df_sub)
  slope <- coef(model)[2]
  r2 <- summary(model)$r.squared

  ggplot(df_sub, aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 6,
             label = paste0("Slope = ", round(slope, 4), 
                            "\nR² = ", round(r2, 3)), 
             size = 5, color = "black", fontface = "bold") +
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
})

# combine the plots into one
combined_scatter_plots_0.99 <- plot_grid(plotlist = scatter_plots_0.9, 
                                        labels= c("a", "b", "c", "d"),
                                        ncol = 2)

combined_scatter_plots_0.9

# Save the combined scatter plots
ggsave(
  filename = "scatter_plots_0.99.png",  
  plot = combined_scatter_plots_0.99,           
  width = 30,                  
  height = 20,    
  units = "cm",                   
  dpi = 1200                       
)
```

#0.999
```{r, echo=FALSE}
# Take 4 randome replicates from the Clonality rate XX
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
scatter_plots_0.999 <- lapply(random_reps_0.999, function(replicate) {
  df_sub <- df_0.999 %>% filter(Replicate == replicate)
  # Fit the linear model on log-transformed Count
  model <- lm(log10(Count) ~ Age, data = df_sub)
  slope <- coef(model)[2]
  r2 <- summary(model)$r.squared

  ggplot(df_sub, aes(x = Age, y = Count)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 5,
             label = paste0("Slope = ", round(slope, 4), 
                            "\nR² = ", round(r2, 3)), 
             size = 5, color = "black", fontface = "bold") +
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
})

# combine the plots into one
combined_scatter_plots_0.999 <- plot_grid(plotlist = scatter_plots_0.999, 
                                          labels= c("a", "b", "c", "d"),
                                          ncol = 2)
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
# Calculate means and standard deviations by Clonality rate
error_data <- merged_slope_parsed_data %>%
  group_by(clonalrate) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_R_Tot = mean(R_Tot),
    sd_R_Tot = sd(R_Tot),
    mean_Pareto_beta = mean(Pareto_beta_Tot),
    sd_Pareto_beta = sd(Pareto_beta_Tot),
    .groups = "drop"
  )


# Plot: Mean R^2 vs. Mean R with error bars no x-asix title and legend on line at top
r2_r <- ggplot(error_data, aes(x = mean_R2, y = mean_R_Tot, color = as.factor(clonalrate))) +
  geom_point(size = 4, shape = 3) +  # Cross-shaped points
  geom_errorbar(aes(ymin = mean_R_Tot - sd_R_Tot, ymax = mean_R_Tot + sd_R_Tot), width = 0.05) +
  geom_errorbarh(aes(xmin = mean_R2 - sd_R2, xmax = mean_R2 + sd_R2), height = 0.05) +
  geom_smooth(method = "lm", se = FALSE,) +
  labs(
    y = expression(italic(R)),
    color = "Clonality rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title.x = element_blank(),  # No x-axis title
    axis.text.x = element_blank(),  # No x-axis text
    axis.title = element_text(face = "italic", size = 12)
  )


# Plot: Mean R^2 vs. Mean Pareto β with error bars no legend
r2_beta <- ggplot(error_data, aes(x = mean_R2, y = mean_Pareto_beta, color = as.factor(clonalrate))) +
  geom_point(size = 4, shape = 3) +
  geom_errorbar(aes(ymin = mean_Pareto_beta - sd_Pareto_beta, ymax = mean_Pareto_beta + sd_Pareto_beta), width = 0.05) +
  geom_errorbarh(aes(xmin = mean_R2 - sd_R2, xmax = mean_R2 + sd_R2), height = 0.05) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = expression(R^2),
    y = expression(italic("Pareto ")*beta),
    color = "Clonality rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position="none",
    axis.title.x = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(face = "italic", size = 12)
  )

    


# Combine the two plots into one figure with 1 column
combined_r2_plots <- plot_grid(
  r2_r, r2_beta,
  labels = c("a", "b"),
  ncol = 1
)

# Add the legend on top center
legend_shared <- get_legend(
  r2_r +
    theme(
      legend.position = "top",
      legend.justification = "center"
    )
)

final_combined_r2_plots <- plot_grid(
  legend_shared, 
  combined_r2_plots, 
  ncol = 1, 
  rel_heights = c(0.1, 1)
)

final_combined_r2_plots

# Save the R^2 vs R plot
ggsave(
  filename = "R2_vs_R_Pareto_Tot.png",
  plot = final_combined_r2_plots,
  width = 15, height = 20, units = "cm", dpi = 1200, bg = "white"
)

```






## Ages per each clone rate 

```{r, echo=FALSE}

# Summarize data to count distinct ages per Clonality rate and replicate
age_counts <- parsed_data %>%
  group_by(clonalrate, Replicate) %>%
  summarize(Age_Count = n_distinct(Age), groups = "drop")  # Count unique ages per replicate within each Clonality rate

# Convert Clonality rate to a factor, so it appears in order on the x-axis
age_counts$clonalrate <- factor(age_counts$clonalrate, levels = sort(unique(age_counts$clonalrate)))

# Violin-Box Plot
age_class_plot <- ggplot(age_counts, aes(x = clonalrate, y = Age_Count, fill = clonalrate)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "#104E8B", size = 0.8) +  # Violin plot with transparency
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA, color = "black") +  # Overlay box plot
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +  # Add median point
  scale_fill_brewer(palette = "Blues") +  # Color palette for fill
  labs(
    x = "Clonality rate",
    y = "Number of age classes in replicate"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  theme(
    legend.position = "none",  # Remove legend for clarity
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  )

#save the plot
ggsave(
  filename = "age_classes_per_clonalrate.png",  # Save the plot
  plot = age_class_plot,  # The plot object
  width = 20,  # Width of the saved image
  height = 15,  # Height of the saved image
  units = "cm",  # Units for width and height
  dpi = 1200,  # Resolution of the saved image
  bg = "white" )

```


## oldest age per clonalrate

```{r, echo=FALSE}

# Summarize data to find the oldest age per Clonality rate and replicate
oldest_ages <- parsed_data %>%
  group_by(clonalrate, Replicate) %>%
  summarize(Max_Age = max(Age), .groups = "drop")  # Get the max age per replicate within each Clonality rate

# Convert Clonality rate to a factor for ordered plotting
oldest_ages$clonalrate <- factor(oldest_ages$clonalrate, levels = sort(unique(oldest_ages$clonalrate)))

#Plot the oldest age per Clonality rate using a box plot
ggplot(oldest_ages, aes(x = clonalrate, y = log(Max_Age), fill = clonalrate)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.color = "red", alpha = 0.7) +  # Box plot for oldest ages
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "darkblue") +  # Add jittered points for individual data
  labs(
    title = "Oldest Age Observed per Clonality rate",
    x = "Clonality rate (C value)",
    y = "Oldest Age in Replicate"
  ) +
  scale_fill_brewer(palette = "Blues") +  # Color palette
  theme_minimal(base_size = 15) +  # Minimal theme with larger font size
  theme(
    legend.position = "none",  # Remove legend for clarity
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  )

#max age at clone rate = 0.9
# Find the maximum age at Clonality rate 0.9
max_age_09 <- parsed_data %>%
  filter(clonalrate == "0.9") %>%
  summarize(Max_Age = max(Age))

```



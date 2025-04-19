###########################################################
########################################
#########
#########     SIMULATION        (Section 6)
#########
########################################
###########################################################


library(dplyr)
library(ggplot2)



#### AUXILIARY FUNCTIONS:

# 1) Generate a probability distribution for a fixed number of classes
#             n_classes (default=5) without ties 
#

generate_prob_dist <- function(n_classes = 5) {
  probs <- runif(n_classes)
  probs <- probs / sum(probs)  # Normalize
  # Slight jitter to avoid exact ties
  probs <- probs + runif(n_classes, 0, 1e-6)
  probs / sum(probs)
}


# 2) MAP (mode) and Ord-MAP (median) criteria for a probability distribution
#
#

get_MAP <- function(probs) {
  which.max(probs)
}

get_OrdMAP <- function(probs) {
  cum_probs <- cumsum(probs)
  which(cum_probs >= 0.5)[1]
}


# 3) Shannon entropy (for a probability distribution) 
#
#


shannon_entropy <- function(probs) {
  -sum(probs * log2(probs))
}


################################################################################
################################################################################
###########################
####    SIMULATION: PART 1
###########################

# Parameters
set.seed(123)
n_classes <- 6
n_sim_steps <- seq(10, 10000, by = 10)


results_list <- list()

# Simulation for each true class

for (true_class in 1:n_classes) {
  for (n_sim in n_sim_steps) {
    diffs <- numeric(n_sim)
    entropies <- numeric(n_sim)
    
    for (i in 1:n_sim) {
      probs <- generate_prob_dist(n_classes)
      pred_map <- get_MAP(probs)
      pred_ordmap <- get_OrdMAP(probs)
      mae_map <- abs(pred_map-true_class)
      mae_ordmap <- abs(pred_ordmap-true_class)
      diffs[i] <- sign(mae_map - mae_ordmap)
      entropies[i] <- shannon_entropy(probs)
    }
    
    diff_score <- (sum(diffs == 1) - sum(diffs == -1)) / n_sim * 100
    avg_entropy <- mean(entropies)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      true_class = factor(true_class),
      n_sim = n_sim,
      entropy = avg_entropy,
      diff_score = diff_score
    )
  }
}

results_df <- bind_rows(results_list)


### Plot

ggplot(results_df, aes(x = n_sim, y = diff_score, color = true_class)) +
  geom_line(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Simulations with 6 classes",
       x = "Number of Simulations",
       y = "% (MAP MAE > Ord-MAP MAE) - % (Ord-MAP MAE > MAP MAE)",
       color = "True Class") +
  theme_minimal()


################################################################################
################################################################################
###########################
####    SIMULATION: PART 2
###########################

# Parameters
set.seed(123)
n_sim <- 10000
n_classes <- 6


results <- data.frame()

# Simulation for each true class


for (true_class in 1:n_classes) {
  for (i in 1:n_sim) {
    p <- runif(n_classes)
    p <- p / sum(p)  # Normalizing to sum up to 1
    
    pred_map <- get_MAP(p)
    pred_ordmap <- get_OrdMAP(p)
    
    results <- rbind(results, data.frame(
      TrueClass = true_class,
      Entropy = shannon_entropy(p),
      MAE_MAP = abs(pred_map - true_class),
      MAE_Ord = abs(pred_ord - true_class),
      MAE_Diff = MAE_MAP - MAE_Ord
    ))
  }
}



### Group by entropy bins and true class

df_binned <- results %>%
  mutate(bin = cut(Entropy, breaks = 20)) %>%
  group_by(TrueClass, bin) %>%
  summarise(
    mean_entropy = mean(Entropy),
    mean_diff = mean(MAE_Diff),
    se = sd(MAE_Diff) / sqrt(n()),
    .groups = "drop"
  )



### Plot

ggplot(df_binned, aes(x = mean_entropy, y = mean_diff, color = factor(TrueClass))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_diff - se, ymax = mean_diff + se), width = 0.02, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Mean MAE difference vs. Shannon entropy per true class. 6 classes",
    x = "Shannon Entropy",
    y = "MAE(MAP) – MAE(Ord-MAP)",
    color = "True class"
  ) +
  theme_minimal()



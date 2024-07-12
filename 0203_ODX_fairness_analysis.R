# Load necessary libraries
library(ggplot2)
library(patchwork)
library(dplyr)
library(survival)
library(survminer)
library(tcltk2)

# Create a directory for saving plots
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Function to process a single sample of data
process_sample <- function(clin_TGx, sample_num, seed) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Prepare data: Create event and event_time columns
  clin_TGx <- clin_TGx %>%
    mutate(event = ifelse(!is.na(days_to_death), 1, 0),
           event_time = coalesce(days_to_death, days_to_last_follow_up)) %>%
    mutate(event_time = as.numeric(event_time),
           event = as.numeric(event))
  
  # Split data into train, calibration, and test sets
  samp <- sample(nrow(clin_TGx), 1092, replace = FALSE)
  train <- clin_TGx[samp[1:364], ]
  calib <- clin_TGx[samp[365:728], ]
  test <- clin_TGx[samp[729:1092], ]
  
  # Fit Cox proportional hazards model
  cox_model <- coxph(Surv(time = event_time, event = event) ~ ODX, data = train)
  summary(cox_model)
  
  # Determine cutpoints for ODX scores
  cutpoint_1 <- surv_cutpoint(data = train, time = "event_time", event = "event", variables = "ODX")
  temp <- train[train$event == 1, ]
  temp$event[temp$event_time > 365] <- 0
  cutpoint_2 <- surv_cutpoint(data = temp, time = "event_time", event = "event", variables = "ODX")
  cutpoints <- unlist(c(cutpoint_1$cutpoint[1], cutpoint_2$cutpoint[1]))
  
  # Assign predicted risk bins in calibration set
  calib <- calib %>%
    mutate(pred_risk_bin = ifelse(ODX <= cutpoints[1], 1, 
                                  ifelse(ODX <= cutpoints[2], 2, 3)))
  
  # Define true risk bins in calibration set
  calib <- calib %>%
    mutate(true_risk_bin = ifelse(event_time <= 365 & event == 1, 1, 
                                  ifelse(event_time <= 5*365 & event == 1, 2,
                                         ifelse(event_time <= 365 & event == 0, 2, 
                                                ifelse(event_time > 5*365 & event == 1, 3, 
                                                       ifelse(event_time > 365 & event == 0, 3, 4))))))
  
  # Calculate conformal scores
  calib <- calib %>%
    mutate(conf_score = ifelse(pred_risk_bin == true_risk_bin, 0,
                               ifelse(pred_risk_bin == 1 & true_risk_bin == 2, abs(ODX - cutpoints[1]), 
                                      ifelse(pred_risk_bin == 1 & true_risk_bin == 3, abs(ODX - cutpoints[2]),
                                             ifelse(pred_risk_bin == 2 & true_risk_bin == 1, abs(ODX - cutpoints[1]),
                                                    ifelse(pred_risk_bin == 2 & true_risk_bin == 3, abs(ODX - cutpoints[2]),
                                                           ifelse(pred_risk_bin == 3 & true_risk_bin == 1, abs(ODX - cutpoints[1]),
                                                                  ifelse(pred_risk_bin == 3 & true_risk_bin == 2, abs(ODX - cutpoints[2]), 4))))))))
  calib <- calib[!is.na(calib$conf_score), ]
  
  # Plot density of conformal scores
  plot_path <- paste0("plots/sample_", sample_num, "_conformal_density.png")
  png(plot_path)
  plot(density(calib$conf_score), main = paste("Density of Conformal Scores for Sample", sample_num))
  dev.off()
  
  # Calculate quantile for conformal scores
  alpha <- 0.05
  n <- nrow(calib)
  q_level <- ceiling((n + 1) * (1 - alpha)) / n
  qhat <- quantile(calib$conf_score, probs = q_level, type = 6, names = FALSE)
  
  # Calculate prediction intervals for the test set
  lo_bound <- test$ODX - qhat
  up_bound <- test$ODX + qhat
  lo_bound_bin <- ifelse(lo_bound <= cutpoints[1], 1, ifelse(lo_bound < cutpoints[2], 2, 3))
  up_bound_bin <- ifelse(up_bound <= cutpoints[1], 1, ifelse(up_bound < cutpoints[2], 2, 3))
  
  prediction_sets <- Map(list, lo_bound_bin, up_bound_bin)
  prediction_sets <- Map(function(x, y) {
    if ((y - x) == 0) {
      c(y)
    } else {
      if ((y - x) > 1) {
        c(x, y, 2)
      } else {
        c(x, y)
      }
    }
  }, lo_bound_bin, up_bound_bin)
  
  test <- test %>%
    mutate(prediction_sets = prediction_sets,
           pred_risk_bin = ifelse(ODX <= cutpoints[1], 1, ifelse(ODX <= cutpoints[2], 2, 3)),
           true_risk_bin = ifelse(event_time <= 365 & event == 1, 1, 
                                  ifelse(event_time <= 5*365 & event == 1, 2,
                                         ifelse(event_time <= 365 & event == 0, 2, 
                                                ifelse(event_time > 5*365 & event == 1, 3, 
                                                       ifelse(event_time > 365 & event == 0, 3, 4)))))
    )
  
  # Combine true and predicted risk bins for visualization
  combined_data <- data.frame(
    value = c(test$true_risk_bin, test$pred_risk_bin),
    variable = factor(rep(c("True Risk", "Predicted Risk"), each = nrow(test)))
  )
  
  # Plot histogram with dodged bars
  p <- ggplot(combined_data, aes(x = value, fill = variable)) +
    geom_bar(position = "dodge", alpha = 0.5) +
    labs(x = "Risk Bin", y = "Frequency", fill = "Variable") +
    theme_minimal() +
    ggtitle(paste("Risk Bin Distribution for Sample", sample_num))
  
  return(list(cox_model = cox_model, fit = fit, cutpoints = cutpoints, combined_data = combined_data, conf_score = calib$conf_score, test = test))
}

# Main function to run the process for 100 samples
run_analysis <- function(clin_TGx, n_samples = 100) {
  results <- vector("list", n_samples)
  seeds <- sample(1:10000, n_samples) # Generate unique seeds for each sample
  for (i in 1:n_samples) {
    cat("Processing sample", i, "\n")
    results[[i]] <- process_sample(clin_TGx, i, seeds[i])
  }
  return(results)
}

# Run the analysis
results <- run_analysis(clin_TGx)

# Collect combined data from all samples
all_combined_data <- do.call(rbind, lapply(results, function(x) x$combined_data))

# Plot the averages of all samples
p_avg <- ggplot(all_combined_data, aes(x = value, fill = variable)) +
  geom_bar(position = "dodge", alpha = 0.5) +
  labs(x = "Risk Bin", y = "Frequency", fill = "Variable") +
  theme_minimal() +
  ggtitle("Average Risk Bin Distribution Across All Samples")

ggsave(filename = "plots/average_risk_bins.png", plot = p_avg)

# Collect conformal scores from all samples
all_conf_scores <- unlist(lapply(results, function(x) x$conf_score))

# Plot the density of conformal scores
p_conf_scores <- ggplot(data.frame(conf_score = all_conf_scores), aes(x = conf_score)) +
  geom_density() +
  labs(x = "Conformal Score", y = "Density") +
  theme_minimal() +
  ggtitle("Density of Conformal Scores Across All Samples")

ggsave(filename = "plots/density_conformal_scores.png", plot = p_conf_scores)

# Display the plots
print(p_avg)
print(p_conf_scores)

# Analysis code: Compare accuracy, prediction set lengths, and conformal scores by race group



### FAIRNESS ANALYSIS --------------------------------------------------

# Create accuracy_df from results list
accuracy_list <- lapply(results, function(res) {
  test <- res$test
  data.frame(
    race = test$race,
    accuracy = ifelse(test$pred_risk_bin == test$true_risk_bin, 1, 0),
    lengths_prediction_sets = lengths(test$prediction_sets),
    random_risk_bin = sample(1:3, nrow(test), replace = TRUE)  # Random risk allocation
  )
})

# Combine all the data into a single data frame
accuracy_df <- do.call(rbind, accuracy_list)

# Save directory for plots
plot_dir <- "plots"

# Create the directory if it doesn't exist
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Plot 1: Mean Lengths of Prediction Sets by Race
mean_by_subgroup <- data.frame(aggregate(lengths_prediction_sets ~ race, data = accuracy_df, FUN = mean))
colnames(mean_by_subgroup)[2] <- "mean_length"

plot_mean_lengths <- ggplot(mean_by_subgroup, aes(x = race, y = mean_length, fill = race)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Lengths of Prediction Sets by Race",
       x = "Race",
       y = "Mean Length",
       fill = "Race") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(plot_dir, "mean_lengths_by_race.png"), plot = plot_mean_lengths)

# Plot 2: Occurrences of Unique Lengths of Prediction Sets by Race
occurrences_by_race <- accuracy_df %>%
  group_by(race, lengths_prediction_sets) %>%
  count() %>%
  ungroup() %>%
  rename(occurrences = n)

plot_occurrences <- ggplot(occurrences_by_race, aes(x = lengths_prediction_sets, y = occurrences, fill = race)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Occurrences of Unique Lengths of Prediction Sets by Race",
       x = "Length of Prediction Sets",
       y = "Occurrences",
       fill = "Race") +
  theme_minimal()

ggsave(filename = file.path(plot_dir, "occurrences_by_race.png"), plot = plot_occurrences)

# Plot 3: Comparison of Accuracy to Lengths of Conformal Sets by Subgroup
plot_accuracy_vs_conformal <- ggplot(accuracy_df, aes(x = lengths_prediction_sets, y = accuracy, color = race)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of Accuracy to Lengths of Conformal Sets by Subgroup",
       x = "Length of Conformal Sets",
       y = "Accuracy",
       color = "Race") +
  theme_minimal()

ggsave(filename = file.path(plot_dir, "accuracy_vs_conformal.png"), plot = plot_accuracy_vs_conformal)


# plot 4:
# Calculate whether the predicted bin included the true bin
accuracy_df <- accuracy_df %>%
  mutate(includes_true_bin = ifelse(pred_risk_bin == true_risk_bin, "Included", "Not Included"))

# Plot the average length of prediction sets by subgroups with inclusion of true bin
plot_avg_length_with_inclusion <- accuracy_df %>%
  group_by(race, includes_true_bin) %>%
  summarise(avg_length = mean(lengths_prediction_sets)) %>%
  ggplot(aes(x = race, y = avg_length, fill = includes_true_bin)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Average Length of Prediction Sets by Subgroups with True Bin Inclusion",
       x = "Race",
       y = "Average Length",
       fill = "Includes True Bin") +
  theme_minimal()

ggsave(filename = file.path(plot_dir, "avg_length_with_inclusion.png"), plot = plot_avg_length_with_inclusion, width = 10, height = 8, dpi = 300)

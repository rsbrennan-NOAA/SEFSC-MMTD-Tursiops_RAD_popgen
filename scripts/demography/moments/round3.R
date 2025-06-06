# analyze initial run of moments, to choose params for next round

library(tidyverse)

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble()

# then add all results to this
for (k in 1:8){
  data <- read_delim(paste0("analysis/moments/round3/model_0",k,"_fmin_downSample_summary.csv"),
                     show_col_types = FALSE,
                     delim=",")
  data$model <- paste0("0", k)
  tmpdat <- data[,c("model", "Run", "LogLikelihood")]
  all_data <- rbind(all_data,tmpdat)
  
}

p <- ggplot(all_data, aes(x=model, y=LogLikelihood)) +
  geom_boxplot(outliers = TRUE) +
  ylim(min(all_data$LogLikelihood), 0)
p

ggsave("figures/moments/round-3_logLikelihoods.pdf", p, h=4, w=5)

# find highest ll for each model. 
best_models <- all_data %>%
  group_by(model) %>%
  slice_max(LogLikelihood, n = 1) %>%
  arrange(model)



# then add all results to this
for (k in 1:8){
  data <- read_delim(paste0("analysis/moments/round3/model_0",k,"_fmin_downSample_summary.csv"),
                     show_col_types = FALSE,
                     delim=",")
  best_run <- best_models[k, ]
  model_num <- best_models$model[k]
  
  
  best_params <- data[data$Run == best_run$Run,]
  #  write_csv(best_params, paste0("analysis/moments/round2/best_params_model_", model_num, ".csv"))
  
  cat("Saved best parameters for model", model_num, "from run", best_run$Run, "\n")
  cat("  Log-likelihood:", best_params$LogLikelihood, "\n")
  
}


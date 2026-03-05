# ... [previous code in dataset loop]

  # Store the results for the current dataset
  LASSO_results[[dataStr]] <- dataset_results

  # Export CSVs for LASSO_improved results
  if (csv) {
    dataset_csv_dir <- file.path("./results", dataStr, "LASSO_improved")
    if (!dir.exists(dataset_csv_dir)) {
      dir.create(dataset_csv_dir, recursive = TRUE)
    }
    
    for (outcome_name in names(dataset_results)) {
      outcome_results <- dataset_results[[outcome_name]]
      
      # Create combined dataframe for all subsets
      combined_df <- data.frame()
      all_predictors <- metrics  # From your metrics vector
      
      for (subset_name in names(outcome_results)) {
        subset_res <- outcome_results[[subset_name]]
        coefs <- subset_res$selected_predictors
        
        # Create row with all predictors
        row <- data.frame(matrix(NA, ncol = length(all_predictors), 
                               dimnames = list(NULL, all_predictors)))
        row$subset <- subset_name
        
        # Fill coefficients
        if (nrow(coefs) > 0) {
          for (i in 1:nrow(coefs)) {
            pred <- as.character(coefs$predictor[i])
            row[1, pred] <- coefs$coefficient[i]
          }
        }
        
        combined_df <- rbind(combined_df, row)
      }
      
      # Reorder columns
      combined_df <- combined_df[, c("subset", all_predictors)]
      
      # Save CSV
      write.csv(combined_df, 
              file.path(dataset_csv_dir, paste0(outcome_name, "_coefficients.csv")), 
              row.names = FALSE)
    }
  }

# ... [rest of original script]
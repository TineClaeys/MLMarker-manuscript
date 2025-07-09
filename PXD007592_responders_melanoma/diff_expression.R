# ---- Install & Load Required Packages ----
if (!requireNamespace("MSqRob2", quietly = TRUE)) {
    install.packages("MSqRob2", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("MSqRob2", "QFeatures", "limma", "tidyverse")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msqrob2")
library(MSqRob2)
library(tidyverse)
library(limma)  
library(QFeatures)
library(readr)

# ---- Step 1: Load Data ----
# Replace 'your_dataset.csv' with the actual dataset path
data <- read_csv("/home/compomics/git/Projects/MLMarker_paper/PXD007592_responders_melanoma/NSAF_values.csv")

# Ensure correct column names
colnames(data)[1] <- "filename"

# ---- Step 2: Transform Data to Long Format ----
data_long <- data %>%
  pivot_longer(cols = -c(filename, response_state, cluster), 
               names_to = "Protein", 
               values_to = "NSAF")

# Log2 transform NSAF values to avoid log(0) issues
data_long <- data_long %>%
  mutate(logNSAF = log2(NSAF + 1e-6))

# ---- Step 3: Convert to QFeatures Object ----
# Reshape data for QFeatures (wide format with proteins as rows)
assay_data <- data_long %>%
  pivot_wider(names_from = filename, values_from = logNSAF)

# Extract feature metadata
rowdata <- data.frame(Protein = assay_data$Protein)
rownames(assay_data) <- assay_data$Protein
assay_data <- assay_data %>% select(-Protein)

# Create colData for sample annotations
coldata <- data %>%
  select(response_state, cluster) %>%
  distinct()
rownames(coldata) <- coldata$cluster

# Create QFeatures object
qf <- QFeatures(
  assays = list(logNSAF = as.matrix(assay_data)),
  colData = coldata
)

# ---- Step 4: Perform Differential Expression Analysis ----
run_msqrob2 <- function(qfeatures_obj, grouping_col, output_file) {
  # Apply MSqRob2 with the selected grouping factor
  qfeatures_obj <- addAssay(qfeatures_obj, name = "msqrob", assay = assay(qfeatures_obj, "logNSAF"))

  # Define the model formula
  formula <- as.formula(paste("~", grouping_col))
  
  # Fit the MSqRob2 model
  qfeatures_obj <- msqrob(qfeatures_obj, i = "msqrob", formula = formula)

  # Extract results
  results <- rowData(qfeatures_obj[["msqrob"]])$coefficients
  results_df <- as.data.frame(results) %>%
    arrange(p.value)  # Sort by significance

  # Save results
  write.csv(results_df, output_file, row.names = FALSE)

  # ---- Step 5: Generate Volcano Plot ----
  ggplot(results_df, aes(x = logFC, y = -log10(p.value))) +
    geom_point(aes(color = p.value < 0.05)) +
    scale_color_manual(values = c("black", "red")) +
    labs(title = paste("Volcano Plot for", grouping_col, "Comparison"),
         x = "Log2 Fold Change",
         y = "-Log10 P-value") +
    theme_minimal()
}

# ---- Run Analysis for Response State (Good vs Poor) ----
run_msqrob2(qf, "response_state", "msqrob2_response_state_results.csv")

# ---- Run Analysis for Cluster (4 Groups) ----
run_msqrob2(qf, "cluster", "msqrob2_cluster_results.csv")

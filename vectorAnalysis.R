# Load required libraries
library(dunn.test)
library(tidyverse)
library(reshape2)

# Load the CSV
df <- read.csv("vector_lengths_maxConds.csv")  
#df <- read.csv("vector_lengths_halfAcs.csv") 

# Reshape to long format
df_long <- melt(df, variable.name = "Condition", value.name = "Value")

# Kruskal-Wallis test
kruskal_test <- kruskal.test(Value ~ Condition, data = df_long)
print(kruskal_test)

# Post-hoc Dunn's test
dunn_result <- dunn.test(df_long$Value, df_long$Condition, method = "bonferroni")
print(dunn_result)


# Extract results
comparison_labels <- dunn_result$comparisons
adjusted_pvals <- dunn_result$P.adjusted

# Combine into a data frame
dunn_df <- data.frame(Comparison = comparison_labels,
                      P_adj = adjusted_pvals)

# Filter significant comparisons (adjusted p < 0.05)
significant_results <- dunn_df %>% filter(P_adj < 0.05)

# Print only significant pairwise comparisons
cat("\nSignificant pairwise comparisons (Bonferroni-adjusted p < 0.05):\n")
print(significant_results)



#############################################################################################################
#                                                                                                           #
#                                 Script for downstream analysis of MaxQuant SILAC data                     #
#                                                                                                           #
#############################################################################################################

# Setting up the environment
install.packages("ggplot2")
library("ggplot2")

#############################################################################################################

#1. Load MaxQuant (MQ) results
#Load the results from MQ together with it’s headers.
proteinGroups <- "proteinGroups.txt"

proteinGroups <- read.csv(proteinGroupsproteinGroups, sep = "\t")

#############################################################################################################

#2. Filter
#Remove all contaminants. The table proteinGroups.txt from MQ has a column Contaminant.
proteinGroups_filtered <- subset(proteinGroups, Potential.contaminant != "+")
nrow(proteinGroups_filtered)

#We also want to filter out all proteins that are identified by modification side and those that are derived from the reverse database
proteinGroups_filtered <- subset(proteinGroups_filtered, Only.identified.by.site != "+" & Reverse != "+")
nrow(proteinGroups_filtered)
#Now we should have 6893 protein groups left. Check the number of proteins with nrow(proteinGroups_filtered).

#############################################################################################################

#3. Analysis 1
#3.1. Histogram of H/L ratios
#There are 6 different columns containing H/L ratios. Plot a histogram for each of them individually by adapting the variable H_L_Ratio_column:
H_L_Ratio_1F <- "Ratio.H.L.1F"
hist_H_L_Ratio_1F <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_1F)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L 1F",
       y = "Counts")
print(hist_H_L_Ratio_1F)

H_L_Ratio_1R <- "Ratio.H.L.1R"
hist_H_L_Ratio_1R <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_1R)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L 1R",
       y = "Counts")
print(hist_H_L_Ratio_1R)

H_L_Ratio_2F <- "Ratio.H.L.2F"
hist_H_L_Ratio_2F <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_2F)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L 2F",
       y = "Counts")
print(hist_H_L_Ratio_2F)

H_L_Ratio_2R <- "Ratio.H.L.2R"
hist_H_L_Ratio_2R <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_2R)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L 2R",
       y = "Counts")
print(hist_H_L_Ratio_2R)

H_L_Ratio_3F <- "Ratio.H.L.3F"
hist_H_L_Ratio_3F <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_3F)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L 3F",
       y = "Counts")
print(hist_H_L_Ratio_3F)

H_L_Ratio_3R <- "Ratio.H.L.3R"
hist_H_L_Ratio_3R <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_3R)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L 3R",
       y = "Counts")
print(hist_H_L_Ratio_3R)

# Arrange plots in a grid
combined_histograms <- grid.arrange(hist_H_L_Ratio_1F, hist_H_L_Ratio_1R, hist_H_L_Ratio_2F,
                                    hist_H_L_Ratio_2R, hist_H_L_Ratio_3F, hist_H_L_Ratio_3R,
                                    ncol = 2, nrow = 3)

# Save the combined plot as a PNG file
ggsave("combined_histograms.png", plot = combined_histograms, width = 10, height = 8)

#3.2. Histogram of normalized H/L ratios
H_L_Ratio_normalized_1F <- "Ratio.H.L.normalized.1F"
hist_H_L_Ratio_normalized_1F <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_normalized_1F)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L normalized 1F",
       y = "Counts")
print(hist_H_L_Ratio_normalized_1F)

H_L_Ratio_normalized_1R <- "Ratio.H.L.normalized.1R"
hist_H_L_Ratio_normalized_1R <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_normalized_1R)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L normalized 1R",
       y = "Counts")
print(hist_H_L_Ratio_normalized_1R)

H_L_Ratio_normalized_2F <- "Ratio.H.L.normalized.2F"
hist_H_L_Ratio_normalized_2F <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_normalized_2F)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L normalized 2F",
       y = "Counts")
print(hist_H_L_Ratio_normalized_2F)

H_L_Ratio_normalized_2R <- "Ratio.H.L.normalized.2R"
hist_H_L_Ratio_normalized_2R <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_normalized_2R)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L normalized 2R",
       y = "Counts")
print(hist_H_L_Ratio_normalized_2R)

H_L_Ratio_normalized_3F <- "Ratio.H.L.normalized.3F"
hist_H_L_Ratio_normalized_3F <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_normalized_3F)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L normalized 3F",
       y = "Counts")
print(hist_H_L_Ratio_normalized_3F)

H_L_Ratio_normalized_3R <- "Ratio.H.L.normalized.3R"
hist_H_L_Ratio_normalized_3R <- ggplot(proteinGroups_filtered, aes_string(x = H_L_Ratio_normalized_3R)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = "Ratio H/L normalized 3R",
       y = "Counts")
print(hist_H_L_Ratio_normalized_3R)

# Arrange plots in a grid
combined_normalized_histograms <- grid.arrange(hist_H_L_Ratio_normalized_1F, hist_H_L_Ratio_normalized_1R, hist_H_L_Ratio_normalized_2F,
                                               hist_H_L_Ratio_normalized_2R, hist_H_L_Ratio_normalized_3F, hist_H_L_Ratio_normalized_3R,
                                               ncol = 2, nrow = 3)

# Save the combined plot as a PNG file
ggsave("combined_normalized_histograms.png", plot = combined_normalized_histograms, width = 10, height = 8)

#3.3. Summary statistics
#To get a basic summary statistics of our results table proteinGroups_filtered we can simply call:

summary(proteinGroups_filtered)

#But we might want to call it only on the columns containing H/L Ratios:
H_L_Ratio_columns <- grep("Ratio.H.L", colnames(proteinGroups_filtered), value = TRUE)
summary(proteinGroups_filtered[, H_L_Ratio_columns])

#3.4. Scatterplot of H/L ratios
#Let’s create a scatterplot of H/L ratios for all normalized H/L ratios:

scatterplot_H_L_Ratio_normalized_1F <- ggplot(proteinGroups_filtered, aes_string(x="id", y = H_L_Ratio_normalized_1F)) +
  geom_point(alpha = 0.3)
print(scatterplot_H_L_Ratio_normalized_1F)

scatterplot_H_L_Ratio_normalized_1R <- ggplot(proteinGroups_filtered, aes_string(x="id", y = H_L_Ratio_normalized_1R)) +
  geom_point(alpha = 0.3)
print(scatterplot_H_L_Ratio_normalized_1R)

scatterplot_H_L_Ratio_normalized_2F <- ggplot(proteinGroups_filtered, aes_string(x="id", y = H_L_Ratio_normalized_2F)) +
  geom_point(alpha = 0.3)
print(scatterplot_H_L_Ratio_normalized_2F)

scatterplot_H_L_Ratio_normalized_2R <- ggplot(proteinGroups_filtered, aes_string(x="id", y = H_L_Ratio_normalized_2R)) +
  geom_point(alpha = 0.3)
print(scatterplot_H_L_Ratio_normalized_2R)

scatterplot_H_L_Ratio_normalized_3F <- ggplot(proteinGroups_filtered, aes_string(x="id", y = H_L_Ratio_normalized_3F)) +
  geom_point(alpha = 0.3)
print(scatterplot_H_L_Ratio_normalized_3F)

scatterplot_H_L_Ratio_normalized_3R <- ggplot(proteinGroups_filtered, aes_string(x="id", y = H_L_Ratio_normalized_3R)) +
  geom_point(alpha = 0.3)
print(scatterplot_H_L_Ratio_normalized_3R)

#Note the distribution of ratios. Is it easy to evaluate the data?

#############################################################################################################

#4. Log transformation
# Define the patterns you want to search for
patterns <- c("Ratio.H.L.1F", "Ratio.H.L.1R", "Ratio.H.L.2F", "Ratio.H.L.2R", "Ratio.H.L.3F", "Ratio.H.L.3R")

# Use grepl to check for matches
matches <- grepl(paste(patterns, collapse = "|"), colnames(proteinGroups_filtered))

# Extract the matching column names
H_L_Ratio_columns <- colnames(proteinGroups_filtered)[matches]

# Display the matches
print(H_L_Ratio_columns)

#Make log2 transformation of the H/L Ratios columns:
H_L_Ratio_log2 <- log2(proteinGroups_filtered[, H_L_Ratio_columns])

# add a "log2" to the column names and add them back to the table
colnames(H_L_Ratio_log2) <- paste0(colnames(H_L_Ratio_log2), ".log2")
proteinGroups_filtered_logtransformed <- cbind(proteinGroups_filtered, H_L_Ratio_log2)

#4.1. Histogram of logtransformed H/L ratios
H_L_Ratio_1F_log2 <- "Ratio.H.L.1F.log2"
hist_H_L_Ratio_1F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_1F_log2)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L 1F)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_1F_log2)

H_L_Ratio_1R_log2 <- "Ratio.H.L.1R.log2"
hist_H_L_Ratio_1R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_1R_log2)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L 1R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_1R_log2)

H_L_Ratio_2F_log2 <- "Ratio.H.L.2F.log2"
hist_H_L_Ratio_2F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_2F_log2)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L 2F)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_2F_log2)

H_L_Ratio_2R_log2 <- "Ratio.H.L.2R.log2"
hist_H_L_Ratio_2R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_2R_log2)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L 2R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_2R_log2)

H_L_Ratio_3F_log2 <- "Ratio.H.L.3F.log2"
hist_H_L_Ratio_3F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_3F_log2)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L 3F)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_3F_log2)

H_L_Ratio_3R_log2 <- "Ratio.H.L.3R.log2"
hist_H_L_Ratio_3R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_3R_log2)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L 3R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_3R_log2)

# Arrange plots in a grid
combined_log2_histograms <- grid.arrange(hist_H_L_Ratio_1F_log2, hist_H_L_Ratio_1R_log2, hist_H_L_Ratio_2F_log2,
                                         hist_H_L_Ratio_2R_log2, hist_H_L_Ratio_3F_log2, hist_H_L_Ratio_3R_log2,
                                         ncol = 2, nrow = 3)

# Save the combined plot as a PNG file
ggsave("combined_log2_histograms.png", plot = combined_log2_histograms, width = 10, height = 8)

#Make log2 transformation of the normalized H/L Ratios columns:
H_L_Ratio_normalized_columns <- grep('Ratio.H.L.normalized', colnames(proteinGroups_filtered_logtransformed), value = TRUE)

H_L_Ratio_normalized_log2 <- log2(proteinGroups_filtered_logtransformed[, H_L_Ratio_normalized_columns])

# add a "log2" to the column names and add them back to the table
colnames(H_L_Ratio_normalized_log2) <- paste0(colnames(H_L_Ratio_normalized_log2), ".log2")
proteinGroups_filtered_logtransformed <- cbind(proteinGroups_filtered_logtransformed, H_L_Ratio_normalized_log2)

colnames(proteinGroups_filtered_logtransformed)

#4.2. Histogram of normalized and logtransformed H/L ratios
H_L_Ratio_normalized_1F_log2 <- "Ratio.H.L.normalized.1F.log2"
hist_H_L_Ratio_normalized_1F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_normalized_1F_log2)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L normalized 1F)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_1F_log2)

H_L_Ratio_normalized_1R_log2 <- "Ratio.H.L.normalized.1R.log2"
hist_H_L_Ratio_normalized_1R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_normalized_1R_log2)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L normalized 1R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_1R_log2)

H_L_Ratio_normalized_2F_log2 <- "Ratio.H.L.normalized.2F.log2"
hist_H_L_Ratio_normalized_2F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_normalized_2F_log2)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L normalized 2F)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_2F_log2)

H_L_Ratio_normalized_2R_log2 <- "Ratio.H.L.normalized.2R.log2"
hist_H_L_Ratio_normalized_2R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_normalized_2R_log2)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L normalized 2R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_2R_log2)

H_L_Ratio_normalized_3F_log2 <- "Ratio.H.L.normalized.3F.log2"
hist_H_L_Ratio_normalized_3F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_normalized_3F_log2)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L normalized 3F)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_3F_log2)

H_L_Ratio_normalized_3R_log2 <- "Ratio.H.L.normalized.3R.log2"
hist_H_L_Ratio_normalized_3R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = H_L_Ratio_normalized_3R_log2)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("Log"[2]*"(Ratio H/L normalized 3R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_3R_log2)

# Arrange plots in a grid
combined_log2_normalized_histograms <- grid.arrange(hist_H_L_Ratio_normalized_1F_log2, hist_H_L_Ratio_normalized_1R_log2, hist_H_L_Ratio_normalized_2F_log2,
                                                    hist_H_L_Ratio_normalized_2R_log2, hist_H_L_Ratio_normalized_3F_log2, hist_H_L_Ratio_normalized_3R_log2,
                                                    ncol = 2, nrow = 3)

# Save the combined plot as a PNG file
ggsave("combined_log2_normalized_histograms.png", plot = combined_log2_normalized_histograms, width = 10, height = 8)

#############################################################################################################

#5. Analysis 2
#5.1. Scatterplot of log2
#Create the same scatterplots as before, but this time on the log2 data:

scatterplot_H_L_Ratio_normalized_1F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_1F_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3)
print(scatterplot_H_L_Ratio_normalized_1F_log2)

scatterplot_H_L_Ratio_normalized_1R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_1R_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3)
print(scatterplot_H_L_Ratio_normalized_1R_log2)

scatterplot_H_L_Ratio_normalized_2F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_2F_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3)
print(scatterplot_H_L_Ratio_normalized_2F_log2)

scatterplot_H_L_Ratio_normalized_2R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_2R_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3)
print(scatterplot_H_L_Ratio_normalized_2R_log2)

scatterplot_H_L_Ratio_normalized_3F_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_3F_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3)
print(scatterplot_H_L_Ratio_normalized_3F_log2)

scatterplot_H_L_Ratio_normalized_3R_log2 <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_3R_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3)
print(scatterplot_H_L_Ratio_normalized_3R_log2)

#Select proteins above a certain threshold (e.g. 1) in replicate 1 and look if they’re also high in the other replicates:

H_L_Ratio_normalized_1F_log2_above_threshold_1 <- subset(proteinGroups_filtered_logtransformed, Ratio.H.L.normalized.1F.log2 > 1)

scatterplot_H_L_Ratio_normalized_1F_log2_threshold <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_1F_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3) +
  geom_point(data = H_L_Ratio_normalized_1F_log2_above_threshold_1, aes_string(x = "id", y = H_L_Ratio_normalized_1F_log2), col = "red")  
print(scatterplot_H_L_Ratio_normalized_1F_log2_threshold)

#Now let’s plot the corresponding H/L ratios for the second replicate:

H_L_Ratio_normalized_2F_log2_above_threshold_1 <- subset(proteinGroups_filtered_logtransformed, Ratio.H.L.normalized.2F.log2 > 1)

scatterplot_H_L_Ratio_normalized_2F_log2_threshold <- ggplot(proteinGroups_filtered_logtransformed, aes_string(x = "id", y = H_L_Ratio_normalized_2F_log2)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.3) +
  geom_point(data = H_L_Ratio_normalized_2F_log2_above_threshold_1, aes_string(x = "id", y = H_L_Ratio_normalized_2F_log2), col = "blue")  
print(scatterplot_H_L_Ratio_normalized_2F_log2_threshold)

#5.2. Multiscatterplot and Pearson correlation
#We can look at the correlations between the different replicates using a multiscatterplot:

H_L_Ratio_normalized_log2_columns <- grep("Ratio.H.L.normalized.*.log2", colnames(proteinGroups_filtered_logtransformed), value = TRUE)
H_L_Ratio_normalized_log2_columns

# Iterate over selected columns and remove non-numeric values
for (column_name in H_L_Ratio_normalized_log2_columns) {
  # Remove non-numeric characters and convert to numeric
  proteinGroups_filtered_logtransformed[[column_name]] <- as.numeric(gsub("[^0-9.-]", "", as.character(proteinGroups_filtered_logtransformed[[column_name]])))
}

# Now, you can use the cleaned numeric columns for further analysis
pairs(proteinGroups_filtered_logtransformed[, H_L_Ratio_normalized_log2_columns])

#The Pearson correlation can be computed with the following command:

cor(proteinGroups_filtered_logtransformed[,H_L_Ratio_normalized_log2_columns], method = "pearson", use = "pairwise.complete.obs") 

#############################################################################################################

#6. Filter more
#We filter on the number of peptides observed per protein groups. We’re keeping only proteins which have > 3 Peptides and > 1 Unique peptides:

proteinGroups_strict_filtered_logtransformed <- subset(proteinGroups_filtered_logtransformed, Peptides > 3 & Unique.peptides > 1)
nrow(proteinGroups_strict_filtered_logtransformed)
#Now we only have 5094 protein groups left

#And we make again a multi scatterplot and compute the Pearson correlation:

# multi scatterplot
pairs(proteinGroups_strict_filtered_logtransformed[, H_L_Ratio_normalized_log2_columns])

# Pearson correlation
cor(proteinGroups_strict_filtered_logtransformed[, H_L_Ratio_normalized_log2_columns], method = "pearson", use = "pairwise.complete.obs") 
#Is the correlation better now?

#############################################################################################################

#7.1. Transform H/R ratio of the reverse samples (-x)
# Define the patterns you want to search for
patterns <- c("Ratio.H.L.1R.log2", "Ratio.H.L.2R.log2", "Ratio.H.L.3R.log2")

# Use grepl to check for matches
matches_reverse <- grepl(paste(patterns, collapse = "|"), colnames(proteinGroups_strict_filtered_logtransformed))

# Extract the matching column names
H_L_Ratios_R_log2_negx <- colnames(proteinGroups_strict_filtered_logtransformed)[matches_reverse]

# Negate their values
H_L_Ratios_R_log2_negx <- -(proteinGroups_strict_filtered_logtransformed[,H_L_Ratios_R_log2_negx])

# add a "-x" to the column names and add them back to the table
colnames(H_L_Ratios_R_log2_negx) <- paste0(colnames(H_L_Ratios_R_log2_negx), ".negx")
proteinGroups_strict_filtered_logtransformed_negx <- cbind(proteinGroups_strict_filtered_logtransformed, H_L_Ratios_R_log2_negx)

#7.2. Histogram of logtransformed and (-x) transformed reverse H/L ratios
H_L_Ratio_1R_log2_negx <- "Ratio.H.L.1R.log2.negx"
hist_H_L_Ratio_1R_log2_negx <- ggplot(proteinGroups_strict_filtered_logtransformed_negx, aes_string(x = H_L_Ratio_1R_log2_negx)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("-"*"Log"[2]*"(Ratio H/L 1R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_1R_log2_negx)

H_L_Ratio_2R_log2_negx <- "Ratio.H.L.2R.log2.negx"
hist_H_L_Ratio_2R_log2_negx <- ggplot(proteinGroups_strict_filtered_logtransformed_negx, aes_string(x = H_L_Ratio_2R_log2_negx)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("-"*"Log"[2]*"(Ratio H/L 2R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_2R_log2_negx)

H_L_Ratio_3R_log2_negx <- "Ratio.H.L.3R.log2.negx"
hist_H_L_Ratio_3R_log2_negx <- ggplot(proteinGroups_strict_filtered_logtransformed_negx, aes_string(x = H_L_Ratio_3R_log2_negx)) +
  geom_histogram(fill = "#7DD082", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("-"*"Log"[2]*"(Ratio H/L 3R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_3R_log2_negx)

# Arrange plots in a grid
combined_neg_log2_histograms <- grid.arrange(hist_H_L_Ratio_1F_log2, hist_H_L_Ratio_1R_log2_negx, hist_H_L_Ratio_2F_log2,
                                             hist_H_L_Ratio_2R_log2_negx, hist_H_L_Ratio_3F_log2, hist_H_L_Ratio_3R_log2_negx,
                                             ncol = 2, nrow = 3)

# Save the combined plot as a PNG file
ggsave("combined_neg_log2_histograms.png", plot = combined_neg_log2_histograms, width = 10, height = 8)

#7.3. Transform normalized H/R ratio of the reverse samples (-x)
H_L_Ratio_normalized_R_log2 <- grep('Ratio.H.L.normalized.*R.log2', colnames(proteinGroups_strict_filtered_logtransformed_negx), value=TRUE)

# Extract columns and negate their values
H_L_Ratios_normalized_R_log2_negx <- -(proteinGroups_strict_filtered_logtransformed_negx[,H_L_Ratio_normalized_R_log2])

# add a "-x" to the column names and add them back to the table
colnames(H_L_Ratios_normalized_R_log2_negx) <- paste0(colnames(H_L_Ratios_normalized_R_log2_negx), ".negx")
proteinGroups_strict_filtered_logtransformed_negx <- cbind(proteinGroups_strict_filtered_logtransformed_negx, H_L_Ratios_normalized_R_log2_negx)

#7.4. Histogram of logtransformed and (-x) transformed reverse H/L ratios
H_L_Ratio_normalized_1R_log2_negx <- "Ratio.H.L.normalized.1R.log2.negx"
hist_H_L_Ratio_normalized_1R_log2_negx <- ggplot(proteinGroups_strict_filtered_logtransformed_negx, aes_string(x = H_L_Ratio_normalized_1R_log2_negx)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("-"*"Log"[2]*"(Ratio H/L normalized 1R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_1R_log2_negx)

H_L_Ratio_normalized_2R_log2_negx <- "Ratio.H.L.normalized.2R.log2.negx"
hist_H_L_Ratio_normalized_2R_log2_negx <- ggplot(proteinGroups_strict_filtered_logtransformed_negx, aes_string(x = H_L_Ratio_normalized_2R_log2_negx)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("-"*"Log"[2]*"(Ratio H/L normalized 2R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_2R_log2_negx)

H_L_Ratio_normalized_3R_log2_negx <- "Ratio.H.L.normalized.3R.log2.negx"
hist_H_L_Ratio_normalized_3R_log2_negx <- ggplot(proteinGroups_strict_filtered_logtransformed_negx, aes_string(x = H_L_Ratio_normalized_3R_log2_negx)) +
  geom_histogram(fill = "#2980B9", color = "black", bins = 100) +  # Adjust the fill colour and outline colour
  xlim(-2.5, 2.5) +
  labs(x = expression("-"*"Log"[2]*"(Ratio H/L normalized 3R)"),
       y = expression("Counts"))
print(hist_H_L_Ratio_normalized_3R_log2_negx)

# Arrange plots in a grid
combined_neg_log2_normalized_histograms <- grid.arrange(hist_H_L_Ratio_normalized_1F_log2, hist_H_L_Ratio_normalized_1R_log2_negx, hist_H_L_Ratio_normalized_2F_log2,
                                                        hist_H_L_Ratio_normalized_2R_log2_negx, hist_H_L_Ratio_normalized_3F_log2, hist_H_L_Ratio_normalized_3R_log2_negx,
                                                        ncol = 2, nrow = 3)

# Save the combined plot as a PNG file
ggsave("combined_neg_log2_normalized_histograms.png", plot = combined_neg_log2_normalized_histograms, width = 10, height = 8)

#############################################################################################################

#8. Statistical tests
#Compute the t-test for the normalized log2 ratios:
H_L_Ratio_normalized_after_transformations = c("Ratio.H.L.normalized.1F.log2", "Ratio.H.L.normalized.1R.log2.negx", "Ratio.H.L.normalized.2F.log2", "Ratio.H.L.normalized.2R.log2.negx", "Ratio.H.L.normalized.3F.log2", "Ratio.H.L.normalized.3R.log2.negx")

p_values <- apply(proteinGroups_strict_filtered_logtransformed_negx[, H_L_Ratio_normalized_after_transformations], 1, function(x) {
  # we have to wrap the test into a tryCatch, since the t-test is not working when there are too many NA's
  tryCatch({
    one_t_test <- t.test(x, mu = 0, var.equal = FALSE)
    one_t_test$p.value
  }, error = function(cond){
    NA
  })
})

#number of significant p-values
sum(p_values <= 0.05, na.rm = TRUE)
#But we have to adjust for multiple testing (e.g. with Benjamini-Hochberg):

p_values_adj <- p.adjust(p_values, method = "BH")
#How many proteins pass t-test with p-value < 0.05 now?
sum(p_values_adj <= 0.05, na.rm = TRUE)

protein_names <- proteinGroups_strict_filtered_logtransformed_negx$Protein.names
gene_names <- proteinGroups_strict_filtered_logtransformed_negx$Gene.names
log2_FC_1F = proteinGroups_strict_filtered_logtransformed_negx$Ratio.H.L.normalized.1F.log2
log2_FC_1R = proteinGroups_strict_filtered_logtransformed_negx$Ratio.H.L.normalized.1R.log2.negx
log2_FC_2F = proteinGroups_strict_filtered_logtransformed_negx$Ratio.H.L.normalized.2F.log2
log2_FC_2R = proteinGroups_strict_filtered_logtransformed_negx$Ratio.H.L.normalized.2R.log2.negx
log2_FC_3F = proteinGroups_strict_filtered_logtransformed_negx$Ratio.H.L.normalized.3F.log2
log2_FC_3R = proteinGroups_strict_filtered_logtransformed_negx$Ratio.H.L.normalized.3R.log2.negx

#Create a data frame with protein names, gene names, p-values, and adjusted p-values
df_p_values_adj <- data.frame(Protein = protein_names, Gene = gene_names,
                              log2_FC_1F = log2_FC_1F, log2_FC_1R = log2_FC_1R,
                              log2_FC_2F = log2_FC_2F, log2_FC_2R = log2_FC_2R,
                              log2_FC_3F = log2_FC_3F, log2_FC_3R = log2_FC_3R,
                              p_value = p_values,
                              p_value_adj = p_values_adj)

# Calculate the average of log2_FC values for each row
df_p_values_adj$avg_log2_FC <- rowMeans(df_p_values_adj[, c("log2_FC_1F", "log2_FC_1R", "log2_FC_2F", "log2_FC_2R", "log2_FC_3F", "log2_FC_3R")], na.rm = TRUE)

# Reorder the columns to place avg_log2_FC before p_value_adj
df_p_values_adj <- df_p_values_adj[, c("Protein", "Gene", "log2_FC_1F", "log2_FC_1R", "log2_FC_2F", "log2_FC_2R", "log2_FC_3F", "log2_FC_3R", "avg_log2_FC", "p_value", "p_value_adj")]
df_p_values_adj

significant_proteins <- subset(df_p_values_adj, p_values_adj < 0.05)

# Sort based on descending p_value_adj
significant_proteins_sorted_p_value_adj <- significant_proteins[order(significant_proteins$p_value_adj), ]
significant_proteins_sorted_p_value_adj

# Sort based on descending average of all log2 fold changes
significant_proteins_sorted_log2_FC_desc <- significant_proteins[order(-significant_proteins$avg_log2_FC), ]
significant_proteins_sorted_log2_FC_desc

# Sort based on ascending average of all log2 fold changes
significant_proteins_sorted_log2_FC_asc <- significant_proteins[order(significant_proteins$avg_log2_FC), ]
significant_proteins_sorted_log2_FC_asc

#7.1 Extracting significant hits after Benjamini-Hochberg correction
install.packages("xlsx")
library("xlsx")

# Define the path and filename for the Excel file
significant_proteins_excel <- "significant_proteins.xlsx"
# Write the significant proteins to the Excel file
write.xlsx(significant_proteins, file = significant_proteins_excel, row.names = FALSE)

#7.2 Volcano plot of adjusted p-values
install.packages("ggrepel")
library("ggrepel")
# add a column of NO
df_p_values_adj$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_p_values_adj$diffexpressed[df_p_values_adj$avg_log2_FC > 0.6 & df_p_values_adj$p_value_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_p_values_adj$diffexpressed[df_p_values_adj$avg_log2_FC < -0.6 & df_p_values_adj$p_value_adj < 0.05] <- "DOWN"

#Create a new column "diffexpressedlabel" to the dataframe, that will contain the name of genes differentially expressed (NA in case they are not)
df_p_values_adj$diffexpressedlabel <- NA
df_p_values_adj$diffexpressedlabel[df_p_values_adj$diffexpressed != "NO"] <- df_p_values_adj$Gene[df_p_values_adj$diffexpressed != "NO"]

#Extract the downregulated and upregulated genes in an Excel file
# Filter dataframe to select rows where diffexpressedlabel is not NA
significant_genes <- df_p_values_adj[!is.na(df_p_values_adj$diffexpressedlabel), c("Gene", "Protein", "log2_FC_1F", "log2_FC_1R", "log2_FC_2F", "log2_FC_2R", "log2_FC_3F", "log2_FC_3R", "avg_log2_FC", "p_value_adj", "diffexpressed")]

# Write to Excel file
write.xlsx(significant_genes, "significant_genes.xlsx", row.names = FALSE)

# Create a basic volcano plot with threshold lines
volcanoplotproteinGroups <- ggplot(data = df_p_values_adj, aes(x = avg_log2_FC, y = -log10(p_value_adj), col = diffexpressed, label = diffexpressedlabel)) +
  theme_minimal() +  # Apply a minimal theme
  geom_text_repel() +  # Add text labels with repulsion to avoid overlap
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1, alpha = ifelse(df_p_values_adj$diffexpressed == "NO", 0.2, 1)) + 
  scale_color_manual(values = c("#CC0000", "black", "#00AFBB"), # to set the colours of the variables
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe
  labs(color = "Protein expression", #legend title
       x = expression("Log"[2] * "(fold change)"), y = expression("-Log"[10] * "(p-value)"),
       title = expression("Protein expression in" ~ italic("proteinGroups") ~ "KO vs. WT HEK293T cells")) +    
  theme(plot.title = element_text(hjust = 0.5)) + #Centre the title
  guides(col = guide_legend(override.aes = list(size = 1.5))) # Adjust dot size in legend

# Save the volcano plot as a PDF file
ggsave("volcanoplotproteinGroups.pdf",  # Specify the filename for the PDF
       plot = volcanoplotproteinGroups, width = 8, height = 6)  # Specify the ggplot object to be saved

# Save the volcano plot as a PNG file
ggsave("volcanoplotproteinGroups.png",  # Specify the filename for the PDF
       plot = volcanoplotproteinGroups, width = 8, height = 6)  # Specify the ggplot object to be saved

# You can also save it as png, jpeg, tiff, bmp, svg, ps...
# For more on saving plots in R check https://biocorecrg.github.io/CRG_RIntroduction/with-the-console.html

#7.3 Volcano plot of uncorrected p-values
df_p_values_raw = subset(df_p_values_adj, select = -p_value_adj)
df_p_values_raw

# add a column of NO
df_p_values_raw$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_p_values_raw$diffexpressed[df_p_values_raw$avg_log2_FC > 0.6 & df_p_values_raw$p_value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_p_values_raw$diffexpressed[df_p_values_raw$avg_log2_FC < -0.6 & df_p_values_raw$p_value < 0.05] <- "DOWN"

#Create a new column "diffexpressedlabel" to the dataframe, that will contain the name of genes differentially expressed (NA in case they are not)
df_p_values_raw$diffexpressedlabel <- NA
df_p_values_raw$diffexpressedlabel[df_p_values_raw$diffexpressed != "NO"] <- df_p_values_raw$Gene[df_p_values_raw$diffexpressed != "NO"]

#Extract the downregulated and upregulated genes in an Excel file
# Filter dataframe to select rows where diffexpressedlabel is not NA
significant_genes_raw_p_values <- df_p_values_raw[!is.na(df_p_values_raw$diffexpressedlabel), c("Gene", "Protein", "log2_FC_1F", "log2_FC_1R", "log2_FC_2F", "log2_FC_2R", "log2_FC_3F", "log2_FC_3R", "avg_log2_FC", "p_value", "diffexpressed")]

# Write to Excel file
write.xlsx(significant_genes_raw_p_values, "significant_genes_raw_p_values.xlsx", row.names = FALSE)

# Create a basic volcano plot with threshold lines
volcanoplotproteinGroups_raw_p_values <- ggplot(data = df_p_values_raw, aes(x = avg_log2_FC, y = -log10(p_value), col = diffexpressed, label = diffexpressedlabel)) +
  theme_minimal() +  # Apply a minimal theme
  geom_text_repel() +  # Add text labels with repulsion to avoid overlap
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1, alpha = ifelse(df_p_values_raw$diffexpressed == "NO", 0.1, 1)) + 
  scale_color_manual(values = c("#CC0000", "black", "#00AFBB"), # to set the colours of the variables
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe
  labs(color = "Protein expression", #legend title
       x = expression("Log"[2] * "(fold change)"), y = expression("-Log"[10] * "(p-value)"),
       title = expression("Protein expression in" ~ italic("proteinGroups") ~ "KO vs. WT HEK293T cells")) +    
  theme(plot.title = element_text(hjust = 0.5)) + #Centre the title
  guides(col = guide_legend(override.aes = list(size = 1.5))) # Adjust dot size in legend

# Save the volcano plot as a PDF file
ggsave("volcanoplotproteinGroups_raw_p_values.pdf",  # Specify the filename for the PDF
       plot = volcanoplotproteinGroups_raw_p_values, width = 8, height = 6)  # Specify the ggplot object to be saved

# Save the volcano plot as a PNG file
ggsave("volcanoplotproteinGroups_raw_p_values.png",  # Specify the filename for the PDF
       plot = volcanoplotproteinGroups_raw_p_values, width = 8, height = 6)  # Specify the ggplot object to be saved
# You can also save it as png, jpeg, tiff, bmp, svg, ps...
# For more on saving plots in R check https://biocorecrg.github.io/CRG_RIntroduction/with-the-console.html
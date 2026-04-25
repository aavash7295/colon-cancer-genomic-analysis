## BINF702 FINAL PROJECT
## Colon Cancer Genomic Analysis | NCI60 Dataset


# ---------------- LIBRARIES ----------------

library(ISLR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(Rtsne)
library(factoextra)
library(gridExtra)
library(mgcv)
library(kohonen)
library(randomForest)   # FIX: removed duplicate library(randomForest) line
library(nnet)
library(NeuralNetTools)

# ---------------- LOAD DATA ----------------
data("NCI60")
nci.data <- NCI60$data
nci.labs <- NCI60$labs


# ---------------- FILTER COLON ----------------
colon_indices    <- which(nci.labs == "COLON")
colon_expression <- nci.data[colon_indices, ]

# ---------------- Z-SCORE SCALING ----------------
scaled_colon <- scale(colon_expression)

# ---------------- Z-RANGE COUNTS ----------------
count_1_to_2 <- colSums(scaled_colon > 1 & scaled_colon <= 2, na.rm = TRUE)

# ---------------- MARKET MARKERS (Z in 1-2 for >2 patients) ----------------
high_consistency_list <- colnames(scaled_colon)[which(count_1_to_2 > 2)]
market_numeric        <- scaled_colon[, high_consistency_list, drop = FALSE]

cat("Market Genes found:", length(high_consistency_list), "\n")
print(high_consistency_list)
cat("Market matrix dimensions:", dim(market_numeric), "\n")

# ---------------- BACKGROUND OUTLIERS (|Z| > 2) ----------------
mask_sig         <- abs(scaled_colon) > 2
sig_gene_indices <- which(colSums(mask_sig, na.rm = TRUE) > 0)
background_numeric <- scaled_colon[, sig_gene_indices, drop = FALSE]

cat("Background outlier genes:", ncol(background_numeric), "\n")

# ---------------- OUTLIER SUMMARY ----------------
summary_counts <- data.frame(
  Patient           = rownames(scaled_colon),
  Upregulated       = rowSums(scaled_colon >  2, na.rm = TRUE),
  Downregulated     = rowSums(scaled_colon < -2, na.rm = TRUE),
  Total_Significant = rowSums(abs(scaled_colon) > 2, na.rm = TRUE)
)
print(summary_counts)

# ---------------- STACKED BAR CHART ----------------
summary_long <- summary_counts %>%
  select(Patient, Upregulated, Downregulated) %>%
  pivot_longer(-Patient, names_to = "Direction", values_to = "Count")

ggplot(summary_long, aes(x = Patient, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Upregulated" = "#B22222", "Downregulated" = "#008080")) +
  labs(title = "Clinical Genomic Outlier Identification",
       subtitle = "Genes with |Z| > 2.0 per Patient (Cohort V42-V48)",
       x = "Patient", y = "Gene Count") +
  theme_minimal() + theme(legend.position = "bottom")



# ---------------- HEATMAP ----------------
pheatmap(market_numeric,
         main  = "Market Marker Heatmap (High-Consistency Driver Genes)",
         color = colorRampPalette(c("navy", "white", "red"))(50))

# ---------------- PCA + t-SNE ----------------
pca_res <- prcomp(market_numeric, scale. = TRUE)
pca_df  <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                      Patient = rownames(market_numeric))
p_pca   <- ggplot(pca_df, aes(PC1, PC2, color = Patient, label = Patient)) +
  geom_point(size = 4) + geom_text(vjust = -1) +
  labs(title = "PCA - Market Markers") + theme_minimal()

set.seed(42)
tsne_res <- Rtsne(market_numeric, perplexity = 2, check_duplicates = FALSE)
tsne_df  <- data.frame(X = tsne_res$Y[,1], Y = tsne_res$Y[,2],
                       Patient = rownames(market_numeric))
p_tsne   <- ggplot(tsne_df, aes(X, Y, color = Patient, label = Patient)) +
  geom_point(size = 4) + geom_text(vjust = -1) +
  labs(title = "t-SNE - Market Markers") + theme_minimal()

grid.arrange(p_pca, p_tsne, ncol = 2)


# ---------------- K-MEANS (k=3) ----------------
set.seed(123)
km_market <- kmeans(market_numeric, centers = 3, nstart = 25)
fviz_cluster(km_market, data = market_numeric,
             main = "K-Means Clustering (k=3) - Market Genes")

# ---------------- HIERARCHICAL CLUSTERING (Ward's D2) ----------------
dist_matrix <- dist(market_numeric, method = "euclidean")
hc_model    <- hclust(dist_matrix, method = "ward.D2")
fviz_dend(hc_model, k = 3, rect = TRUE,
          main = "Hierarchical Clustering - Ward's D2",
          sub  = "Euclidean Distance | k=3 groups")

# Concordance: do K-means and HC agree on groupings?
hc_clusters  <- cutree(hc_model, k = 3)
concordance  <- data.frame(Patient      = rownames(market_numeric),
                           KMeans_Clust = km_market$cluster,
                           HC_Clust     = hc_clusters)
cat("\nClustering Concordance (K-Means vs Hierarchical):\n")
print(concordance)
cat("Agreement rate:",
    round(mean(km_market$cluster == hc_clusters) * 100, 1), "%\n")

# ---------------- SELF-ORGANIZING MAP (SOM) ----------------

library(kohonen)

# Define SOM grid (same as your presentation: 2x3 hexagonal)
som_grid <- somgrid(xdim = 2, ydim = 3, topo = "hexagonal")

set.seed(123)

# Train SOM on Market genes (driver genes)
som_model <- som(as.matrix(market_numeric),
                 grid = som_grid,
                 rlen = 500,
                 alpha = c(0.05, 0.01))

# ---- Plot 1: Patient Mapping ----
plot(som_model,
     type = "mapping",
     labels = rownames(market_numeric),
     main = "SOM: Patient Topological Mapping",
     col = "darkblue",
     pch = 19)

# ---- Plot 2: U-Matrix (Cluster Boundaries) ----
plot(som_model,
     type = "dist.neighbours",
     main = "SOM: U-Matrix (Cluster Boundaries)",
     palette.name = grey.colors)

# ---- Plot 3: Codebook Vectors (Gene Influence) ----
plot(som_model,
     type = "codes",
     main = "SOM: Codebook Vectors (Gene Influence)")


# ---------------- GAM (k=3 spline, N=7) ----------------
# FIX: guard already present; use mean of ALL market genes (not just first)
if (ncol(background_numeric) >= 10) {
  gam_data  <- data.frame(
    Market_Mean    = rowMeans(market_numeric),
    Background_Avg = rowMeans(background_numeric[, 1:10])
  )
  gam_model <- gam(Background_Avg ~ s(Market_Mean, k = 3), data = gam_data)
  cat("\nGAM Summary:\n")
  print(summary(gam_model))

  ggplot(gam_data, aes(x = Market_Mean, y = Background_Avg)) +
    geom_point(size = 4, color = "red") +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "blue") +
    labs(title    = "GAM: Market Drivers Predicting Background Outliers",
         subtitle = paste("Deviance explained:",
                          round(summary(gam_model)$dev.expl * 100, 1), "%"),
         x = "Mean Market Gene Z-score", y = "Mean Background Z-score") +
    theme_minimal()
}


# ---------------- RANDOM FOREST ----------------
# Target: High vs Low genomic instability (above/below median outlier count)
colnames(market_numeric) <- make.names(colnames(market_numeric))
instability <- ifelse(summary_counts$Total_Significant >
                        median(summary_counts$Total_Significant), "High", "Low")

rf_data       <- as.data.frame(market_numeric)
rf_data$Class <- as.factor(instability)

set.seed(123)
rf_model <- randomForest(Class ~ ., data = rf_data,
                         ntree = 500, importance = TRUE)
cat("\nRandom Forest Model:\n")
print(rf_model)

# Gini Importance — ranked gene table
importance_df <- data.frame(
  Gene            = rownames(importance(rf_model)),
  Gini_Importance = importance(rf_model, type = 2)[, 1]
) %>% arrange(desc(Gini_Importance))

cat("\nTop Genes by Gini Importance:\n")
print(importance_df)
varImpPlot(rf_model, main = "Random Forest: Variable Importance (Gini)")


## NEURAL NETWORK


# Normalize function
normalize <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-10)

# Create data FIRST
nn_data <- as.data.frame(apply(market_numeric, 2, normalize))

# THEN fix column names
colnames(nn_data) <- make.names(colnames(nn_data))

# Add class
nn_data$Class <- as.factor(instability)

# Train model
set.seed(123)
nn_model <- nnet(Class ~ ., data = nn_data,
                 size = 3, decay = 0.01,
                 maxit = 500, trace = FALSE)

print(nn_model)

# Predictions
nn_pred <- predict(nn_model, nn_data, type = "class")
nn_acc  <- mean(nn_pred == nn_data$Class) * 100

cat(sprintf("Neural Network Training Accuracy: %.1f%%\n", nn_acc))

# Olden's Influence Scores — ranks genes by contribution to classification

## OLDEN (FIXED PROPERLY)


olden_scores <- olden(nn_model, bar_plot = FALSE)

# Convert safely depending on structure
if (is.matrix(olden_scores) || is.data.frame(olden_scores)) {
  olden_df <- data.frame(
    Gene = rownames(olden_scores),
    Olden_Importance = olden_scores[,1]
  )
} else {
  olden_df <- data.frame(
    Gene = names(olden_scores),
    Olden_Importance = unlist(olden_scores)
  )
}

olden_df <- olden_df %>%
  arrange(desc(abs(Olden_Importance)))

print(olden_df)

# ---------------- CROSS-MODEL BIOMARKER VALIDATION ----------------
# Genes ranked highly by BOTH models = high-confidence biomarkers
importance_df$RF_Rank <- rank(-importance_df$Gini_Importance)
olden_df$NN_Rank      <- rank(-abs(olden_df$Olden_Importance))

combined <- merge(
  importance_df[, c("Gene", "Gini_Importance", "RF_Rank")],
  olden_df[,     c("Gene", "Olden_Importance", "NN_Rank")],
  by = "Gene"
) %>% mutate(Mean_Rank = (RF_Rank + NN_Rank) / 2) %>%
  arrange(Mean_Rank)

cat("\n=== HIGH-CONFIDENCE BIOMARKERS (confirmed by both models) ===\n")
print(head(combined, 5))

# Side-by-side importance plot for top confirmed biomarkers
top5 <- head(combined, 5) %>%
  select(Gene, Gini_Importance, Olden_Importance) %>%
  mutate(Olden_Importance = abs(Olden_Importance)) %>%
  pivot_longer(-Gene, names_to = "Model", values_to = "Score") %>%
  mutate(Model = recode(Model,
                        "Gini_Importance"  = "RF (Gini)",
                        "Olden_Importance" = "NN (Olden)"))

ggplot(top5, aes(x = reorder(Gene, Score), y = Score, fill = Model)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("RF (Gini)" = "#2E9FDF", "NN (Olden)" = "#E7B800")) +
  coord_flip() +
  labs(title    = "Cross-Model Biomarker Validation",
       subtitle = "Top genes confirmed as drivers by both RF and Neural Network",
       x = "Gene", y = "Importance Score", fill = "Model") +
  theme_minimal()



# Load required libraries
library(jsonlite)
library(ggplot2)
library(scales)

# Load data from JSON
json_file <- "data/ar53_msa_reps_r220_squared_s10_pruned60_score_dict.json"
data <- fromJSON(json_file)

# Convert list to data frame for ggplot
long_data <- do.call(rbind, lapply(names(data), function(name) {
  data.frame(label = name, chiscore = unlist(data[[name]]))
}))

label_colors <- c(
  "before_real" = "#ffbf4d", #Observed
  "after_real" = "#4d4dff", #chi-square after pruning
  "before_permuted" = "#e6337c" #Permuted
)

pdf("figures/arGTDB_chi_distributions.pdf", width = 5, height = 2.5)
# Plot KDE (density) with log scale
ggplot(long_data, aes(x = chiscore, fill = label)) +
  geom_density(alpha = 0.8, color = "black") +
  scale_x_log10() +
  scale_fill_manual(values = label_colors) +
  labs(title = "Density Plot of Chi-Square Scores",
       x = expression(chi^2 ~ "score"),
       y = "Density",
       fill = "Distributions",
       color = "Distributions") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.background = element_blank(),       
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(),              
    axis.text = element_text(),               
    axis.title = element_text(),              
    legend.position = "right"
  )
 
dev.off()


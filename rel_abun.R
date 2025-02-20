library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
abundance_data <- read.csv("zotu_table_tax.txt", header = TRUE, sep = "\t")

metadata <- read.csv("metadata.txt", header = TRUE, sep = "\t")

# Check the column names
colnames(abundance_data) <- make.names(colnames(abundance_data))

# Convert the data to long format for easier manipulation
abundance_long <- abundance_data %>%
  pivot_longer(cols = -c(X.OTU_ID, Taxonomy), names_to = "Sample", values_to = "Abundance")

# Merge with metadata to include species information
abundance_long <- abundance_long %>%
  left_join(metadata, by = "Sample")

# Calculate relative abundance
abundance_long <- abundance_long %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance) * 100)

# Reorder the Sample factor levels based on Species and then Sample
abundance_long$Sample <- factor(abundance_long$Sample, levels = abundance_long$Sample[order(abundance_long$Species, abundance_long$Sample)])

# Reorder the Sample factor levels based on Species
abundance_long$Sample <- factor(abundance_long$Sample, levels = rev(sort(unique(abundance_long$Sample[order(abundance_long$Species)])))

# Plot the data
ggplot(abundance_long, aes(x = Sample, y = RelativeAbundance, fill = Taxonomy)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Microbial composition  - Myrmica",
       x = "Sample",
       y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1) ,axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12)) + 
  scale_fill_manual(values=c("#8e063b", "darkgreen",  "lightgreen","#ef9708","lightblue")) +
  coord_flip()

####################


library(tidyverse)
library(dplyr)
library(data.table)
library(ggpubr)
library(ggforce)

setwd("/Users/liboxuan/Desktop/Browns_lab/Clade_B_analysis/vfdb")
# Read the CSV file into a data frame
vfdb <- fread("vfdb_labeled_VF_Gene.csv", header = TRUE)
host_trait <- fread ("trait_host.csv", header = TRUE)
vfdb <- host_trait[vfdb, on = "Genome_ID"]
# Modify the vfdb data frame fill the missing "Early CF Expanded"
# vfdb <- vfdb %>%
#   mutate(Genome_ID = as.character(Genome_ID),  # Ensure column 6 is character type
#          Niche = ifelse(grepl("^GCA_02510|^GCA_02512|^GCA_02513", Genome_ID), "Early CF Expanded", Niche))

# Process ---------------------------
processed_data <- vfdb %>%
  group_by(Genome_ID, Host_Related) %>%
  count(VF) %>%
  ungroup() %>%
  group_by(Host_Related, VF) %>%
  summarize(avg_vfdb_perGenome = mean(n))

filtered_data <- processed_data %>%
  filter(Host_Related %in% c("1", "2"))

# Plotting
ggplot(filtered_data, aes(x = as.factor(Host_Related), y = avg_vfdb_perGenome, fill = VF)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(breaks = c(0, 1, 5, 10, 20, 30, 40, 50)) +
  ylab("Average vfdb per Genome") +
  xlab("Host_Related") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1))

# Convert to Pres/Abs ---------------------------

vfdb_dat_PresAbs <- vfdb %>% 
  select(Genome_ID, VF, Host_Related) %>%
  group_by(Genome_ID, VF, Host_Related) %>%
  count() %>%
  ungroup()

# Filtering by Clade
vfdb_dat_PresAbs <- vfdb_dat_PresAbs %>%
  filter(Host_Related %in% c("1", "0"))

vfdb_simple <- vfdb %>%
  select(VF) %>%
  distinct() %>%
  group_by(`VF`) %>%
  slice(1) %>%
  ungroup()

vfdb_dat_PresAbs_plot <- vfdb_dat_PresAbs %>% left_join(vfdb_simple)

vfdb_dat_PresAbs_plot <- vfdb_dat_PresAbs_plot %>%
  mutate(log_n = log10(n)) 


### normal plot

vfdb_dat_PresAbs_plot %>%
  rename(VF = `VF`) %>%
  ggplot(aes(VF, Genome_ID, fill = log_n)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 0.5), 
        axis.text.y = element_text(size = 0.5)) +
  facet_col(~ Host_Related, scales = "free_y", space = "free",
            strip.position = "top") +
  labs(fill = "Log Count")

ggsave("plot.pdf", width = 20, height = 16, units = "in")





# top vf
library(readr)
library(data.table)
library(dplyr)
setwd("/Users/liboxuan/Desktop/Browns_lab/Clade_B_analysis/vfdb")
# Read the CSV file into a data frame
vfdb <- fread("vfdb_labeled_VF_Gene.csv", header = TRUE)
# Get total number of unique genomes
total_genomes <- n_distinct(vfdb$Genome_ID)

# Count in how many genomes each Gene appears
gene_dist <- vfdb %>%
  group_by(Accession) %>%
  summarise(genome_count = n_distinct(Genome_ID)) %>%
  mutate(freq = genome_count / total_genomes * 100) %>%
  filter(freq >= 95) %>%
  arrange(desc(freq))

# View top widely distributed genes
print(gene_dist)

write_csv(gene_dist, "/Users/liboxuan/Desktop/Browns_lab/Clade_B_analysis/vfdb/core_VF_host_related.csv")

# Generate anotehr df with list of accession that kind of evenly distributed between clade 1 and clade 2. Dont need to consider the frequency

# Count number of genomes per clade
genomes_per_clade <- vfdb %>%
  filter(Clade %in% c("1", "2")) %>%
  distinct(Genome_ID, Clade) %>%
  count(Clade, name = "total_genomes")

# Count how many genomes each VF is found in, per clade
vf_freq_by_clade <- vfdb %>%
  filter(Clade %in% c("1", "2")) %>%
  group_by(Accession, Clade) %>%
  summarise(genome_count = n_distinct(Genome_ID), .groups = "drop") %>%
  left_join(genomes_per_clade, by = "Clade") %>%
  mutate(freq = genome_count / total_genomes)

write_csv(vf_freq_by_clade, "/Users/liboxuan/Desktop/Browns_lab/Clade_B_analysis/vfdb/raw_vfdb_frequency.csv")

# Reshape to wide format: Clade_1_freq, Clade_2_freq
vf_freq_wide <- vf_freq_by_clade %>%
  select(Accession, Clade, freq) %>%
  pivot_wider(names_from = Clade, values_from = freq, names_prefix = "Clade_", values_fill = 0)

# Filter VFs that are relatively evenly distributed between clades
evenly_distributed_vfs <- vf_freq_wide %>%
  mutate(
    min_freq = pmin(Clade_1, Clade_2),
    max_freq = pmax(Clade_1, Clade_2),
    freq_ratio = ifelse(max_freq == 0, NA, min_freq / max_freq)
  ) %>%
  filter(freq_ratio >= 0.5) %>%         # Adjust threshold as needed
  arrange(desc(freq_ratio))    # Adjust this threshold as needed

# Preview and save
print(evenly_distributed_vfs)
write_csv(evenly_distributed_vfs, "/Users/liboxuan/Desktop/Browns_lab/Clade_B_analysis/vfdb/even_VF.csv")
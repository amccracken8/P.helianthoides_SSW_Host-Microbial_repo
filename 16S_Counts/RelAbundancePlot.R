library(tidyverse)
library(data.table)
library(foreach)


microbes <- fread("level-7.csv")
meta <- fread("pycno_samples_NE.txt")


df <- microbes[microbes$sampleID %in% meta$sample,]

df <- df %>% dplyr::select(sampleID,colnames(df)[grepl("^(d_)", colnames(df))])
df <- tibble::as_tibble(df)

df_long <- df %>%
  pivot_longer(
    cols = -sampleID,
    names_to = "microbe",
    values_to = "count"
  )

df_long <- df_long %>%
  group_by(sampleID) %>%
  mutate(rel_abundance = count / sum(count)) %>%
  ungroup()


df_long %>%
  group_by(microbe) %>%
  mutate(total_abundance = sum(rel_abundance)) %>%
  ungroup() %>%
  mutate(microbe = forcats::fct_reorder(microbe, total_abundance))


## top abundant microbes
df_long_top <- df_long %>%
  group_by(microbe) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(microbe = fct_lump_n(microbe, n = 15, w = total)) %>%
  group_by(microbe) %>%
  mutate(total_lumped = sum(count)) %>%
  ungroup() %>%
  mutate(microbe = fct_reorder(microbe, total_lumped))

df_long_top <- df_long_top %>%
  group_by(sampleID, microbe) %>%
  summarise(
    count = sum(count),
    rel_abundance = sum(rel_abundance),
    .groups = "drop"
  )


df_long_top$microbe <- as.character(df_long_top$microbe)

df_long_top$microbe <- ifelse(
  grepl("f__[^;]+", df_long_top$microbe),
  
  # If family exists → keep from f__ onward, drop trailing empty s__
  sub(".*?(f__[^;]+(?:;g__[^;]+)?).*", "\\1", df_long_top$microbe),
  
  # If no family → extract class
  ifelse(
    grepl("c__[^;]+", df_long_top$microbe),
    sub(".*?(c__[^;]+).*", "\\1", df_long_top$microbe),
    df_long_top$microbe
  )
)


# Calculate total abundance per microbe
abundance_order <- df_long_top %>%
  group_by(microbe) %>%
  summarise(total = sum(rel_abundance, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%          # most abundant first
  pull(microbe)

# Relevel factor
df_long_top$microbe <- factor(df_long_top$microbe,
                              levels = abundance_order)

r <- ggplot(df_long_top, aes(x = sampleID, y = rel_abundance, fill = microbe)) +
  geom_col(
    width = 0.9,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = rev(c(
    "#000000", "#E69F00", "#56B4E9", "firebrick",
    "#F0E442", "#0072B2", "pink2", "#CC79A7",
    "seagreen", "firebrick2", "violet", "skyblue3",
    "salmon", "grey", "darkred", "cornflowerblue"
  ))) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Relative abundance"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.text.y = element_text(angle = 90, vjust = 0.5,face="bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.6, "cm")
  )

## absoute abundance

a <- ggplot(df_long_top, aes(x = sampleID, y = count, fill = microbe)) +
  geom_col(
    width = 0.9,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = rev(c(
    "#000000", "#E69F00", "#56B4E9", "firebrick",
    "#F0E442", "#0072B2", "pink2", "#CC79A7",
    "seagreen", "firebrick2", "violet", "skyblue3",
    "salmon", "grey", "darkred", "cornflowerblue"
  ))) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Absolute abundance (counts)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.text.y = element_text(angle = 90, vjust = 0.5,face="bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.6, "cm")
  )


library(patchwork)

(a + r) +
  plot_layout(guides = "collect")

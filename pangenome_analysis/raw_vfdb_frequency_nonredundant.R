library(tidyverse)


# ============================================================
# 1. Read frequency table
# ============================================================

vf_freq <- read_csv(
  "raw_vfdb_frequency_with_locus_tag.csv",
  show_col_types = FALSE
) %>%
  mutate(
    Clade = recode(
      as.character(Clade),
      "1" = "A",
      "2" = "B"
    ),
    
    locus_tag = na_if(
      trimws(as.character(locus_tag)),
      ""
    )
  )


# ============================================================
# 2. Calculate total coverage for every VFG
# ============================================================

vfg_coverage <- vf_freq %>%
  group_by(
    Accession,
    locus_tag
  ) %>%
  summarise(
    combined_genome_count = sum(
      genome_count,
      na.rm = TRUE
    ),
    
    number_of_clades_detected = n_distinct(
      Clade
    ),
    
    .groups = "drop"
  )


# ============================================================
# 3. Choose one representative per mapped locus
#
# Selection rule:
#   1. Highest combined genome count
#   2. Accession name as a deterministic tie-breaker
#
# Do not use p-value or frequency difference for selection.
# ============================================================

mapped_representatives <- vfg_coverage %>%
  filter(
    !is.na(locus_tag)
  ) %>%
  arrange(
    locus_tag,
    desc(combined_genome_count),
    desc(number_of_clades_detected),
    Accession
  ) %>%
  group_by(
    locus_tag
  ) %>%
  mutate(
    n_vfg_for_locus = n(),
    representative_rank = row_number()
  ) %>%
  filter(
    representative_rank == 1
  ) %>%
  ungroup()


# ============================================================
# 4. Retain unmapped VFGs individually
# ============================================================

unmapped_representatives <- vfg_coverage %>%
  filter(
    is.na(locus_tag)
  ) %>%
  mutate(
    n_vfg_for_locus = NA_integer_,
    representative_rank = 1L
  )


# ============================================================
# 5. Create representative list
# ============================================================

representative_list <- bind_rows(
  mapped_representatives,
  unmapped_representatives
) %>%
  select(
    Accession,
    locus_tag,
    combined_genome_count,
    number_of_clades_detected,
    n_vfg_for_locus
  )


# ============================================================
# 6. Filter original frequency table
# ============================================================

vf_freq_nonredundant <- vf_freq %>%
  semi_join(
    representative_list %>%
      select(Accession),
    by = "Accession"
  ) %>%
  left_join(
    representative_list %>%
      select(
        Accession,
        combined_genome_count,
        n_vfg_for_locus
      ),
    by = "Accession"
  ) %>%
  arrange(
    locus_tag,
    Accession,
    Clade
  )


# ============================================================
# 7. Record discarded repeated entries
# ============================================================

discarded_vfgs <- vfg_coverage %>%
  filter(
    !is.na(locus_tag)
  ) %>%
  anti_join(
    mapped_representatives %>%
      select(
        Accession,
        locus_tag
      ),
    by = c(
      "Accession",
      "locus_tag"
    )
  ) %>%
  left_join(
    mapped_representatives %>%
      select(
        locus_tag,
        representative_Accession = Accession,
        representative_coverage =
          combined_genome_count
      ),
    by = "locus_tag"
  ) %>%
  arrange(
    locus_tag,
    desc(combined_genome_count)
  )


# ============================================================
# 8. Save outputs
# ============================================================

write_csv(
  representative_list,
  "VF_representative_list_one_per_locus.csv",
  na = ""
)

write_csv(
  vf_freq_nonredundant,
  "raw_vfdb_frequency_nonredundant.csv",
  na = ""
)

write_csv(
  discarded_vfgs,
  "VF_discarded_repeated_locus_entries.csv",
  na = ""
)


# ============================================================
# 9. Summary
# ============================================================

cat(
  "Original unique VFGs:",
  n_distinct(vf_freq$Accession),
  "\n"
)

cat(
  "Mapped unique locus tags:",
  n_distinct(
    vf_freq$locus_tag,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Mapped representatives retained:",
  nrow(mapped_representatives),
  "\n"
)

cat(
  "Unmapped VFGs retained individually:",
  nrow(unmapped_representatives),
  "\n"
)

cat(
  "Final nonredundant units:",
  nrow(representative_list),
  "\n"
)

cat(
  "Repeated VFG entries discarded:",
  nrow(discarded_vfgs),
  "\n"
)
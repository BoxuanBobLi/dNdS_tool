# ============================================================
# Test VF presence/absence enrichment between Clade A and B
# Input: raw_vfdb_frequency_with_locus_tag.csv
# ============================================================

library(tidyverse)


# ------------------------------------------------------------
# 0. Read input and rename clades
# ------------------------------------------------------------

vf_freq_by_clade <- read_csv(
  "raw_vfdb_frequency_nonredundant.csv",
  show_col_types = FALSE
) %>%
  mutate(
    Clade = as.character(Clade),
    
    # Internal values 1 and 2 become biological labels A and B
    Clade = recode(
      Clade,
      "1" = "A",
      "2" = "B"
    ),
    
    locus_tag = na_if(
      trimws(as.character(locus_tag)),
      ""
    )
  )


# Check required columns
required_columns <- c(
  "Accession",
  "Clade",
  "locus_tag",
  "genome_count",
  "total_genomes",
  "freq"
)

missing_columns <- setdiff(
  required_columns,
  colnames(vf_freq_by_clade)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# Confirm only A and B are present
unexpected_clades <- setdiff(
  unique(vf_freq_by_clade$Clade),
  c("A", "B")
)

if (length(unexpected_clades) > 0) {
  stop(
    "Unexpected clade values after recoding: ",
    paste(unexpected_clades, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Prepare genomes-per-clade table
# ------------------------------------------------------------

genomes_per_clade_named <- genomes_per_clade %>%
  mutate(
    Clade = as.character(Clade),
    
    Clade = recode(
      Clade,
      "1" = "A",
      "2" = "B"
    )
  )


# Confirm that genome totals are available for both clades
missing_clade_totals <- setdiff(
  c("A", "B"),
  unique(genomes_per_clade_named$Clade)
)

if (length(missing_clade_totals) > 0) {
  stop(
    "Missing genome totals for clade(s): ",
    paste(missing_clade_totals, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2. Create one Accession-to-locus_tag lookup
# ------------------------------------------------------------

vf_locus_lookup <- vf_freq_by_clade %>%
  select(
    Accession,
    locus_tag
  ) %>%
  arrange(
    Accession,
    is.na(locus_tag)
  ) %>%
  distinct(
    Accession,
    .keep_all = TRUE
  )


# Check whether any accession has multiple different locus tags
multiple_locus_tags <- vf_freq_by_clade %>%
  filter(
    !is.na(locus_tag)
  ) %>%
  distinct(
    Accession,
    locus_tag
  ) %>%
  count(
    Accession,
    name = "n_locus_tags"
  ) %>%
  filter(
    n_locus_tags > 1
  )

if (nrow(multiple_locus_tags) > 0) {
  warning(
    nrow(multiple_locus_tags),
    " accessions have more than one non-missing locus tag. ",
    "The first locus tag was retained."
  )
}


# ------------------------------------------------------------
# 3. Include VFs absent from one clade
# ------------------------------------------------------------

all_combinations <- vf_locus_lookup %>%
  crossing(
    Clade = c("A", "B")
  )


vf_freq_complete <- all_combinations %>%
  left_join(
    vf_freq_by_clade %>%
      select(
        -locus_tag
      ),
    by = c(
      "Accession",
      "Clade"
    )
  ) %>%
  left_join(
    genomes_per_clade_named,
    by = "Clade",
    suffix = c("", "_from_clade")
  ) %>%
  mutate(
    total_genomes = coalesce(
      total_genomes,
      total_genomes_from_clade
    ),
    
    genome_count = replace_na(
      genome_count,
      0L
    ),
    
    absent_count =
      total_genomes - genome_count,
    
    freq =
      genome_count / total_genomes
  ) %>%
  select(
    -total_genomes_from_clade
  )


# ------------------------------------------------------------
# 4. Reshape to one row per VF
# ------------------------------------------------------------

vf_freq_wide_test <- vf_freq_complete %>%
  select(
    Accession,
    locus_tag,
    Clade,
    genome_count,
    absent_count,
    total_genomes,
    freq
  ) %>%
  pivot_wider(
    id_cols = c(
      Accession,
      locus_tag
    ),
    
    names_from = Clade,
    
    values_from = c(
      genome_count,
      absent_count,
      total_genomes,
      freq
    ),
    
    names_glue = "{.value}_Clade_{Clade}"
  )


# ------------------------------------------------------------
# 5. Fisher exact test for each VF
# ------------------------------------------------------------

run_fisher <- function(
    present_A,
    absent_A,
    present_B,
    absent_B
) {
  
  contingency <- matrix(
    c(
      present_A,
      absent_A,
      present_B,
      absent_B
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Clade = c(
        "Clade_A",
        "Clade_B"
      ),
      Status = c(
        "Present",
        "Absent"
      )
    )
  )
  
  test <- fisher.test(
    contingency,
    alternative = "two.sided"
  )
  
  tibble(
    odds_ratio_CladeA_vs_CladeB =
      unname(test$estimate),
    
    p_value =
      test$p.value
  )
}


# ------------------------------------------------------------
# 6. Run enrichment analysis
# ------------------------------------------------------------

vf_clade_enrichment <- vf_freq_wide_test %>%
  rowwise() %>%
  mutate(
    fisher_result = list(
      run_fisher(
        present_A =
          genome_count_Clade_A,
        
        absent_A =
          absent_count_Clade_A,
        
        present_B =
          genome_count_Clade_B,
        
        absent_B =
          absent_count_Clade_B
      )
    )
  ) %>%
  unnest(
    fisher_result
  ) %>%
  ungroup() %>%
  mutate(
    freq_difference_CladeA_minus_CladeB =
      freq_Clade_A - freq_Clade_B,
    
    absolute_freq_difference =
      abs(
        freq_difference_CladeA_minus_CladeB
      ),
    
    freq_ratio = case_when(
      freq_Clade_A > 0 &
        freq_Clade_B > 0 ~
        pmax(
          freq_Clade_A,
          freq_Clade_B
        ) /
        pmin(
          freq_Clade_A,
          freq_Clade_B
        ),
      
      freq_Clade_A > 0 |
        freq_Clade_B > 0 ~
        Inf,
      
      TRUE ~
        NA_real_
    ),
    
    fdr_bh = p.adjust(
      p_value,
      method = "BH"
    ),
    
    # FDR significance alone
    significant_fdr =
      fdr_bh <= 0.05,
    
    # Final biological threshold:
    # FDR <= 0.05 and at least 5 percentage points different
    significant_and_meaningful =
      fdr_bh <= 0.05 &
      absolute_freq_difference >= 0.05,
    
    enriched_clade = case_when(
      significant_and_meaningful &
        freq_difference_CladeA_minus_CladeB > 0 ~
        "Clade_A",
      
      significant_and_meaningful &
        freq_difference_CladeA_minus_CladeB < 0 ~
        "Clade_B",
      
      significant_fdr &
        absolute_freq_difference < 0.05 ~
        "Significant_small_difference",
      
      TRUE ~
        "Not_significant"
    )
  ) %>%
  arrange(
    fdr_bh,
    desc(
      absolute_freq_difference
    )
  )


# ------------------------------------------------------------
# 7. Arrange output columns
# ------------------------------------------------------------

vf_clade_enrichment <- vf_clade_enrichment %>%
  relocate(
    locus_tag,
    .after = Accession
  ) %>%
  relocate(
    genome_count_Clade_A,
    genome_count_Clade_B,
    absent_count_Clade_A,
    absent_count_Clade_B,
    total_genomes_Clade_A,
    total_genomes_Clade_B,
    freq_Clade_A,
    freq_Clade_B,
    .after = locus_tag
  )


# ------------------------------------------------------------
# 8. Save all results
# ------------------------------------------------------------

write_csv(
  vf_clade_enrichment,
  "VF_presence_absence_enrichment_Clade_A_vs_B_all.csv",
  na = ""
)


write_csv(
  vf_clade_enrichment %>%
    filter(
      enriched_clade == "Clade_A"
    ),
  "VF_presence_enriched_in_Clade_A.csv",
  na = ""
)


write_csv(
  vf_clade_enrichment %>%
    filter(
      enriched_clade == "Clade_B"
    ),
  "VF_presence_enriched_in_Clade_B.csv",
  na = ""
)

write_csv(
  vf_clade_enrichment %>%
    filter(
      enriched_clade ==
        "Significant_small_difference"
    ),
  "VF_statistically_significant_small_difference.csv",
  na = ""
)

# ------------------------------------------------------------
# 9. Save unmatched locus-tag rows
# ------------------------------------------------------------

write_csv(
  vf_clade_enrichment %>%
    filter(
      is.na(locus_tag)
    ),
  "VF_enrichment_without_locus_tag.csv",
  na = ""
)


# ------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------

cat(
  "Total VFs tested:",
  nrow(vf_clade_enrichment),
  "\n"
)

cat(
  "VFs with locus tags:",
  sum(
    !is.na(
      vf_clade_enrichment$locus_tag
    )
  ),
  "\n"
)

cat(
  "VFs without locus tags:",
  sum(
    is.na(
      vf_clade_enrichment$locus_tag
    )
  ),
  "\n"
)

cat(
  "VFs enriched in Clade A:",
  sum(
    vf_clade_enrichment$enriched_clade ==
      "Clade_A"
  ),
  "\n"
)

cat(
  "VFs enriched in Clade B:",
  sum(
    vf_clade_enrichment$enriched_clade ==
      "Clade_B"
  ),
  "\n"
)

cat(
  "Not significant:",
  sum(
    vf_clade_enrichment$enriched_clade ==
      "Not_significant"
  ),
  "\n"
)

cat(
  "Statistically significant but <5% difference:",
  sum(
    vf_clade_enrichment$enriched_clade ==
      "Significant_small_difference"
  ),
  "\n"
)
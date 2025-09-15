nextflow.enable.dsl=2

// ----------- Inputs & folders -----------
params.out_root = "${params.work_dir}"

// Source channels
Channel
  .fromPath(params.vf_list)
  .set { vf_list_ch }

Channel
  .fromPath(params.label_csv)
  .set { label_csv_ch }

// ---------------- Processes ----------------

// Step 1: per-VF CSVs
process GREP_VF {
  tag { vf_list.baseName }    // or: tag "grep"
  cpus 16
  publishDir "${params.out_root}/VFs", mode: 'copy', overwrite: true

  input:
  path vf_list
  path label_csv

  output:
  path "VFs/*.csv", emit: vf_csvs

  script:
  """
  set -e
  mkdir -p VFs
  # Extract the VF IDs (1st col, skipping header) and emit one CSV per VF
  tail -n +2 ${vf_list} | cut -d',' -f1 | while read vf; do
    vf_safe=\$(echo "\$vf" | tr -d '\\r\\n')
    grep "\$vf_safe" ${label_csv} > VFs/\${vf_safe}.csv || true
  done
  """
}

// Step 2: CSV -> unaligned FASTA
process CSV_TO_FASTA {
  tag { csv.baseName }
  cpus 2
  publishDir "${params.out_root}/fastas", mode: 'copy', overwrite: true

  input:
  path csv

  output:
  path "*.fasta", emit: raw_fastas

  script:
  """
  python ${params.align_script} -i ${csv} -o .
  """
}

// Step 3: remove gaps
process REMOVE_GAPS {
  tag { fasta.baseName }
  cpus 1
  publishDir "${params.out_root}/fastas_nogap", mode: 'copy', overwrite: true

  input:
  path fasta

  output:
  // optional:true lets the task succeed even if no output (dropped item)
  path "*.nogap.fasta", optional: true, emit: nogap_fastas

  script:
  """
  set -euo pipefail
  base=\$(basename "${fasta}" .fasta)
  out="\${base}.nogap.fasta"
  tmp="\${out}.tmp"

  # Strip '-' only from sequence lines, keep headers as-is.
  # Track if any sequence remains non-empty after stripping.
  awk '
    /^>/ { print; next }
    {
      gsub(/[- \\t\\r]/,"",\$0);
      print;
      if (length(\$0) > 0) nonempty=1
    }
    END {
      if (!nonempty) exit 3
    }
  ' "${fasta}" > "\${tmp}" || status=\$?

  if [ "\${status:-0}" -eq 3 ]; then
    echo "[REMOVE_GAPS] ${fasta} -> no residues after removing gaps; skipping" >&2
    # do not create output => Nextflow will drop this item
    exit 0
  fi

  mv "\${tmp}" "\${out}"
  """
}



// Step 4: TranslatorX codon-aware alignment
process TRANSLATORX {
  tag { fasta.baseName }
  cpus 4
  publishDir "${params.out_root}/final_aln", mode: 'copy', overwrite: true

  input:
  path fasta

  output:
  path "${fasta.baseName}/*", emit: final_aln_dirs

  script:
  """
  set -euo pipefail
  base=\$(basename "${fasta}" .fasta)
  mkdir -p "\${base}"
  translatorx_vLocal.pl -i "${fasta}" -o "\${base}/\${base}" -p F -t T
  """
}


// Step 5: split alignments by trait
process SPLIT_ALIGNMENT {
  tag "split"
  cpus 2
  publishDir "${params.out_root}", mode: 'copy', overwrite: true

  input:
  path aln_dir  // staged here so the script sees all inputs in CWD

  output:
  path "0_splitted_aln/*.fasta", optional: true, emit: group0_fastas
  path "1_splitted_aln/*.fasta", optional: true, emit: group1_fastas

  script:
  """
  set -euo pipefail
  # The script expects -i as the root of final alignments (they are staged in CWD)
  python ${params.split_script} -m ${params.trait_file} -i . -o .
  """
}


// Step 6a: sanitize group0
process SANITIZE_G0 {
  tag { fasta.baseName }
  cpus 1
  publishDir "${params.out_root}/sanitized/g0", mode: 'copy', overwrite: true
  input:
    path fasta
  output:
    path "san/*.fasta", optional: true, emit: sanitized_fastas
  script:
  """
  set -euo pipefail
  mkdir -p san
  in="${fasta}"
  out="san/\$(basename "${fasta}")"

  # skip if the split FASTA is missing or empty
  if [ ! -s "\$in" ]; then
    echo "[SANITIZE_G0] \$in is empty; skipping" >&2
    exit 0
  fi

  bash "${params.sanitize_script}" "\$in" "\$out"
  """
}

// Step 6b: sanitize group1
process SANITIZE_G1 {
  tag { fasta.baseName }
  cpus 1
  publishDir "${params.out_root}/sanitized/g1", mode: 'copy', overwrite: true
  input:
    path fasta
  output:
    path "san/*.fasta", optional: true, emit: sanitized_fastas
  script:
  """
  set -euo pipefail
  mkdir -p san
  in="${fasta}"
  out="san/\$(basename "${fasta}")"

  if [ ! -s "\$in" ]; then
    echo "[SANITIZE_G1] \$in is empty; skipping" >&2
    exit 0
  fi

  bash "${params.sanitize_script}" "\$in" "\$out"
  """
}





// ---- Helper: pair filenames between groups by basename ----
def pairByBase(g0, g1) {
  g0.map { f -> tuple(f.baseName, f) }
    .join( g1.map { f -> tuple(f.baseName, f) } )
    .map  { id, a, b -> tuple(id, a, b) }
}


process DNDS_0v1 {
  tag { id }
  cpus 32
  publishDir(
    "${params.out_root}/dnds_output/0_vs_1",
    mode: 'copy',
    overwrite: true,
    pattern: "*_groupwise_dnds.csv",
    saveAs: { fname -> "${id}_0v1_${fname}" }   // <- unique per task
  )

  input:
  tuple val(id),
        path(a, stageAs: 'G0.fasta'),
        path(b, stageAs: 'G1.fasta')

  output:
  path "*_groupwise_dnds.csv"

  script:
  """
  set -euo pipefail
  python ${params.dnds_script} groupwise \\
      G0.fasta \\
      G1.fasta \\
      -o . -t 32 --consensus
  """
}

process DNDS_1v0 {
  tag { id }
  cpus 32
  publishDir(
    "${params.out_root}/dnds_output/1_vs_0",
    mode: 'copy',
    overwrite: true,
    pattern: "*_groupwise_dnds.csv",
    saveAs: { fname -> "${id}_1v0_${fname}" }
  )

  input:
  tuple val(id),
        path(a, stageAs: 'G0.fasta'),
        path(b, stageAs: 'G1.fasta')

  output:
  path "*_groupwise_dnds.csv"

  script:
  """
  set -euo pipefail
  python ${params.dnds_script} groupwise \\
      G1.fasta \\
      G0.fasta \\
      -o . -t 32 --consensus
  """
}

process CLEANUP_OUTPUTS {
  tag 'final-cleanup'
  cpus 1
  // nothing to publish — this only cleans intermediates
  input:
    val _barrier   // forces this to run once, after all DNDS files are produced

  output:
    path "CLEANUP_DONE"  // tiny marker so Nextflow knows this ran

  script:
  """
  set -euo pipefail

  # Keep ONLY dnds_output under ${params.out_root}
  ROOT="${params.out_root}"

  # Safety: ensure ROOT is not empty and exists
  if [ -z "\$ROOT" ] || [ ! -d "\$ROOT" ]; then
    echo "Invalid ROOT: '\$ROOT'" >&2
    exit 1
  fi

  # Remove intermediate folders if they exist
  rm -rf "\$ROOT/VFs" \
         "\$ROOT/fastas" \
         "\$ROOT/fastas_nogap" \
         "\$ROOT/final_aln" \
         "\$ROOT/0_splitted_aln" \
         "\$ROOT/1_splitted_aln" \
         "\$ROOT/sanitized"

  # Optionally remove reports if you don't want them
  # rm -rf "\$ROOT/.reports"

  # Leave: \$ROOT/dnds_output (your finals)
  : > CLEANUP_DONE
  """
}


// ---------------- Workflow wiring ----------------
workflow {
  // 1) GREP
  def vf_csvs         = GREP_VF(vf_list_ch, label_csv_ch)

  // 2) CSV -> FASTA
  def raw_fastas      = CSV_TO_FASTA(vf_csvs.flatten())

  // 3) Remove gaps
  def nogap_fastas    = REMOVE_GAPS(raw_fastas)

  // 4) TranslatorX
  def final_aln_dirs  = TRANSLATORX(nogap_fastas)

  // 5) Split
  def (group0_fastas, group1_fastas) = SPLIT_ALIGNMENT(final_aln_dirs.collect())

  // 6) Sanitize each group with distinct processes
  def g0_clean = SANITIZE_G0(group0_fastas.flatten())
  def g1_clean = SANITIZE_G1(group1_fastas.flatten())

  // 7) Pair & run dN/dS
  def paired_fastas = pairByBase(g0_clean, g1_clean)
  DNDS_0v1(paired_fastas)
  DNDS_1v0(paired_fastas)

  // ---- Barrier: wait for all DNDS files, then cleanup
  // Mix the two output channels and collect to a single token
  def barrier = out_0v1.mix(out_1v0).collect()
  CLEANUP_OUTPUTS(barrier)

}

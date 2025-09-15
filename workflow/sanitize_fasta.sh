#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <in.fasta> <out.fasta> [header_map.tsv]" >&2
  exit 1
fi

in="$1"
out="$2"
map="${3-}"

[[ -s "$in" ]] || { echo "Input file not found or empty: $in" >&2; exit 2; }

# If a map file was requested, start it fresh
if [[ -n "${map}" ]]; then
  : > "${map}"
fi

awk -v map_out="${map:-}" '
  BEGIN{
    OFS=""
    have_map=(map_out!="")
  }
  # Header line
  /^>/{
    h=$0
    sub(/^>/,"",h)
    split(h,a,/[ \t]+/)         # take first token
    id=a[1]
    gsub(/[^A-Za-z0-9._-]/,"_",id)   # gubbins-like sanitization
    cnt[id]++
    if(cnt[id]>1) new_id=sprintf("%s__dup%03d", id, cnt[id]-1)
    else          new_id=id
    print ">", new_id
    next
  }
  # Sequence line: replace any char not in ACGTNacgtn- with n
  {
    gsub(/[^ACGTNacgtn-]/,"n")
    print
  }
  END{
    if(have_map){
      for(k in cnt){
        if(cnt[k]==1) printf("%s\t%s\n", k, k) >> map_out
        else {
          for(i=1;i<=cnt[k];i++)
            printf("%s\t%s__dup%03d\n", k, k, i-1) >> map_out
        }
      }
    }
  }
' "${in}" > "${out}"

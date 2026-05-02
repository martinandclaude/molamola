#!/usr/bin/env bash
# Quick BND audit using bcftools - sanity check before plotting.
#
# Usage:
#   ./summarize_bnds.sh path/to/sample.sniffles.vcf[.gz]

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <vcf>" >&2
    exit 1
fi
VCF="$1"

if ! command -v bcftools >/dev/null; then
    echo "bcftools not found - run: conda activate molamola" >&2
    exit 1
fi

echo "=== Source VCF ==="
echo "$VCF"
echo

echo "=== SV type counts ==="
bcftools view -h "$VCF" >/dev/null  # validates header
bcftools view "$VCF" 2>/dev/null \
    | awk -F'\t' '/^[^#]/ { for (i=1;i<=NF;i++) if (match($i, /SVTYPE=[A-Z]+/)) {print substr($i, RSTART+7, RLENGTH-7); break} }' \
    | sort | uniq -c | sort -rn
echo

echo "=== BND filter status ==="
bcftools view -i 'INFO/SVTYPE="BND"' "$VCF" 2>/dev/null \
    | awk -F'\t' '/^[^#]/ {print $7}' | sort | uniq -c | sort -rn
echo

echo "=== BND inter- vs intra-chromosomal ==="
bcftools view -i 'INFO/SVTYPE="BND"' "$VCF" 2>/dev/null \
    | awk -F'\t' '
        /^[^#]/ {
            split($5, p, /[\[\]]/); split(p[2], mc, ":")
            if (mc[1] == $1) print "intra"; else print "inter"
        }' \
    | sort | uniq -c
echo

echo "=== Top 10 BND chromosome pairs (PASS only) ==="
bcftools view -f PASS -i 'INFO/SVTYPE="BND"' "$VCF" 2>/dev/null \
    | awk -F'\t' '
        /^[^#]/ {
            split($5, p, /[\[\]]/); split(p[2], mc, ":")
            a=$1; b=mc[1]
            if (a > b) { t=a; a=b; b=t }
            print a "\t" b
        }' \
    | sort | uniq -c | sort -rn | head -10

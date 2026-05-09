#!/bin/bash
# Download missing files: athal gff3, bnap gff3, ptri gff3
set -euo pipefail
OUTDIR="/data5/qiulei/PPI/data"
HTTPS_BASE="https://ftp.ensemblgenomes.org/pub/plants/release-57"
mkdir -p "$OUTDIR/logs"

download_file() {
    local url="$1" outfile="$2" logfile="$3" desc="$4"
    echo "Downloading $desc..."
    wget -c --no-check-certificate -t 3 -T 120 "$url" -O "$outfile" 2>> "$logfile" && {
        echo "  OK: $desc ($(du -h "$outfile" | cut -f1))"
    } || echo "  FAILED: $desc"
}

ATHAL_DIR="$OUTDIR/arabidopsis_thaliana"
mkdir -p "$ATHAL_DIR"
download_file "${HTTPS_BASE}/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gff3.gz" \
    "$ATHAL_DIR/arabidopsis_thaliana.gff3.gz" "$OUTDIR/logs/athal_download.log" "A. thaliana gff3"

BNAP_DIR="$OUTDIR/brassica_napus"
rm -f "$BNAP_DIR/brassica_napus.gff3.gz"
download_file "${HTTPS_BASE}/gff3/brassica_napus/Brassica_napus.AST_PRJEB5043_v1.57.gff3.gz" \
    "$BNAP_DIR/brassica_napus.gff3.gz" "$OUTDIR/logs/bnap_download.log" "B. napus gff3"

PTRI_DIR="$OUTDIR/populus_trichocarpa"
rm -f "$PTRI_DIR/populus_trichocarpa.gff3.gz"
download_file "${HTTPS_BASE}/gff3/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v4.57.gff3.gz" \
    "$PTRI_DIR/populus_trichocarpa.gff3.gz" "$OUTDIR/logs/ptri_download.log" "P. trichocarpa gff3"
echo "All missing files downloaded."

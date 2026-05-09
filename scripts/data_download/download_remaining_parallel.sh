#!/bin/bash
# Parallel download of remaining species files
OUTDIR="/data5/qiulei/PPI/data"
HTTPS_BASE="https://ftp.ensemblgenomes.org/pub/plants/release-57"
mkdir -p "$OUTDIR/logs"

download_file() {
    local sp="$1" url="$2" outfile="$3" label="$4"
    if [ -s "$outfile" ]; then echo "[$sp] $label exists"; return 0; fi
    wget -c --no-check-certificate -t 5 -T 180 "$url" -O "$outfile" 2>> "$OUTDIR/logs/${sp}.log"
    if [ -s "$outfile" ]; then echo "[$sp] $label done ($(du -h $outfile | cut -f1))"
    else rm -f "$outfile"; return 1; fi
}

# musa_acuminata gff3
sp="musa_acuminata"; mkdir -p "$OUTDIR/$sp"
download_file "$sp" "${HTTPS_BASE}/gff3/musa_acuminata/Musa_acuminata.Musa_acuminata_v2.57.gff3.gz" "$OUTDIR/$sp/musa_acuminata.gff3.gz" "gff3" &

# asparagus_officinalis
sp="asparagus_officinalis"; mkdir -p "$OUTDIR/$sp"
download_file "$sp" "${HTTPS_BASE}/fasta/asparagus_officinalis/pep/Asparagus_officinalis.Aspof.V1.pep.all.fa.gz" "$OUTDIR/$sp/asparagus_officinalis.pep.fa.gz" "pep" &
download_file "$sp" "${HTTPS_BASE}/fasta/asparagus_officinalis/cds/Asparagus_officinalis.Aspof.V1.cds.all.fa.gz" "$OUTDIR/$sp/asparagus_officinalis.cds.fa.gz" "cds" &
download_file "$sp" "${HTTPS_BASE}/gff3/asparagus_officinalis/Asparagus_officinalis.Aspof.V1.57.gff3.gz" "$OUTDIR/$sp/asparagus_officinalis.gff3.gz" "gff3" &

# nymphaea_colorata
sp="nymphaea_colorata"
mkdir -p "$OUTDIR/$sp"
download_file "$sp" "${HTTPS_BASE}/fasta/nymphaea_colorata/pep/Nymphaea_colorata.ASM883128v1.pep.all.fa.gz" "$OUTDIR/$sp/nymphaea_colorata.pep.fa.gz" "pep" &
download_file "$sp" "${HTTPS_BASE}/fasta/nymphaea_colorata/cds/Nymphaea_colorata.ASM883128v1.cds.all.fa.gz" "$OUTDIR/$sp/nymphaea_colorata.cds.fa.gz" "cds" &
download_file "$sp" "${HTTPS_BASE}/gff3/nymphaea_colorata/Nymphaea_colorata.ASM883128v1.57.gff3.gz" "$OUTDIR/$sp/nymphaea_colorata.gff3.gz" "gff3" &

# marchantia_polymorpha
sp="marchantia_polymorpha"
mkdir -p "$OUTDIR/$sp"
download_file "$sp" "${HTTPS_BASE}/fasta/marchantia_polymorpha/pep/Marchantia_polymorpha.Marchanta_polymorpha_v1.pep.all.fa.gz" "$OUTDIR/$sp/marchantia_polymorpha.pep.fa.gz" "pep" &
download_file "$sp" "${HTTPS_BASE}/fasta/marchantia_polymorpha/cds/Marchantia_polymorpha.Marchanta_polymorpha_v1.cds.all.fa.gz" "$OUTDIR/$sp/marchantia_polymorpha.cds.fa.gz" "cds" &
download_file "$sp" "${HTTPS_BASE}/gff3/marchantia_polymorpha/Marchantia_polymorpha.Marchanta_polymorpha_v1.57.gff3.gz" "$OUTDIR/$sp/marchantia_polymorpha.gff3.gz" "gff3" &

wait
echo "All parallel downloads complete."

#!/bin/bash
# Download species NOT on Ensembl Plants via NCBI Datasets CLI.
set -euo pipefail
OUTDIR="/data5/qiulei/PPI/data"
TMPDIR="/tmp/ncbi_download"
mkdir -p "$OUTDIR/logs" "$TMPDIR"

SPECIES=(
  "ginkgo_biloba:Ginkgo biloba:GCF_000287395.1:Gymnosperm"
  "pinus_taeda:Pinus taeda:GCF_000404065.3:Gymnosperm"
  "picea_abies:Picea abies:GCF_900067695.1:Gymnosperm"
  "ceratopteris_richardii:Ceratopteris richardii:GCF_001661285.1:Fern"
  "azolla_filiculoides:Azolla filiculoides:GCF_000148685.1:Fern"
  "persea_americana:Persea americana:GCF_025297015.1:Magnoliid"
  "piper_nigrum:Piper nigrum:GCA_004796395.1:Magnoliid"
  "aquilegia_coerulea:Aquilegia coerulea:GCF_000745375.1:Basal_Eudicot"
  "nelumbo_nucifera:Nelumbo nucifera:GCF_000365195.2:Basal_Eudicot"
  "phoenix_dactylifera:Phoenix dactylifera:GCF_009389715.2:Monocot"
  "phalaenopsis_equestris:Phalaenopsis equestris:GCF_001263595.1:Monocot"
  "elaeis_guineensis:Elaeis guineensis:GCF_000442705.1:Monocot"
  "spirodela_polyrhiza:Spirodela polyrhiza:GCF_002804225.1:Monocot"
  "zostera_marina:Zostera marina:GCF_001185155.1:Monocot"
  "fragaria_vesca:Fragaria vesca:GCF_000184155.2:Rosid"
  "carica_papaya:Carica papaya:GCF_000150535.2:Rosid"
  "ricinus_communis:Ricinus communis:GCF_000151685.1:Rosid"
  "arachis_hypogaea:Arachis hypogaea:GCF_003086295.2:Rosid"
  "petunia_axillaris:Petunia axillaris:GCF_000223135.1:Asterid"
  "cuscuta_australis:Cuscuta australis:GCA_002007495.1:Asterid"
  "catharanthus_roseus:Catharanthus roseus:GCA_001951915.1:Asterid"
  "salvia_miltiorrhiza:Salvia miltiorrhiza:GCF_002888555.2:Asterid"
)

download_ncbi() {
    local sp="$1" sci="$2" acc="$3" clade="$4"
    local sdir="$OUTDIR/$sp"
    mkdir -p "$sdir" "$TMPDIR/$sp"
    local existing=$(ls "$sdir"/*.gz 2>/dev/null | wc -l)
    [ "$existing" -ge 3 ] && echo "$sp already has $existing files" && return 0
    echo "Downloading: $sci ($sp) [$clade]"
    datasets download genome accession "$acc" --include gff3,cds,protein --filename "$TMPDIR/${sp}.zip" 2>> "$OUTDIR/logs/${sp}.log"
    [ ! -s "$TMPDIR/${sp}.zip" ] && echo "FAILED" && return 1
    rm -rf "$TMPDIR/$sp" && unzip -q -o "$TMPDIR/${sp}.zip" -d "$TMPDIR/$sp" 2>> "$OUTDIR/logs/${sp}.log"
    local pep=$(find "$TMPDIR/$sp" -name "*protein.faa" | head -1)
    local cds=$(find "$TMPDIR/$sp" -name "*cds_from_genomic.fna" | head -1)
    local gff=$(find "$TMPDIR/$sp" -name "*.gff" | head -1)
    [ -n "$pep" ] && gzip -c "$pep" > "$sdir/${sp}.pep.fa.gz"
    [ -n "$cds" ] && gzip -c "$cds" > "$sdir/${sp}.cds.fa.gz"
    [ -n "$gff" ] && gzip -c "$gff" > "$sdir/${sp}.gff3.gz"
    rm -f "$TMPDIR/${sp}.zip"
}

for entry in "${SPECIES[@]}"; do
    IFS=: read -r sp sci acc clade <<< "$entry"
    download_ncbi "$sp" "$sci" "$acc" "$clade"
    sleep 3
done
echo "All NCBI downloads complete."

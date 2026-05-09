#!/bin/bash
# Download additional species to fill phylogenetic gaps (Ensembl Plants release 57)
set -euo pipefail
OUTDIR="/data5/qiulei/PPI/data"
HTTPS_BASE="https://ftp.ensemblgenomes.org/pub/plants/release-57"
mkdir -p "$OUTDIR/logs"

download_species() {
    local sp="$1" assembly="$2" desc="$3"
    local sdir="$OUTDIR/$sp"
    mkdir -p "$sdir"
    local genus_upper=$(echo "$sp" | sed 's/_.*//' | sed 's/^./\U&/')
    local species_lower=$(echo "$sp" | sed 's/^[^_]*_//')
    local gsu="${genus_upper}_${species_lower}"
    echo "Downloading: $desc ($sp)"
    for type in pep cds gff3; do
        local suffix ext url out
        case $type in
            pep)  suffix="pep.all.fa.gz"; ext="pep.fa.gz"
                  url="${HTTPS_BASE}/fasta/${sp}/pep/${gsu}.${assembly}.${suffix}"
                  out="${sdir}/${sp}.${ext}" ;;
            cds)  suffix="cds.all.fa.gz"; ext="cds.fa.gz"
                  url="${HTTPS_BASE}/fasta/${sp}/cds/${gsu}.${assembly}.${suffix}"
                  out="${sdir}/${sp}.${ext}" ;;
            gff3) suffix="57.gff3.gz"; ext="gff3.gz"
                  url="${HTTPS_BASE}/gff3/${sp}/${gsu}.${assembly}.${suffix}"
                  out="${sdir}/${sp}.${ext}" ;;
        esac
        if [ -s "$out" ]; then echo "  $type exists"; else
            wget -c --no-check-certificate -t 3 -T 120 "$url" -O "$out" 2>> "$OUTDIR/logs/${sp}.log"
            [ -s "$out" ] && echo "  $type done" || echo "  $type FAILED"
        fi
    done
}

download_species "manihot_esculenta" "Manihot_esculenta_v6" "Cassava"
download_species "citrus_clementina" "Citrus_clementina_v1.0" "Clementine"
download_species "pistacia_vera" "PisVer_v2" "Pistachio"
download_species "malus_domestica_golden" "ASM211411v1" "Apple"
download_species "prunus_persica" "Prunus_persica_NCBIv2" "Peach"
download_species "rosa_chinensis" "RchiOBHm-V2" "Rose"
download_species "ficus_carica" "UNIPI_FiCari_1.0" "Fig"
download_species "medicago_truncatula" "MedtrA17_4.0" "Barrelclover"
download_species "phaseolus_vulgaris" "PhaVulg1_0" "Common bean"
download_species "theobroma_cacao" "Theobroma_cacao_20110822" "Cacao"
download_species "corylus_avellana" "CavTom2PMs-1.0" "Hazel"
download_species "quercus_lobata" "ValleyOak3.0" "Oak"
download_species "eucalyptus_grandis" "Egrandis1_0" "Eucalyptus"
download_species "ananas_comosus" "F153" "Pineapple"
download_species "musa_acuminata" "Musa_acuminata_v2" "Banana"
download_species "asparagus_officinalis" "Aspof.V1" "Asparagus"
download_species "nymphaea_colorata" "ASM883128v1" "Water lily"
download_species "marchantia_polymorpha" "Marchanta_polymorpha_v1" "Liverwort"
echo "All downloads complete."

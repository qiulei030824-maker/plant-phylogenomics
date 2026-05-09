#!/bin/bash
# ============================================================================
# download_remaining_robust.sh
#
# Robust download for remaining species using aria2c (HTTPS) + tmux.
# Features: checkpoint, retry, continue, resume.
#
# Ensembl Plants + NCBI Datasets
# ============================================================================
set -euo pipefail

OUTDIR="/data5/qiulei/PPI/data"
TMPDIR="/tmp/download_robust"
CHECKPOINT="$TMPDIR/checkpoint.txt"
LOGDIR="$OUTDIR/logs"
mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

HTTPS_BASE="https://ftp.ensemblgenomes.org/pub/plants/release-57"

# Remaining Ensembl species
declare -A ENSEMBL_REMAINING=(
  ["actinidia_chinensis"]=""
  ["aegilops_tauschii"]=""
  ["cajanus_cajan"]=""
  ["chenopodium_quinoa"]=""
  ["coffea_canephora"]=""
  ["corchorus_capsularis"]=""
  ["corymbia_citriodora"]=""
  ["cynara_cardunculus"]=""
  ["daucus_carota"]=""
  ["digitaria_exilis"]=""
  ["dioscorea_rotundata"]=""
  ["eragrostis_tef"]=""
  ["fraxinus_excelsior"]=""
  ["hordeum_vulgare"]=""
  ["ipomoea_triloba"]=""
  ["juglans_regia"]=""
  ["lactuca_sativa"]=""
  ["leersia_perrieri"]=""
  ["lolium_perenne"]=""
  ["lupinus_angustifolius"]=""
  ["panicum_hallii"]=""
  ["pisum_sativum"]=""
  ["prunus_avium"]=""
  ["prunus_dulcis"]=""
  ["quercus_suber"]=""
  ["secale_cereale"]=""
  ["sesamum_indicum"]=""
  ["setaria_viridis"]=""
  ["theobroma_cacao_criollo"]=""
  ["trifolium_pratense"]=""
  ["triticum_aestivum"]=""
  ["triticum_turgidum"]=""
  ["triticum_urartu"]=""
  ["vigna_angularis"]=""
  ["vigna_radiata"]=""
  ["vigna_unguiculata"]=""
)

# NCBI remaining species
declare -A NCBI_REMAINING=(
  ["picea_abies"]="GCF_900067695.1"
  ["ceratopteris_richardii"]="GCF_001661285.1"
  ["azolla_filiculoides"]="GCF_000148685.1"
  ["persea_americana"]="GCF_025297015.1"
  ["piper_nigrum"]="GCA_004796395.1"
  ["aquilegia_coerulea"]="GCF_000745375.1"
  ["nelumbo_nucifera"]="GCF_000365195.2"
  ["phoenix_dactylifera"]="GCF_009389715.2"
  ["phalaenopsis_equestris"]="GCF_001263595.1"
  ["elaeis_guineensis"]="GCF_000442705.1"
  ["spirodela_polyrhiza"]="GCF_002804225.1"
  ["zostera_marina"]="GCF_001185155.1"
  ["fragaria_vesca"]="GCF_000184155.2"
  ["carica_papaya"]="GCF_000150535.2"
  ["ricinus_communis"]="GCF_000151685.1"
  ["arachis_hypogaea"]="GCF_003086295.2"
  ["petunia_axillaris"]="GCF_000223135.1"
  ["cuscuta_australis"]="GCA_002007495.1"
  ["catharanthus_roseus"]="GCA_001951915.1"
  ["salvia_miltiorrhiza"]="GCF_002888555.2"
)

already_have() {
    local sp="$1"
    local c=$(ls "$OUTDIR/$sp/"*.gz 2>/dev/null | wc -l)
    [ "$c" -ge 3 ]
}

checkpoint_done() {
    local sp="$1"; local src="$2"
    echo "$src:$sp:done" >> "$CHECKPOINT"
}
is_done() {
    local sp="$1"; local src="$2"
    grep -q "^$src:$sp:done$" "$CHECKPOINT" 2>/dev/null
}

aria2_download() {
    local url="$1"; local outfile="$2"; local logfile="$3"; local desc="$4"
    if [ -s "$outfile" ]; then
        echo "  $desc already exists ($(du -h "$outfile" | cut -f1))"
        return 0
    fi
    echo "  Downloading $desc..."
    aria2c -c -x 4 -s 4 --retry-wait=5 --max-tries=10 --check-certificate=false \
        --connect-timeout=30 --timeout=300 \
        --dir="$(dirname "$outfile")" --out="$(basename "$outfile")" \
        "$url" 2>> "$logfile"
    if [ -s "$outfile" ]; then
        echo "  $desc done ($(du -h "$outfile" | cut -f1))"
        return 0
    fi
    return 1
}

download_ensembl_remaining() {
    echo "Phase 1: Ensembl Plants (HTTPS + aria2c)"
    for sp in "${!ENSEMBL_REMAINING[@]}"; do
        is_done "$sp" "ensembl" && echo "  $sp checkpoint done" && continue
        already_have "$sp" && checkpoint_done "$sp" "ensembl" && continue
        mkdir -p "$OUTDIR/$sp"
        local pep_file=$(wget --no-check-certificate -q -O- "${HTTPS_BASE}/fasta/${sp}/pep/" 2>/dev/null | grep -oP 'href="[^"]+\.pep\.all\.fa\.gz"' | head -1 | sed 's/href="//;s/"//')
        local cds_file=$(wget --no-check-certificate -q -O- "${HTTPS_BASE}/fasta/${sp}/cds/" 2>/dev/null | grep -oP 'href="[^"]+\.cds\.all\.fa\.gz"' | head -1 | sed 's/href="//;s/"//')
        local gff_file=$(wget --no-check-certificate -q -O- "${HTTPS_BASE}/gff3/${sp}/" 2>/dev/null | grep -oP 'href="[^"]+\.gff3\.gz"' | head -1 | sed 's/href="//;s/"//')
        local ok=0
        [ -n "$pep_file" ] && aria2_download "${HTTPS_BASE}/fasta/${sp}/pep/${pep_file}" "$OUTDIR/$sp/${sp}.pep.fa.gz" "$LOGDIR/${sp}.log" "pep" && ok=$((ok+1))
        [ -n "$cds_file" ] && aria2_download "${HTTPS_BASE}/fasta/${sp}/cds/${cds_file}" "$OUTDIR/$sp/${sp}.cds.fa.gz" "$LOGDIR/${sp}.log" "cds" && ok=$((ok+1))
        [ -n "$gff_file" ] && aria2_download "${HTTPS_BASE}/gff3/${sp}/${gff_file}" "$OUTDIR/$sp/${sp}.gff3.gz" "$LOGDIR/${sp}.log" "gff3" && ok=$((ok+1))
        [ "$ok" -ge 3 ] && checkpoint_done "$sp" "ensembl"
    done
}

download_ncbi_small() {
    echo "Phase 2: NCBI small genomes (datasets CLI)"
    for sp in "${!NCBI_REMAINING[@]}"; do
        local acc="${NCBI_REMAINING[$sp]}"
        is_done "$sp" "ncbi" && continue
        already_have "$sp" && checkpoint_done "$sp" "ncbi" && continue
        mkdir -p "$OUTDIR/$sp" "$TMPDIR/$sp"
        datasets download genome accession "$acc" --include gff3,cds,protein --filename "$TMPDIR/${sp}.zip" 2>> "$LOGDIR/${sp}.log" || {
            sleep 10 && datasets download genome accession "$acc" --include gff3,cds,protein --filename "$TMPDIR/${sp}.zip" 2>> "$LOGDIR/${sp}.log" || continue
        }
        [ ! -s "$TMPDIR/${sp}.zip" ] && continue
        rm -rf "$TMPDIR/$sp" && unzip -q -o "$TMPDIR/${sp}.zip" -d "$TMPDIR/$sp" 2>> "$LOGDIR/${sp}.log"
        local pep_src=$(find "$TMPDIR/$sp" -name "*protein.faa" | head -1)
        local cds_src=$(find "$TMPDIR/$sp" -name "*cds_from_genomic.fna" | head -1)
        local gff_src=$(find "$TMPDIR/$sp" -name "*.gff" | head -1)
        local ok=0
        [ -n "$pep_src" ] && gzip -c "$pep_src" > "$OUTDIR/$sp/${sp}.pep.fa.gz" && ok=$((ok+1))
        [ -n "$cds_src" ] && gzip -c "$cds_src" > "$OUTDIR/$sp/${sp}.cds.fa.gz" && ok=$((ok+1))
        [ -n "$gff_src" ] && gzip -c "$gff_src" > "$OUTDIR/$sp/${sp}.gff3.gz" && ok=$((ok+1))
        [ "$ok" -ge 3 ] && checkpoint_done "$sp" "ncbi" && rm -f "$TMPDIR/${sp}.zip"
    done
}

download_pinus() {
    echo "Phase 3: Pinus taeda (large genome ~22Gb, pep+gff only)"
    is_done "pinus_taeda" "ncbi" && return
    mkdir -p "$OUTDIR/pinus_taeda"
    datasets download genome accession "GCF_000404065.3" --include gff3,protein --filename "$TMPDIR/pinus.zip" 2>> "$LOGDIR/pinus_taeda.log" || {
        sleep 30 && datasets download genome accession "GCF_000404065.3" --include gff3,protein --filename "$TMPDIR/pinus.zip" 2>> "$LOGDIR/pinus_taeda.log" || return
    }
    [ -s "$TMPDIR/pinus.zip" ] && {
        rm -rf "$TMPDIR/pinus" && unzip -q -o "$TMPDIR/pinus.zip" -d "$TMPDIR/pinus"
        local pep=$(find "$TMPDIR/pinus" -name "*protein.faa" | head -1)
        local gff=$(find "$TMPDIR/pinus" -name "*.gff" | head -1)
        [ -n "$pep" ] && gzip -c "$pep" > "$OUTDIR/pinus_taeda/pinus_taeda.pep.fa.gz"
        [ -n "$gff" ] && gzip -c "$gff" > "$OUTDIR/pinus_taeda/pinus_taeda.gff3.gz"
        checkpoint_done "pinus_taeda" "ncbi" && rm -f "$TMPDIR/pinus.zip"
    }
}

download_ensembl_remaining
download_ncbi_small
download_pinus
echo "All downloads complete."

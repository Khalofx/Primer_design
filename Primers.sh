#!/bin/bash

echo "🧬 AutoPrimer3 Bash: Internal + External Split-Location Primer Design with FASTA & BED"

# Step 1: User Inputs
read -p "Enter genome assembly (e.g., hg38): " GENOME
read -p "Enter chromosome (e.g., chr1): " CHR
read -p "Enter start position: " START
read -p "Enter end position: " END
read -p "Enter INTERNAL product size range (e.g., 300-1000): " PRODUCT_SIZE_INTERNAL

# Fixed rules for external primers
EXTERNAL_FLANK=300
EXTERNAL_TARGET=700
PRODUCT_SIZE_EXTERNAL="400-1000"

# Output dir
TIMESTAMP=$(date +"%Y-%m-%d_%H%M")
OUTDIR="primer3_output_${TIMESTAMP}"
mkdir -p "$OUTDIR"

# Fetch from UCSC
fetch_seq() {
  local REGION=$1
  local FILE=$2
  local URL="http://genome.ucsc.edu/cgi-bin/das/${GENOME}/dna?segment=${CHR}:${REGION}"
  echo "🔄 Fetching $REGION from UCSC..."
  curl -s "$URL" | grep -v "dna" | tr -d '\n' | sed 's/.*<DNA.*>\(.*\)<\/DNA>.*/\1/' | tr -d ' \t\n\r' | tr 'a-z' 'A-Z' > "$FILE"
}

# Make Primer3 input
make_input() {
  local SEQ=$1
  local FILE=$2
  local ID=$3
  local LEFT=$4
  local RIGHT=$5
  local LEFT_EXCL=$6
  local RIGHT_EXCL=$7
  local RANGE=$8

  cat > "$FILE" <<EOF
SEQUENCE_ID=$ID
SEQUENCE_TEMPLATE=$SEQ
PRIMER_TASK=pick_pcr_primers
PRIMER_PICK_LEFT_PRIMER=$LEFT
PRIMER_PICK_RIGHT_PRIMER=$RIGHT
SEQUENCE_PRIMER_LEFT_EXCL=$LEFT_EXCL
SEQUENCE_PRIMER_RIGHT_EXCL=$RIGHT_EXCL
PRIMER_NUM_RETURN=5
PRIMER_PRODUCT_SIZE_RANGE=$RANGE
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=23
PRIMER_MIN_TM=57.0
PRIMER_OPT_TM=60.0
PRIMER_MAX_TM=63.0
PRIMER_MAX_POLY_X=4
PRIMER_SALT_MONOVALENT=50.0
PRIMER_DNA_CONC=50.0
PRIMER_EXPLAIN_FLAG=1
=
EOF
}

# Parse Primer3 output to TSV
parse_output() {
  local FILE=$1
  local TSV=$2
  local OFFSET=$3
  local TYPE=$4

  awk -v start=$OFFSET -v label="$TYPE" '
  BEGIN {
    print "Primer Pair\tLeft Primer\tRight Primer\tProduct Size\tLeft Tm\tRight Tm\tLeft GC%\tRight GC%\tAmplicon Start\tAmplicon End"
  }
  /^PRIMER_LEFT_[0-9]+_SEQUENCE=/ {split($0,a,"="); left=a[2]}
  /^PRIMER_RIGHT_[0-9]+_SEQUENCE=/ {split($0,a,"="); right=a[2]}
  /^PRIMER_LEFT_[0-9]+=/ {split($0,a,"="); split(a[2],pos,","); left_start=pos[1]}
  /^PRIMER_RIGHT_[0-9]+=/ {split($0,a,"="); split(a[2],pos,","); right_start=pos[1]}
  /^PRIMER_LEFT_[0-9]+_TM=/ {split($0,a,"="); left_tm=a[2]}
  /^PRIMER_RIGHT_[0-9]+_TM=/ {split($0,a,"="); right_tm=a[2]}
  /^PRIMER_LEFT_[0-9]+_GC_PERCENT=/ {split($0,a,"="); left_gc=a[2]}
  /^PRIMER_RIGHT_[0-9]+_GC_PERCENT=/ {split($0,a,"="); right_gc=a[2]}
  /^PRIMER_PAIR_[0-9]+_PRODUCT_SIZE=/ {
    split($0,a,"=");
    size=a[2];
    amplicon_start = start + left_start;
    amplicon_end = amplicon_start + size - 1;
    printf("Pair %d\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\n", ++pair, left, right, size, left_tm, right_tm, left_gc, right_gc, amplicon_start, amplicon_end)
  }' "$FILE" > "$TSV"
}

# FASTA output
generate_fasta() {
  local TSV=$1
  local OUT_FASTA=$2
  local TYPE=$3
  awk -v label="$TYPE" 'BEGIN{FS="\t"} NR>1 {
    print ">"label"_Pair" NR-1 "_left\n" $2
    print ">"label"_Pair" NR-1 "_right\n" $3
  }' "$TSV" > "$OUT_FASTA"
}

# BED output
generate_bed() {
  local TSV=$1
  local OUT_BED=$2
  local TYPE=$3
  local CHR=$4
  awk -v chr="$CHR" -v label="$TYPE" 'BEGIN{FS="\t"} NR>1 {
    left_start=$9
    left_end=left_start+length($2)
    right_end=$10
    right_start=right_end-length($3)
    print chr"\t"left_start"\t"left_end"\t"label"_Pair" NR-1 "_left"
    print chr"\t"right_start"\t"right_end"\t"label"_Pair" NR-1 "_right"
  }' "$TSV" > "$OUT_BED"
}

### 🧬 INTERNAL PRIMERS
fetch_seq "${START},${END}" "$OUTDIR/internal_seq.txt"
SEQ_INTERNAL=$(<"$OUTDIR/internal_seq.txt")
make_input "$SEQ_INTERNAL" "$OUTDIR/internal_input.txt" "internal" 1 1 "" "" "$PRODUCT_SIZE_INTERNAL"
primer3_core < "$OUTDIR/internal_input.txt" > "$OUTDIR/internal_output.txt"
parse_output "$OUTDIR/internal_output.txt" "$OUTDIR/internal.tsv" "$START" "internal"

### 🔹 UPSTREAM FLANKING PRIMERS
UP_FLANK_START=$((START - EXTERNAL_FLANK))
UP_TARGET_END=$((START + EXTERNAL_TARGET - 1))
fetch_seq "${UP_FLANK_START},${UP_TARGET_END}" "$OUTDIR/upstream_seq.txt"
SEQ_UP=$(<"$OUTDIR/upstream_seq.txt")
EXCL_LEFT_UP="${EXTERNAL_FLANK},${EXTERNAL_TARGET}"
EXCL_RIGHT_UP="0,${EXTERNAL_FLANK}"
make_input "$SEQ_UP" "$OUTDIR/upstream_input.txt" "upstream" 1 1 "$EXCL_LEFT_UP" "$EXCL_RIGHT_UP" "$PRODUCT_SIZE_EXTERNAL"
primer3_core < "$OUTDIR/upstream_input.txt" > "$OUTDIR/upstream_output.txt"
parse_output "$OUTDIR/upstream_output.txt" "$OUTDIR/upstream.tsv" "$UP_FLANK_START" "upstream"

### 🔹 DOWNSTREAM FLANKING PRIMERS
DOWN_TARGET_START=$((END - EXTERNAL_TARGET + 1))
DOWN_FLANK_END=$((END + EXTERNAL_FLANK))
fetch_seq "${DOWN_TARGET_START},${DOWN_FLANK_END}" "$OUTDIR/downstream_seq.txt"
SEQ_DOWN=$(<"$OUTDIR/downstream_seq.txt")
EXCL_LEFT_DOWN="${EXTERNAL_TARGET},${EXTERNAL_FLANK}"
EXCL_RIGHT_DOWN="0,${EXTERNAL_TARGET}"
make_input "$SEQ_DOWN" "$OUTDIR/downstream_input.txt" "downstream" 1 1 "$EXCL_LEFT_DOWN" "$EXCL_RIGHT_DOWN" "$PRODUCT_SIZE_EXTERNAL"
primer3_core < "$OUTDIR/downstream_input.txt" > "$OUTDIR/downstream_output.txt"
parse_output "$OUTDIR/downstream_output.txt" "$OUTDIR/downstream.tsv" "$DOWN_TARGET_START" "downstream"

### 🔄 FASTA & BED EXPORT
generate_fasta "$OUTDIR/internal.tsv" "$OUTDIR/internal_primers.fa" "internal"
generate_bed "$OUTDIR/internal.tsv" "$OUTDIR/internal_primers.bed" "internal" "$CHR"

generate_fasta "$OUTDIR/upstream.tsv" "$OUTDIR/upstream_primers.fa" "upstream"
generate_bed "$OUTDIR/upstream.tsv" "$OUTDIR/upstream_primers.bed" "upstream" "$CHR"

generate_fasta "$OUTDIR/downstream.tsv" "$OUTDIR/downstream_primers.fa" "downstream"
generate_bed "$OUTDIR/downstream.tsv" "$OUTDIR/downstream_primers.bed" "downstream" "$CHR"

echo "✅ All done. Output saved in: $OUTDIR"
ls "$OUTDIR"

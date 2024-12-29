#!/bin/bash
while getopts "i:o:a:r:w:" opt; do
  case $opt in
    i) input=$OPTARG ;;
    o) output=$OPTARG ;;
    a) adapter=$OPTARG ;;
    r) references=$OPTARG ;;
    w) wetlab=$OPTARG ;;  # Ensure this is defined
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1;;
  esac
done


# Paths to scripts
sampleSheet="/home/lannguyen/Documents/wgs_suran/process/update/1.read_sampleSheet.R"
QC="/home/lannguyen/Documents/wgs_suran/process/update/2.QC.sh"
assembly="/home/lannguyen/Documents/wgs_suran/process/update/3.assembly.sh"
samTobam="/home/lannguyen/Documents/wgs_suran/process/update/4.samTobam.sh"
correctBase="/home/lannguyen/Documents/wgs_suran/process/update/5.correct_base.sh"
callVariant="/home/lannguyen/Documents/wgs_suran/process/update/6.call_variant.sh"

# Run the R script for sample validation
Rscript --vanilla "$sampleSheet" "$input" "$output"

# Run QC script
bash "$QC" -i "$input" -o "$output" -a "$adapter"

# Run assembly script
bash "$assembly" -o "$output" -r "$references"

# Run samTobam script
bash "$samTobam" -o "$output" 

# Run correct base script
bash "$correctBase" -o "$output" -r "$references" -w "$wetlab"

# Call variant script
bash "$callVariant" -o "$output" -r "$references"

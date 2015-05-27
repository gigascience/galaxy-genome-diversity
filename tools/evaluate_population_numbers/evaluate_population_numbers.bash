#!/usr/bin/env bash

if [ $# -ne 3 ]; then
    echo "usage"
    exit 1
fi

input_ped_file="$1"
output_file="$2"
max_populations="$3"

ADMIXTURE=admixture

for (( i=1; $i <= $max_populations; i++ )); do
    $ADMIXTURE --cv "$input_ped_file" $i 2>&1 | grep CV | perl -ne 's/CV error/CVE/; print;' >> "$output_file"
done


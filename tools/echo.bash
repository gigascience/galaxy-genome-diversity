#!/usr/bin/env bash

if [ $# -lt 3 ]; then
    echo "usage"
    exit 1
fi

input="$1"
output="$2"
shift 2

for individual in "$@"; do
    echo "$individual" >> "$output"
done

exit 0


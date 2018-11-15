#!/bin/bash

expected_file=$1
computed_file=$2

# For the oracle, check only that the final column (the error) matches
# Paste the table file and the expected table file together as to do a line by line diff
IFS=$'\n'
rownum=1    # Keep track of rows
for l in $(paste <(grep -v AIMGM ${computed_file} | grep -v BD-ADI) <(grep -v AIMGM ${expected_file} | grep -v BD-ADI)); do
    # Do separate diff for each row (treat each row as a test)
    actual=$(echo ${l} | cut -d$'\t' -f1)
    expected=$(echo ${l} | cut -d$'\t' -f2)
    # Fail if last column does not match, exit with error code on fail
    actualval=$(echo ${actual} | cut -d',' -f13)
    expectedval=$(echo ${expected} | cut -d',' -f13)
    if [[ $(diff <(echo ${actualval}) <(echo ${expectedval})) != "" ]]; then
        echo "${computed_file} ROW ${rownum}: FAIL, expected ${expectedval}, got ${actualval}"
        exit 1
    fi
    rownum=$((rownum + 1))
done

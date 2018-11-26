#!/bin/bash

expected_file=$1
computed_file=$2

# Helper function for rounding numbers to same place in comparison
round() {
    echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
}
export -f round

# For the oracle, check only that the final column (the error) matches
# Paste the table file and the expected table file together as to do a line by line diff
IFS=$'\n'
rownum=1    # Keep track of rows
exitcode=0  # Keep track of error code for reporting final failure or not
for l in $(paste <(grep -v AIMGM ${computed_file} | grep -v BD-ADI) <(grep -v AIMGM ${expected_file} | grep -v BD-ADI)); do
    # Do separate diff for each row (treat each row as a test)
    actual=$(echo ${l} | cut -d$'\t' -f1)
    expected=$(echo ${l} | cut -d$'\t' -f2)
    # Fail if last column does not match, exit with error code on fail,
    # but round to ignore final two digits, which matters when comparing against paper as it does not go as precise
    actualval=$(echo ${actual} | cut -d',' -f13 | cut -c 1-8 | xargs -I{} bash -c 'round {} 4')$(echo ${actual} | cut -d',' -f13 | cut -c 9-)
    expectedval=$(echo ${expected} | cut -d',' -f13 | cut -c 1-8 | xargs -I{} bash -c 'round {} 4')$(echo ${actual} | cut -d',' -f13 | cut -c 9-)
    if [[ $(diff <(echo ${actualval}) <(echo ${expectedval})) != "" ]]; then
        echo "${computed_file} ROW ${rownum}: FAIL, expected ${expectedval}, got ${actualval}"
        exitcode=1
    fi
    rownum=$((rownum + 1))
done

exit ${exitcode}

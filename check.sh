#!/bin/bash

# For the oracle, check only that the final column (the error) matches
for i in 1 2 3 4; do
    # Paste the table file and the expected table file together as to do a line by line diff
    IFS=$'\n'
    rownum=1    # Keep track of rows
    for l in $(paste <(grep -v AIMGM table${i}.csv | grep -v BD-ADI) <(grep -v AIMGM expected_table${i}.csv | grep -v BD-ADI)); do
        # Do separate diff for each row (treat each row as a test)
        actual=$(echo ${l} | cut -d$'\t' -f1)
        expected=$(echo ${l} | cut -d$'\t' -f2)
        # Fail if last column does not match
        if [[ $(diff <(echo ${actual} | cut -d',' -f13) <(echo ${expected} | cut -d',' -f13)) != "" ]]; then
            echo "Table ${i} Row ${rownum}: FAIL"
        else
            echo "Table ${i} ROW ${rownum}: PASS"
        fi
        rownum=$((rownum + 1))
    done
done

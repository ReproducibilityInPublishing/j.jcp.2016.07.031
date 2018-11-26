#!/bin/bash
set -e

tabnum=$1
paper=$2

# Navigate into this directory
cd "$( dirname "${BASH_SOURCE[0]}" )"

echo "RUNNING TABLE ${tabnum}" && bash table${tabnum}_smoke.sh

# Check their results, compare against paper or not
if [[ ${paper} == 1 ]]; then
    bash check.sh expected_table${tabnum}_smoke_paper.csv table${tabnum}_smoke.csv
else
    bash check.sh expected_table${tabnum}_smoke.csv table${tabnum}_smoke.csv
fi

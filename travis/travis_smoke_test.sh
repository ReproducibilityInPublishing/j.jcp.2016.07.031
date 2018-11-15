#!/bin/bash
set -e

# Navigate into this directory
cd "$( dirname "${BASH_SOURCE[0]}" )"

echo "RUNNING TABLE 1" && bash table1_smoke.sh
echo "RUNNING TABLE 2" && bash table2_smoke.sh

# Check their results
bash check.sh expected_table1_smoke.csv table1_smoke.csv
bash check.sh expected_table2_smoke.csv table2_smoke.csv

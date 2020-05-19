#!/bin/bash

# Run all tables
bash table1_short.sh
bash table2_short.sh
bash table3_short.sh
bash table4_short.sh

# Check against the expected tables
python check_short.py

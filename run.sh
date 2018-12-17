make clean; make

# Run all tables
bash table1.sh
bash table2.sh
bash table3.sh
bash table4.sh

# Check against the expected tables
bash check.sh

echo "AIMGM" > table2_main.csv

../example1 -l 4 -N 8192 -a 0.01 2>/dev/null >> table2_main.csv
../example1 -l 5 -N 8192 -a 0.01 2>/dev/null >> table2_main.csv
../example1 -l 6 -N 8192 -a 0.01 2>/dev/null >> table2_main.csv

echo "BD-ADI" >> table2_main.csv

../example1_BDADI -l 4 -N 8192 -a 0.01 2>/dev/null >> table2_main.csv
../example1_BDADI -l 5 -N 8192 -a 0.01 2>/dev/null >> table2_main.csv
../example1_BDADI -l 6 -N 8192 -a 0.01 2>/dev/null >> table2_main.csv

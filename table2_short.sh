echo "AIMGM" > table2_short.csv

./example1 -l 2 -N 8192 -a 0.01 2>/dev/null >> table2_short.csv
./example1 -l 3 -N 8192 -a 0.01 2>/dev/null >> table2_short.csv
./example1 -l 4 -N 8192 -a 0.01 2>/dev/null >> table2_short.csv

echo "BD-ADI" >> table2_short.csv

./example1_BDADI -l 2 -N 8192 -a 0.01 2>/dev/null >> table2_short.csv
./example1_BDADI -l 3 -N 8192 -a 0.01 2>/dev/null >> table2_short.csv
./example1_BDADI -l 4 -N 8192 -a 0.01 2>/dev/null >> table2_short.csv

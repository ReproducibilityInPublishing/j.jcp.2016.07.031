echo "AIMGM" > table2.csv

./example1 -l 2 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1 -l 3 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1 -l 4 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1 -l 5 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1 -l 6 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1 -l 7 -N 8192 -a 0.01 2>/dev/null >> table2.csv

echo "BD-ADI" >> table2.csv

./example1_BDADI -l 2 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1_BDADI -l 3 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1_BDADI -l 4 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1_BDADI -l 5 -N 8192 -a 0.01 2>/dev/null >> table2.csv
./example1_BDADI -l 6 -N 8192 -a 0.01 2>/dev/null >> table2.csv
#./example1_BDADI -l 7 -N 8192 -a 0.01 2>/dev/null >> table2.csv

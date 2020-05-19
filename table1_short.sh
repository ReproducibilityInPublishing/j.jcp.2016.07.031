echo "AIMGM" > table1_short.csv

./example1 -l 9 -a 0.1666667 -N 8 2>/dev/null >> table1_short.csv
./example1 -l 9 -a 0.1666667 -N 16 2>/dev/null >> table1_short.csv

./example1 -l 9 -a 0.5 -N 32 2>/dev/null >> table1_short.csv
./example1 -l 9 -a 0.5 -N 64 2>/dev/null >> table1_short.csv

./example1 -l 9 -a 0.99 -N 100 2>/dev/null >> table1_short.csv
./example1 -l 9 -a 0.99 -N 200 2>/dev/null >> table1_short.csv

echo "BD-ADI" >> table1_short.csv

./example1_BDADI -l 9 -a 0.1666667 -N 8 2>/dev/null >> table1_short.csv
./example1_BDADI -l 9 -a 0.1666667 -N 16 2>/dev/null >> table1_short.csv

./example1_BDADI -l 9 -a 0.5 -N 32 2>/dev/null >> table1_short.csv
./example1_BDADI -l 9 -a 0.5 -N 64 2>/dev/null >> table1_short.csv

./example1_BDADI -l 9 -a 0.99 -N 100 2>/dev/null >> table1_short.csv
./example1_BDADI -l 9 -a 0.99 -N 200 2>/dev/null >> table1_short.csv

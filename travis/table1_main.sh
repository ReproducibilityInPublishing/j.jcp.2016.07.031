echo "AIMGM" > table1_main.csv

../example1 -l 9 -a 0.1666667 -N 16 2>/dev/null >> table1_main.csv
../example1 -l 9 -a 0.1666667 -N 32 2>/dev/null >> table1_main.csv

../example1 -l 9 -a 0.5 -N 64 2>/dev/null >> table1_main.csv
../example1 -l 9 -a 0.5 -N 128 2>/dev/null >> table1_main.csv

../example1 -l 9 -a 0.99 -N 200 2>/dev/null >> table1_main.csv
../example1 -l 9 -a 0.99 -N 400 2>/dev/null >> table1_main.csv

echo "BD-ADI" >> table1_main.csv

../example1_BDADI -l 9 -a 0.1666667 -N 16 2>/dev/null >> table1_main.csv
../example1_BDADI -l 9 -a 0.1666667 -N 32 2>/dev/null >> table1_main.csv

../example1_BDADI -l 9 -a 0.5 -N 64 2>/dev/null >> table1_main.csv
../example1_BDADI -l 9 -a 0.5 -N 128 2>/dev/null >> table1_main.csv

../example1_BDADI -l 9 -a 0.99 -N 200 2>/dev/null >> table1_main.csv
../example1_BDADI -l 9 -a 0.99 -N 400 2>/dev/null >> table1_main.csv

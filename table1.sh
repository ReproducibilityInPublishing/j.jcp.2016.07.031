"AIMGM" > table1.csv

./example1 -l 9 -a 0.1666667 -N 8 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.1666667 -N 16 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.1666667 -N 32 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.1666667 -N 64 2>/dev/null >> table1.csv

./example1 -l 9 -a 0.5 -N 32 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.5 -N 64 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.5 -N 128 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.5 -N 256 2>/dev/null >> table1.csv

./example1 -l 9 -a 0.99 -N 100 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.99 -N 200 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.99 -N 400 2>/dev/null >> table1.csv
./example1 -l 9 -a 0.99 -N 800 2>/dev/null >> table1.csv

"BD-ADI" >> table1.csv

./example1_BDADI -l 9 -a 0.1666667 -N 8 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.1666667 -N 16 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.1666667 -N 32 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.1666667 -N 64 2>/dev/null >> table1.csv

./example1_BDADI -l 9 -a 0.5 -N 32 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.5 -N 64 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.5 -N 128 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.5 -N 256 2>/dev/null >> table1.csv

./example1_BDADI -l 9 -a 0.99 -N 100 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.99 -N 200 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.99 -N 400 2>/dev/null >> table1.csv
./example1_BDADI -l 9 -a 0.99 -N 800 2>/dev/null >> table1.csv

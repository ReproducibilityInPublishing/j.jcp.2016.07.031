echo "AIMGM" > table1_extended.csv

../example1 -l 9 -a 0.1666667 -N 64 2>/dev/null >> table1_extended.csv

../example1 -l 9 -a 0.5 -N 256 2>/dev/null >> table1_extended.csv

../example1 -l 9 -a 0.99 -N 800 2>/dev/null >> table1_extended.csv

echo "BD-ADI" >> table1_extended.csv

../example1_BDADI -l 9 -a 0.1666667 -N 64 2>/dev/null >> table1_extended.csv

../example1_BDADI -l 9 -a 0.5 -N 256 2>/dev/null >> table1_extended.csv

../example1_BDADI -l 9 -a 0.99 -N 800 2>/dev/null >> table1_extended.csv

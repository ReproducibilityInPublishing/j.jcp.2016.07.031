echo "AIMGM" > table2_extended.csv

../example1 -l 7 -N 8192 -a 0.01 2>/dev/null >> table2_extended.csv

echo "BD-ADI" >> table2_extended.csv

../example1_BDADI -l 7 -N 8192 -a 0.01 2>/dev/null >> table2_extended.csv

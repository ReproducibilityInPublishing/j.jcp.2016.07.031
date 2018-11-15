echo "AIMGM" > table2_smoke.csv

../example1 -l 2 -N 8192 -a 0.01 2>/dev/null >> table2_smoke.csv
../example1 -l 3 -N 8192 -a 0.01 2>/dev/null >> table2_smoke.csv

echo "BD-ADI" >> table2_smoke.csv

../example1_BDADI -l 2 -N 8192 -a 0.01 2>/dev/null >> table2_smoke.csv
../example1_BDADI -l 3 -N 8192 -a 0.01 2>/dev/null >> table2_smoke.csv

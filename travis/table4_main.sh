echo "AIMGM" > table4_main.csv

../example2 -l 4 -N 15000 --alpha 0.01 2>/dev/null >> table4_main.csv
../example2 -l 5 -N 15000 --alpha 0.01 2>/dev/null >> table4_main.csv
../example2 -l 6 -N 15000 --alpha 0.01 2>/dev/null >> table4_main.csv

echo "BFSMGM" >> table4_main.csv

../example2_BFSMGM -l 4 -N 15000 --alpha 0.01 2>/dev/null >> table4_main.csv
../example2_BFSMGM -l 5 -N 15000 --alpha 0.01 2>/dev/null >> table4_main.csv
../example2_BFSMGM -l 6 -N 15000 --alpha 0.01 2>/dev/null >> table4_main.csv

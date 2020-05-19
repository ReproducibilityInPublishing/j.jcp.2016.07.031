echo "AIMGM" > table4_short.csv

./example2 -l 2 -N 15000 --alpha 0.01 2>/dev/null >> table4_short.csv
./example2 -l 3 -N 15000 --alpha 0.01 2>/dev/null >> table4_short.csv

echo "BFSMGM" >> table4_short.csv

./example2_BFSMGM -l 2 -N 15000 --alpha 0.01 2>/dev/null >> table4_short.csv
./example2_BFSMGM -l 3 -N 15000 --alpha 0.01 2>/dev/null >> table4_short.csv

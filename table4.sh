echo "AIMGM" > table4.csv

./example2 -l 2 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2 -l 3 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2 -l 4 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2 -l 5 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2 -l 6 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2 -l 7 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv

echo "BFSMGM" >> table4.csv

./example2_BFSMGM -l 2 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2_BFSMGM -l 3 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2_BFSMGM -l 4 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2_BFSMGM -l 5 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2_BFSMGM -l 6 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv
./example2_BFSMGM -l 7 -N 15000 --alpha 0.01 2>/dev/null >> table4.csv

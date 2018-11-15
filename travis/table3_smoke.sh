echo "AIMGM" > table3_smoke.csv

./example2 -l 9 --alpha 0.16666667 -N 32 2>/dev/null >> table3_smoke.csv
./example2 -l 9 --alpha 0.16666667 -N 64 2>/dev/null >> table3_smoke.csv

./example2 -l 9 --alpha 0.5 -N 50 2>/dev/null >> table3_smoke.csv
./example2 -l 9 --alpha 0.5 -N 100 2>/dev/null >> table3_smoke.csv

./example2 -l 9 --alpha 0.99 -N 50 2>/dev/null >> table3_smoke.csv
./example2 -l 9 --alpha 0.99 -N 100 2>/dev/null >> table3_smoke.csv

echo "BFSMGM" >> table3_smoke.csv

./example2_BFSMGM -l 9 --alpha 0.16666667 -N 32 2>/dev/null >> table3_smoke.csv
./example2_BFSMGM -l 9 --alpha 0.16666667 -N 64 2>/dev/null >> table3_smoke.csv

./example2_BFSMGM -l 9 --alpha 0.5 -N 50 2>/dev/null >> table3_smoke.csv
./example2_BFSMGM -l 9 --alpha 0.5 -N 100 2>/dev/null >> table3_smoke.csv

./example2_BFSMGM -l 9 --alpha 0.99 -N 50 2>/dev/null >> table3_smoke.csv
./example2_BFSMGM -l 9 --alpha 0.99 -N 100 2>/dev/null >> table3_smoke.csv

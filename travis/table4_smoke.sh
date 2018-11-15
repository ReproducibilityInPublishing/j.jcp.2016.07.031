echo "AIMGM" > table4_smoke.csv

../example2 -l 2 -N 15000 --alpha 0.01 2>/dev/null >> table4_smoke.csv
../example2 -l 3 -N 15000 --alpha 0.01 2>/dev/null >> table4_smoke.csv
../example2 -l 4 -N 15000 --alpha 0.01 2>/dev/null >> table4_smoke.csv

echo "BFSMGM" >> table4_smoke.csv

../example2_BFSMGM -l 2 -N 15000 --alpha 0.01 2>/dev/null >> table4_smoke.csv
../example2_BFSMGM -l 3 -N 15000 --alpha 0.01 2>/dev/null >> table4_smoke.csv
../example2_BFSMGM -l 4 -N 15000 --alpha 0.01 2>/dev/null >> table4_smoke.csv

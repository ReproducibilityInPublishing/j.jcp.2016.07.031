echo "AIMGM" > table4_extended.csv

../example2 -l 7 -N 15000 --alpha 0.01 2>/dev/null >> table4_extended.csv

echo "BFSMGM" >> table4_extended.csv

../example2_BFSMGM -l 7 -N 15000 --alpha 0.01 2>/dev/null >> table4_extended.csv

echo "AIMGM" > table3_extended.csv

../example2 -l 9 --alpha 0.16666667 -N 512 2>/dev/null >> table3_extended.csv

../example2 -l 9 --alpha 0.5 -N 800 2>/dev/null >> table3_extended.csv

../example2 -l 9 --alpha 0.99 -N 800 2>/dev/null >> table3_extended.csv

echo "BFSMGM" >> table3_extended.csv

../example2_BFSMGM -l 9 --alpha 0.16666667 -N 512 2>/dev/null >> table3_extended.csv

../example2_BFSMGM -l 9 --alpha 0.5 -N 800 2>/dev/null >> table3_extended.csv

../example2_BFSMGM -l 9 --alpha 0.99 -N 800 2>/dev/null >> table3_extended.csv

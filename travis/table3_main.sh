echo "AIMGM" > table3_main.csv

../example2 -l 9 --alpha 0.16666667 -N 64 2>/dev/null >> table3_main.csv
../example2 -l 9 --alpha 0.16666667 -N 128 2>/dev/null >> table3_main.csv
../example2 -l 9 --alpha 0.16666667 -N 256 2>/dev/null >> table3_main.csv

../example2 -l 9 --alpha 0.5 -N 100 2>/dev/null >> table3_main.csv
../example2 -l 9 --alpha 0.5 -N 200 2>/dev/null >> table3_main.csv
../example2 -l 9 --alpha 0.5 -N 400 2>/dev/null >> table3_main.csv

../example2 -l 9 --alpha 0.99 -N 100 2>/dev/null >> table3_main.csv
../example2 -l 9 --alpha 0.99 -N 200 2>/dev/null >> table3_main.csv
../example2 -l 9 --alpha 0.99 -N 400 2>/dev/null >> table3_main.csv

echo "BFSMGM" >> table3_main.csv

../example2_BFSMGM -l 9 --alpha 0.16666667 -N 64 2>/dev/null >> table3_main.csv
../example2_BFSMGM -l 9 --alpha 0.16666667 -N 128 2>/dev/null >> table3_main.csv
../example2_BFSMGM -l 9 --alpha 0.16666667 -N 256 2>/dev/null >> table3_main.csv

../example2_BFSMGM -l 9 --alpha 0.5 -N 100 2>/dev/null >> table3_main.csv
../example2_BFSMGM -l 9 --alpha 0.5 -N 200 2>/dev/null >> table3_main.csv
../example2_BFSMGM -l 9 --alpha 0.5 -N 400 2>/dev/null >> table3_main.csv

../example2_BFSMGM -l 9 --alpha 0.99 -N 100 2>/dev/null >> table3_main.csv
../example2_BFSMGM -l 9 --alpha 0.99 -N 200 2>/dev/null >> table3_main.csv
../example2_BFSMGM -l 9 --alpha 0.99 -N 400 2>/dev/null >> table3_main.csv

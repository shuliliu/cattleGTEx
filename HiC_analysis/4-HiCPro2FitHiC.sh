#!/bin/bash


for bin in 20000 40000 150000 500000 1000000
do
python ~u/bin/fithic/fithic/utils/HiCPro2FitHiC.py -i ./hic_results/matrix/Sample1/raw/${bin}/Sample1_${bin}.matrix \
-b ./hic_results/matrix/Sample1/raw/${bin}/Sample1_${bin}_abs.bed \
-s ./hic_results/matrix/Sample1/iced/${bin}/Sample1_${bin}_iced.matrix.biases \
-o ./hic_results/matrix/Sample1/iced/${bin} \
-r ${bin} 
done



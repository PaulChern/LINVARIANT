#!/usr/bin/sh

dir0=`pwd`
Lx=41; Ly=41; Lz=2; NT=10;
for type in random 1 2 3 4 5 6
do
cp mma.run mma.${type}.run
sed -i "s/case/${type}/g" mma.${type}.run
sed -i "s/Lx/${Lx}/g" mma.${type}.run
sed -i "s/Ly/${Ly}/g" mma.${type}.run
sed -i "s/Lz/${Lz}/g" mma.${type}.run
sed -i "s/NT/${NT}/g" mma.${type}.run
sbatch mma.${type}.run
done

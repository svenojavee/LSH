
logdir='/mnt/beegfs/robingrp/sojavee/LSHResults/logs'

for h2 in 0.15 0.25 0.35 0.45  0.55 0.65
do
for K in 0.01 0.02 0.05
do 
for P in 0.9 0.95 1
do 
for n in {1..100}
do

output=s${n}_M10k_N20k_h${h2}_K${K}_P${P}

echo '#!/bin/bash' > $logdir/$output.sh
echo "
#SBATCH --job-name=s${n}_M10k_N20k_h${h2}_K${K}_P${P}
#SBATCH --partition=defaultp
#SBATCH --reservation=robingrp_106
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 3
#SBATCH --mem=15G
#SBATCH --time 0-01:30:00
#SBATCH --output=$logdir/${output}.log

module add R

cmd=\"/nfs/scistore13/robingrp/sojavee/gcta64 --grm /mnt/beegfs/robingrp/sojavee/LSHResults/simM10k_N20k --pheno /mnt/beegfs/robingrp/sojavee/LSHResults/${output}.phen  --reml --threads 3 --out /mnt/beegfs/robingrp/sojavee/LSHResults/results/${output} \"
echo
echo \$cmd
echo 
\$cmd" >> $logdir/$output.sh
sbatch $logdir/$output.sh



#/nfs/scistore13/robingrp/sojavee/gcta64 --grm /mnt/beegfs/robingrp/sojavee/LSHResults/simM10k_N20k --pheno /mnt/beegfs/robingrp/sojavee/LSHResults/s${n}_M10k_N20k_h${h2}_K${K}_P${P}.phen  --reml --threads 2 --out /mnt/beegfs/robingrp/sojavee/LSHResults/results/s${n}_M10k_N20k_h${h2}_K${K}_P${P}

#/nfs/scistore13/robingrp/sojavee/gcta64 --grm /mnt/beegfs/robingrp/sojavee/LSHResults/simM10k_N20k --pheno /mnt/beegfs/robingrp/sojavee/LSHResults/s${n}_M10k_N20k_h${h2}_K${K}_P${P}.phen  --reml --threads 2 --out /mnt/beegfs/robingrp/sojavee/LSHResults/results/s${n}_M10k_N20k_h${h2}_K${K}_P${P}


done
done
done
done
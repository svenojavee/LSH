
logdir='/mnt/beegfs/robingrp/sojavee/LSHResults/PCGClogs'

for h2 in 0.15 0.25 0.35 0.45 0.55 0.65
do
for K in 0.005 0.01 0.02  # 0.05
do 
for P in 0.25 # 0.5 0.75 0.9 1 #1.25 1.5 2 4 
do 
for n in  {1..100}
do

output=calc_s${n}_M10k_N20k_h${h2}_K${K}_P${P}

echo '#!/bin/bash' > $logdir/$output.sh
echo "
#SBATCH --job-name=s${n}_M10k_N20k_h${h2}_K${K}_P${P}
#SBATCH --partition=defaultp
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=5G
#SBATCH --time 0-00:30:00
#SBATCH --output=$logdir/${output}.log

source activate /nfs/scistore13/robingrp/sojavee/.conda/envs/pcgc

cmd=\"python /nfs/scistore13/robingrp/sojavee/S-PCGC/pcgc_main.py --sumstats /mnt/beegfs/robingrp/sojavee/LSHResults/sumStatsPCGC/s${n}_M10k_N20k_h${h2}_K${K}_P${P}. --prodr2 /mnt/beegfs/robingrp/sojavee/LSHResults/sumStatsPCGC/simM10k_N20k. --M 10000 --out /mnt/beegfs/robingrp/sojavee/LSHResults/sumStatsPCGC/s${n}_M10k_N20k_h${h2}_K${K}_P${P} \"

echo
echo \$cmd
echo 
\$cmd" >> $logdir/$output.sh
sbatch $logdir/$output.sh


done
done
done
done



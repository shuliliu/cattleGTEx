#!/bin/bash



prog_dir="~"
##########################
use_slurm=1 # use slurm (1: use, 0: don't use)
slurm_partition=short
source $prog_dir/slurm_util_sbams.sh

run_script()
{
    # if slurm is not used
  # this function expects only one argument: 
  # the script file.
  
  # if slurm is used,
  # this function expects 6 arguments: 
  # 1) script file, 2) partition, 3) no. of nodes, 
  # 4) no. of tasks, 5) time to run, 6) memory
  
  if [ $# -lt 1 ]; then 
    echo "no script to run"; 
    return 1; 
  fi
  script_fn=$1
  
  if [ $use_slurm -eq 1 ]; then 
    submit_slurm_job $@
  else
    sh $script_fn
  fi
  
  return 0

}

###info
tissue=$1
dap="~/bin/dap-1.0.0/dap-master/dap_src/dap-g"
aim_folder="~/eQTL/vcf_tissues"
###################################
cd ${aim_folder}/${tissue}
split -l 9600 ${tissue}.assemble.cmd ${tissue}.cmd
for cmd in `ls ${tissue}.cmd*`
do
cd ${prog_dir}
mkdir log_${cmd}
script_fn=$prog_dir/generate_sbams_dap-g.sh
echo -e "cd ${aim_folder}\ncd ${tissue}\n \
        cat ${cmd} | head -n \$SLURM_ARRAY_TASK_ID |tail -n 1 > script.tmp.\$SLURM_ARRAY_TASK_ID.txt\n\
        sh script.tmp.\$SLURM_ARRAY_TASK_ID.txt\n\
        sbams=\`cat script.tmp.\$SLURM_ARRAY_TASK_ID.txt | awk '{print \$8}'\`\n\
        gene=\`echo \${sbams} | awk -F "[/.]" '{print \$2}'\`\n\
        ${dap} -d \${sbams} -p ./${tissue}/prior/\${gene}.prior -ld_control 0.5 --all -t 2 > ./dap_output/\${gene}.dap.tabix\n\
        rm script.tmp.\$SLURM_ARRAY_TASK_ID.txt\n" > "$script_fn"

Num_gene=`wc -l ${aim_folder}/${tissue}/${cmd} | awk '{print $1}'`
times=`echo $((Num_gene/600))`
last=`echo $((Num_gene%600))`
for i in $(seq 1 $times) ##or for ((i=1;i<=${times};i++))
do
array_1=`echo $((600*(i-1)+1))`
array_2=`echo $((600*i))`
echo ${array_1} 
echo ${array_2}
run_script ${script_fn} $slurm_partition 1 1 "1:0:0" "4GB" ${array_1}-${array_2} ${cmd}

##check if it's out of submission limitation; check whether the command has been submitted or not; if not, wait for 15m, and submit it again.
cd ${prog_dir}
tail=`tail -n 1 nominal_muscle.out | cut -d " " -f 1`
submitted="Submitted"
if [ ${tail} != ${submitted} ]; then
sleep 20m;
run_script ${script_fn} $slurm_partition 1 1 "1:0:0" "4GB" ${array_1}-${array_2} ${cmd}
fi

sleep 8m;
done


if [ ${last} -gt 0 ]; then
array_1=`echo $((600*times+1))`
array_2=`echo $((600*times+last))`
echo ${array_1}
echo ${array_2}
run_script ${script_fn} $slurm_partition 1 1 "1:0:0" "4GB" ${array_1}-${array_2} ${cmd}
sleep 8m;
fi
done

echo "DONE"




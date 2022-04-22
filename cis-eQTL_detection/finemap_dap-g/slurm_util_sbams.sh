get_slurm_header()
{
  if [ $# -lt 10 ]; then 
         echo "#not enough arguments in ${FUNCNAME}."; 
       return 1; 
       fi
       
       partition=$1
       workdir=$2
       time=$3
       nodes=$4
       ntasks=$5
       mem=$6
       job_name=$7
       out_fn=$8
       err_fn=$9
       array=${10}


       
       cmd_header="#"'!'"/bin/sh"
       cmd_header="$cmd_header\n#SBATCH --partition=${partition}"
       cmd_header="$cmd_header\n#SBATCH --chdir=${workdir}"
       cmd_header="$cmd_header\n#SBATCH --time=${time}"
       cmd_header="$cmd_header\n#SBATCH --nodes=${nodes}"
       cmd_header="$cmd_header\n#SBATCH --ntasks=${ntasks}"
       cmd_header="$cmd_header\n#SBATCH --mem=${mem}"
       cmd_header="$cmd_header\n#SBATCH --job-name=${job_name}"
       cmd_header="$cmd_header\n#SBATCH --output=${out_fn}"
       cmd_header="$cmd_header\n#SBATCH --error=${err_fn}"
       cmd_header="$cmd_header\n#SBATCH --array=${array}"
       cmd_header="$cmd_header\n"
       
       echo $cmd_header
       return 0
}

submit_slurm_job()
{
  if [ $# -lt 1 ]; then 
    echo "No script to submit"; 
    return 1; 
  fi
  script_fn=$1
  
  # default settings
  partition="debug"
  nodes=1
  ntasks=1
  time="1:0:0"
  mem="3GB"

  if [ $# -ge 2 ]; then partition=$2; fi
  if [ $# -ge 3 ]; then nodes=$3; fi
  if [ $# -ge 4 ]; then ntasks=$4; fi
  if [ $# -ge 5 ]; then time=$5; fi
  if [ $# -ge 6 ]; then mem=$6; fi
  if [ $# -ge 7 ]; then array=$7; fi
  if [ $# -ge 8 ]; then cmd=$8; fi
  # slurm header
  script_dir="$(dirname "$script_fn")"
  job_name="$(basename "$script_fn")"
  slurm_job_fn="${script_fn}.slurm"
  cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "./log_${cmd}/sbams_%%A_%%a.out" "./log_${cmd}/sbams_%%A_%%a.err" ${array})

  # slurm modules
  cmd_modules="module load gsl"
  cmd_modules="$cmd_modules\nmodule load r"
  cmd_modules="$cmd_modules\nmodule load bcftools"
  cmd_modules="$cmd_modules\nmodule load zlib"
  cmd_modules="$cmd_modules\nmodule load boost"
  cmd_modules="$cmd_modules\nmodule load tabix"
  cmd_modules="$cmd_modules\n"
  
  # target script
  cmd_body="sh \"$script_fn\""
  cmd_body="$cmd_body\necho DONE"

  # save slurm job
  printf "$cmd_header\n" > $slurm_job_fn
  printf "$cmd_modules\n" >> $slurm_job_fn
  printf "$cmd_body\n" >> $slurm_job_fn

  # submit slurm job
  sbatch $slurm_job_fn
  
  return 0
}


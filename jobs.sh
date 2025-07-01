#!/usr/bin/bash

maxJ=1.2
minJ=0.0
NJ=20

maxD=3.5
minD=0.0
ND=50

maxB=1.5
minB=1.5
NB=1

outputDir="./output-dir"
outfilesDir="./outfiles" # Teh directory where stdout will be directed to
maxTime="10:00:00"
maxConcurrentJobs=250
###########################################################################################

module load python openmpi fftw
N=$((NJ*ND*NB - 1)) # Indiced through 0 to N-1

mkdir -p $outfilesDir

script=$(cat <<EOF
#!/usr/bin/bash
#SBATCH -p shared
#SBATCH -n 1 
#SBATCH -c 1 
#SBATCH -J parameter-scan
#SBATCH -t ${maxTime}
#SBATCH -o ${outputfilesDir}/\${SLURM_ARRAY_TASK_ID} 
NJ=${NJ}
ND=${ND}
NB=${NB}
N=${N}
index=\${SLURM_ARRAY_TASK_ID}

i=\$(python -c "import math; x = int(math.floor(\${index}/(\${NB}*\${ND}))); print(x)")
rem=\$(python -c "x = int(\${index}) % int(\${NB}*\${ND}); print(x)")
j=\$(python -c "import math; x = int(math.floor(\${rem} / int(\${NB}))); print(x)")
k=\$(python -c "x = int(\${rem}) % int(\${NB}); print(x)")

Jval=\$(python -c "x = (\${i}/\${NJ})*(${maxJ} - ${minJ}) + ${minJ}; print(x)")
Dval=\$(python -c "x = (\${j}/\${ND})*(${maxD} - ${minD}) + ${minD}; print(x)")
Bval=\$(python -c "x = (\${k}/\${NB})*(${maxB} - ${minB}) + ${minB}; print(x)")

mpirun -np 1 ./bin/PT \${Jval} \${Jval} 1 \${Dval} \${Dval} 1 \${Bval} \${Bval} 1 ${outputDir} 
EOF
)

echo "${script}" | sbatch --array=0-${N}%${maxConcurrentJobs}  


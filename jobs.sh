#!/usr/bin/bash

maxJ=3.0
minJ=0.0
NJ=2

maxD=3.0
minD=0.0
ND=2

maxB=1.5
minB=1.0
NB=1

outputDir="./output-dir"
maxTime="00:15:00"
maxConcurrentJobs=4
###########################################################################################

module load python openmpi
N=$((NJ*ND*NB - 1)) # Indiced through 0 to N-1




script=$(cat <<EOF
#!/usr/bin/bash
#SBATCH -p test
#SBATCH -n 1 
#SBATCH -c 1 
#SBATCH -J parameter-scan
#SBATCH -t ${maxTime}

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


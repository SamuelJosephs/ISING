import numpy as np 
import subprocess 
J_arr = np.linspace(-5,5,10)
D_arr = np.linspace(-5,5,10)
B = 1.5
time_limit = "00:45:00"
num_threads=10
submit = False 
for J in J_arr:
    for D in D_arr:
        script = f"""#!/bin/bash
#SBATCH --job-name={J:.3f}_{D:.3f}_{B:.3f}          # Job name
#SBATCH --time={time_limit}             # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                  # Number of tasks (1 task = 1 process)
#SBATCH --cpus-per-task={num_threads}           # Number of CPU cores per task (OpenMP threads)
#SBATCH --output=./outfiles/omp_job_{J:.3f}_{D:.3f}_{B:.3f}.out        # Standard output
#SBATCH --error=./errfiles/omp_job_{J:.3f}_{D:.3f}_{B:.3f}.err         # Standard error
#SBATCH -p shared
# Set the number of OpenMP threads
export OMP_NUM_THREADS={num_threads}
./bin/ISING {J} {D} 1.5 ./output_dir_{J:.3f}_{D:.3f}_{B:.3f} && python3 skyrmion_evolution_simple.py ./output_dir_{J:.3f}_{D:.3f}_{B:.3f} ./visualisations _{J:.3f}_{D:.3f}_{B:.3f}_
        """
        with open(f"script_{J:.3f}_{D:.3f}_{B:.3f}.sh", "w") as f:
            f.write(script)
        if submit:
            result = subprocess.run([f"sbatch",f"script_{J:.3f}_{D:.3f}_{B:.3f}.sh"])

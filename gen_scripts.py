import numpy as np 
import subprocess 
J_arr = np.linspace(-5,5,10)
D_arr = np.linspace(-5,5,10)
J_prime_arr = np.array([0.0])
D_prime_array = np.array([0.0])
B = 1.5
time_limit = "00:45:00"
num_threads=10
submit = False 



for J_prime in J_prime_arr:
    for D_prime in D_prime_array: 
        for J in J_arr:
            for D in D_arr:

                identifier = f"_J={J:.2f}_J'={J_prime:.2f}_D={D:.2f}_D'={D_prime:.2f}_B={B:.2f}_"
                script = f"""#!/bin/bash
        #SBATCH --job-name={identifier}          # Job name
        #SBATCH --time={time_limit}             # Time limit (hh:mm:ss)
        #SBATCH --ntasks=1                  # Number of tasks (1 task = 1 process)
        #SBATCH --cpus-per-task={num_threads}           # Number of CPU cores per task (OpenMP threads)
        #SBATCH --output=./outfiles/omp_job_{identifier}.out        # Standard output
        #SBATCH --error=./errfiles/omp_job_{identifier}.err         # Standard error
        #SBATCH -p shared
        # Set the number of OpenMP threads
        export OMP_NUM_THREADS={num_threads}
        ./bin/ISING {J} {J_prime} {D} {D_prime} 1.5 ./output_dir_{identifier} && python3 skyrmion_evolution_simple.py ./output_dir_{identifier} ./visualisations {identifier}
                """
                with open(f"script_{identifier}.sh", "w") as f:
                    f.write(script)
                if submit:
                    result = subprocess.run([f"sbatch",f"script_{identifier}.sh"])

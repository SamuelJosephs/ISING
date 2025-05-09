
import os
import numpy as np
import subprocess
import textwrap

# Define the directories that need to be created
directories = ['outfiles', 'errfiles', 'output_dir', 'visualisations']

# Create each directory if it does not already exist
for d in directories:
    os.makedirs(d, exist_ok=True)

J_arr = np.linspace(-3, 3, 10)
D_arr = np.linspace(-3, 3, 10)
J_prime_arr = np.array([0.0])
D_prime_array = np.array([0.0])
B_arr = np.array([1.5])
time_limit = "04:00:00"
num_threads = 1
submit = True

VISSCRIPT = 'skyrmion_evolution_simple.py'
ISING = 'bin/ISING'

for B in B_arr:
    for J_prime in J_prime_arr:
        for D_prime in D_prime_array:
            for J in J_arr:
                for D in D_arr:
                    identifier = f"J_{J:.2f}-Jpr_{J_prime:.2f}-D_{D:.2f}-Dpr_{D_prime:.2f}-B_{B:.2f}"
                    script = textwrap.dedent(f"""\
                    #!/bin/bash
                    #SBATCH --job-name={identifier}          # Job name
                    #SBATCH --time={time_limit}             # Time limit (hh:mm:ss)
                    #SBATCH --ntasks=1                  # Number of tasks (1 task = 1 process)
                    #SBATCH --cpus-per-task={num_threads}           # Number of CPU cores per task (OpenMP threads)
                    #SBATCH --output=./outfiles/job_{identifier}.out        # Standard output
                    #SBATCH --error=./errfiles/job_{identifier}.err         # Standard error
                    #SBATCH -p shared

                    module load ffmpeg
                    # Set the number of OpenMP threads
                    export OMP_NUM_THREADS={num_threads}
                    {ISING} {J} {J_prime} {D} {D_prime} 1.5 ./od_{identifier} && python3 {VISSCRIPT} ./od_{identifier} ./visualisations {identifier}
                    """)
                    script_filename = f"script_{identifier}.sh"
                    with open(script_filename, "w") as f:
                        f.write(script)
                    if submit:
                        subprocess.run(["sbatch", script_filename])

#!/bin/bash


NJ=2
ND=2
NB=2
N=$((NJ*ND*NB - 1)) # Indiced through 0 to N-1

script=$(cat <<EOF

NJ=$(NJ)
ND=$(ND)
NB=$(NB)
N=$(N)
index=$(SLURM_ARRAY_TASK_ID)


EOF
)

#!/bin/bash

# List available groups
h5ls optimization.h5

echo "First pre-optimization"
pysisplot --opt --h5_group first_pre

echo "Last pre-optimization"
pysisplot --opt --h5_group last_pre

echo "NEB optimization, energies"
pysisplot --cosens --h5_group opt

echo "NEB optimization, forces"
pysisplot --cosforces --h5_group opt

echo "TS optimization"
pysisplot --opt --h5_group tsopt

echo "IRC"
pysisplot --irc

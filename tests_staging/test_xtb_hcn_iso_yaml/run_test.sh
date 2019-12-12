#!/bin/bash

yaml_inp=test_xtb_hcn_iso.yaml

function set_max_cycles {
    sed -i -r "s/max_cycles: [0-9]+/max_cycles: $1/" $yaml_inp
}

echo "yes" | pysisrun --clean

function run_restart {
    set_max_cycles 5
    pysisrun test_xtb_hcn_iso.yaml | tee 01_first5.out
    set_max_cycles 10
    pysisrun test_xtb_hcn_iso.yaml --restart | tee 02_second5.out
}

function run_full {
    set_max_cycles 50
    pysisrun test_xtb_hcn_iso.yaml | tee 001_full.out
}

#run_full
#run_restart

$1

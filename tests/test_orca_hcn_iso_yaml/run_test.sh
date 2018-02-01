#!/bin/bash

yaml_inp=test_orca_hcn_iso_yaml.yaml

function set_max_cycles {
    sed -i -r "s/max_cycles: [0-9]+/max_cycles: $1/" $yaml_inp
}

function run_restart {
    echo "yes" | pysisrun --clean
    set_max_cycles 5
    pysisrun $yaml_inp | tee 01_first5.out
    set_max_cycles 20
    pysisrun $yaml_inp --restart | tee 02_second5.out
}

function only_restart {
    set_max_cycles 9
    pysisrun $yaml_inp --restart | tee 02_second5.out
}

function run_full {
    set_max_cycles 20
    pysisrun $yaml_inp | tee 00_full.out
}

#run_full
#run_restart

$1

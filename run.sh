#!/bin/bash
clear
export PYTHONPATH=`pwd`:$PYTHONPATH
rm cycle*.trj
python tests/anapot.py
#module load orca/4.0.0.2
#python tests/h2o.py
#python tests/hcn_iso/hcn_iso.py
#python tests/mullerbrownpot.py
#python tests/i18i19/i18i19_pt.py

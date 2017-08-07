#!/bin/bash
clear
export PYTHONPATH=`pwd`:$PYTHONPATH
rm *.trj
#python tests/test_anapot.py
#python tests/test_mullerbrownpot.py
#module load orca/4.0.0.2
#python tests/h2o.py
#python tests/hcn_iso/hcn_iso.py
#python tests/i18i19/i18i19_pt.py
python irc/IMK.py

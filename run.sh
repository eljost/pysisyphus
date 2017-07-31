#!/bin/bash
clear
export PYTHONPATH=`pwd`:$PYTHONPATH
#python tests/anapot.py
#module load orca/4.0.0.2
python tests/h2o.py
#python tests/mullerbrownpot.py

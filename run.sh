#!/bin/bash
clear
export PYTHONPATH=`pwd`:$PYTHONPATH
python tests/anapot.py
#python tests/h2o.py
python tests/mullerbrownpot.py

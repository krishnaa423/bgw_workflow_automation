#!/bin/bash

# The ampersand is at the end so it does not block the shell.
# stderr (2) is directed to stdout(1), which is directed to log file. 
python run_workflow.py &> run_workflow.log  &
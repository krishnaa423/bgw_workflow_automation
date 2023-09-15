# bgw_workflow_automation

Can use the following scripts to manage workflow lifecycle. 

- *python create_workflow.py*
- *python run_workflow.py* or *./script_run_workflow.sh* to launch a non-blocking process that manages the workflow.  
- *remove_workflow.py*

Need to add a CIF file to the project folder and reference at in *create_workflow.py*. Additionally can manually set the struct object in the code. All workflow parameters are set in the *input* dictionary at the start of *create_workflow.py*.  

Take a look at the templates to see examples on how to modify the workflow scripts for your system. The following folders are created for now. 

- 1.1-scf
- 1.2-ph
- 1.3-epw
- 1.5-bands
- 1.6-dos
- 1.7-pdos
- 1.8-pp
- 2.1-wfn
- 2.2-wfnq
- 3.1-wfn_fi
- 3.2-wfnq_fi
- 4-epsilon
- 5.1-sigma
- 5.2-inteqp
- 6-kernel
- 7-absorption
- 8-esf

You can always add or remove folders by copying code templates for other folders. Usually you would need link scripts, input files, jobscripts, and optional utilities for a folder. 

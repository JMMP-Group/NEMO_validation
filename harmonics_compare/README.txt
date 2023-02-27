
Directory containing script which perform validation of harmonic data. Uses pregenerated harmonic obs and harmonic model output.

Overview:



Workflow:

Select machine in config.sh
Ensure paths are correctly set in <machine>_config.sh

Execute with
. ./iter_preprocess_harm.sh

This reads <machine>_config.sh to set paths
Then queues monthly jobs using <machine>_pre_process_harm.sh
 Which runs the machine independent script: pre_process_harm.py
 which uses the obs_operator to extract tidegauge objects for model points nearest observations. These are then saved
 to file.




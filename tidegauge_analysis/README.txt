Directory containing script which perform validation against hourly ssh data.

Workflow:

Pitch machine in config.sh
Ensure paths are correctly set in <machine>_config.sh

Execute with
. ./iter_tg_analysis.sh

This reads <machine>_config.sh to set paths
Then queues a job using <machine>_rs_hourly_ssh.pbs
Which runs the machine independent script: rs_hourly_ssh.py with methods in validate_ssh_tg_hourly.py

Figures are stored in FIGS. Logs are stored in LOGS.


Other files:
============
analyse_ssh_hourly.py
	A partial clone and maybe parent of validate_ssh_tg_hourly.py
	Should be removed when it is thoroughly borrowed from

harmonic_variation.py
	WIP: computing error estimates on harmonic amplitudes generated from observational data of unknown nodal phase.


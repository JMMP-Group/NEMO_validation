
Directory containing script which perform validation against hourly ssh data.

Overview:
There are 3 elements, extract (SSH from model at obs locations); analyse (concurrent ssh and threshld analysis 
on model and obs data); plot. For short timeseries these can and where done as a single (slow) job. However with 10 years of data it makes sense to slice the stages to spread them across parallel tasks.
The extract stage can be sliced across monthly jobs
The analysis stage can be sliced across the port dimension.
This slices makes the entire process fast on JASMIN, though a little harder to follow.


Workflow:

Set machine in config.sh
Ensure paths are correctly set in <machine>_config.sh

Execute with
. ./iter_tg_preprocess.sh

This reads <machine>_config.sh to set paths
Then queues monthly jobs using <machine>_pre_process_tg_monthly.sh
 Which runs the machine independent script: pre_process_tg_monthly.py with methods in validate_ssh_tg_hourly.py
 which creates monthly files of extracted NEMO data at gauge locations

Execute
. ./iter_tg_analysis.sh

This queues N_PORTS jobs <machine>_process_tg_monthly.sh
which extracts the port of interest from all the monthly data.
Harmonic and threshold analyses are applied and save in one file per port.


Then in the terminal that executed the above:

python postprocess_tg_monthly.py
python plot_tg.py

Concat all the port data into a single file
then plots outputs into figure directory: DN_OUT



Might need to make some directories if they are not there. E.g.:
mkdir /gws/nopw/j04/jmmp/CO9_AMM15_validation/P1.5c/
mkdir /gws/nopw/j04/jmmp/CO9_AMM15_validation/P1.5c/tg_analysis/
mkdir /gws/nopw/j04/jmmp/CO9_AMM15_validation/P1.5c/tg_analysis/FIGS/

Other files:
============
analyse_ssh_hourly.py
	A partial clone and maybe parent of validate_ssh_tg_hourly.py
	Should be removed when it is thoroughly borrowed from

harmonic_variation.py
	WIP: computing error estimates on harmonic amplitudes generated from observational data of unknown nodal phase.


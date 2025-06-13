"""
Config file for python, to store directory paths, variable choices etc that are common to multiple scripts.
Reads data set in config.sh
 Enables work flow to be machine independent


"""

from socket import gethostname
from os import environ

class config:
    """
    from config import config
    c = config()  # initialise
    thing=c.THING
    """
    def __init__(self):
        # read SHELL variables. Set in config.sh
        #self.machine = environ.get('MACHINE')
        self.coast_repo = self.get_shell_var('COAST_REPO', True)
        self.dn_fes = self.get_shell_var('DN_FES', True)
        self.dn_obs      = self.get_shell_var('DN_OBS', True)
        self.fn_out_dir      = self.get_shell_var('FN_OUT_DIR', True)
        #self.dn_out = self.get_shell_var('DN_OUT', True)
        self.fn_nemo_data = self.get_shell_var('FN_NEMO_DATA', True)
        self.fn_nemo_domain  = self.get_shell_var('FN_NEMO_DOMAIN', True)
        #self.fn_analysis_out      = self.get_shell_var('FN_ANALYSIS_OUT', True)

        self.grid_obs_rad  = self.get_shell_var('GRID_OBS_RAD', True, int)
        self.thin_obs_rad  = self.get_shell_var('THIN_OBS_RAD', True, int)

    def get_shell_var(self, var:str, debug=False, as_type=str):
        try:
            val = environ.get(var.upper())
            if debug: print(f"{var}: {val}")
            return as_type(val)
        except:
            print(f"Problem with getting shell variable: {var}")

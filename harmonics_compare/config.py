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
        self.run_name = self.get_shell_var('RUN_NAME', True)
        self.dn_fes = self.get_shell_var('DN_FES', True)
        #self.fn_fes_amp = self.get_shell_var('FN_FES_AMP', True)
        #self.fn_fes_pha = self.get_shell_var('FN_FES_PHA', True)
        self.fn_harm_obs = self.get_shell_var('FN_HARM_OBS', True)
        self.dn_out = self.get_shell_var('DN_OUT', True)
        self.fn_nemo_data = self.get_shell_var('FN_NEMO_DATA', True)
        self.fn_nemo_domain  = self.get_shell_var('FN_NEMO_DOMAIN', True)
        self.fn_nemo_cfg  = self.get_shell_var('FN_NEMO_CFG', True)
        self.fn_analysis_out      = self.get_shell_var('FN_ANALYSIS_OUT', True)

    def get_shell_var(self, var:str, debug=False):
        try:
            if debug: print(f"{var}: {environ.get(var.upper())}")
            return environ.get(var.upper())
        except:
            print(f"Problem with getting shell variable: {var}")

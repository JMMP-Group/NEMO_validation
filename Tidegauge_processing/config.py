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
        self.mass     = self.get_shell_var('MASS', True)
        self.fn_nemo_data = self.get_shell_var('FN_NEMO_DATA', True)
        self.fn_nemo_domain  = self.get_shell_var('FN_NEMO_DOMAIN', True)
        self.fn_nemo_cfg  = self.get_shell_var('FN_NEMO_CFG', True)
        self.fn_out = self.get_shell_var('FN_OUT', True)
        self.fn_obs = self.get_shell_var('FN_OBS', True)
        self.dn_out = self.get_shell_var('DN_OUT', True)
        self.fn_extr_out         = self.get_shell_var('FN_EXTR_OUT', True)
        self.fn_analyse_out      = self.get_shell_var('FN_ANALYSE_OUT', True)
        self.fn_ssh_hourly_stats = self.get_shell_var('FN_SSH_HOURLY_STATS', True)

    def get_shell_var(self, var:str, debug=False):
        try:
            if debug: print(f"{var}: {environ.get(var.upper())}")
            return environ.get(var.upper())
        except:
            print(f"Problem with getting shell variable: {var}")

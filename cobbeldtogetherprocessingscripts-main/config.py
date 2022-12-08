"""
Config file for python, to store directory paths, variable choices etc that are common to multiple scripts.
Reads data set in config.sh
 Enables work flow to be machine independent

from config import config
thing=config.THING
where:
    less config.py
    THING = "4b6...5ea"

"""

from socket import gethostname
from os import environ

class config:
    """
    # main.py
    from config import config

    print(config.coast_repo)
    """
    def get_shell_var(self, var:str, debug=False):
        try:
            if debug: print(f"{var}: {environ.get(var.upper())}")
            return environ.get(var.upper())
        except:
            print(f"Problem with getting shell variable: {var}")

    # read SHELL variables. Set in config.sh
    #machine = environ.get('MACHINE')
    region = self.get_shell_var('REGION', True)
    machine = self.get_shell_var('MACHINE', True)
    coast_repo = self.get_shell_var('COAST_REPO', True)
    din_en4 = self.get_shell_var('DIN_EN4', True)
    dout_en4 = self.get_shell_var('DOUT_EN4', True)
    fn_cfg_prof = self.get_shell_var('FN_CFG_PROF', True)


class bounds:
    """ bounds class to define cartesian box for preprocesses EN4 data"""
    def get_bounds(self, region):
        if region == "AMM15":
            lonbounds=[-25.47, 16.25]
            latbounds=[43, 64.5]
            return lonbounds, latbounds

        else:
            print(f"Not expecting region {region}")
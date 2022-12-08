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
        self.region      = self.get_shell_var('REGION', True)
        self.machine     = self.get_shell_var('MACHINE', True)
        self.coast_repo  = self.get_shell_var('COAST_REPO', True)
        self.din_en4     = self.get_shell_var('DIN_EN4', True)
        self.dout_en4    = self.get_shell_var('DOUT_EN4', True)
        self.fn_cfg_prof = self.get_shell_var('FN_CFG_PROF', True)

    def get_shell_var(self, var:str, debug=False):
        try:
            if debug: print(f"{var}: {environ.get(var.upper())}")
            return environ.get(var.upper())
        except:
            print(f"Problem with getting shell variable: {var}")



class bounds:
    """ bounds class to define cartesian box for preprocesses EN4 data"""
    def __init__(self, region):
        self.lonbounds, self.latbounds = self.get_bounds(region)


    def get_bounds(self, region):
        if region == "AMM15":
            lonbounds=[-25.47, 16.25]
            latbounds=[43, 64.5]
            return lonbounds, latbounds

        else:
            print(f"Not expecting region {region}")

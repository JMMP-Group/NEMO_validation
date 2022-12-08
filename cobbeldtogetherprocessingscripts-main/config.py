"""
config file to store directory paths, common to multiple scripts. Is user/machine dependent

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
    machine = environ.get('MACHINE')  # read SHELL variable. Set in config.sh
    print(f"machine:{machine}")
    if "LOTUS" in machine.upper():
        # location of COAsT repo, if using a particular branch
        coast_repo = "/home/users/jelt/GitHub/COAsT"
        din_en4 = "/home/users/jelt/EN4/"
        dout_en4 = "/home/users/jelt/tmp/"
        fn_cfg_prof = "/home/users/jelt/GitHub/COAsT/config/example_en4_profiles.json"

    elif "SPICE" in machine.upper():
        # location of COAsT repo, if using a particular branch
        coast_repo = "/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY"
        din_en4 = "/scratch/fred/EN4/"
        dout_en4 = "/scratch/fred/EN4/"
        fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"

    else:
        print(f"Do not recognise machine: {machine}")



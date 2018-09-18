import pandas as pd

from .core import ShootingPoint


def log(sp_name, CVs, result, filename='log'):
    """Logs information after a shooting point run to log file.

    Parameters
    ----------
    sp_name : str
    CVs : str
        Formatted string containing the values of all CVs.
    result : str
        Outcome of shooting point evaluation. Options are 'accept',
        'reject', or 'inconclusive'."""

    with open(filename, 'a') as f:
        line = "{}\t{}\t{}".format(sp_name, CVs, result)
        f.write(line)
    return


def get_colvar_header(colvar_file):
    """Grabs the names of the collective variables from COLVAR file.
    
    Parameters
    ----------
    colvar_file : str
    
    Returns
    -------
    header : list of str
        List containing all the names of the collective variables."""
    with open(colvar_file, 'r') as f:
        header = f.readline()
    header = header.split(' ')[3:]
    return header

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
    with open(colvar_file, 'r') as f:
        header = f.readline()

    header = header.split(' ')[3:]
    return header

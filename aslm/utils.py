

def log_run(sp_name, CVs, result, filename='log'):
    """Logs information after a shooting point run to log file.

    Parameters
    ----------
    sp_name : str
    CVs : list of int or float
        Formatted string containing the values of all CVs.
    result : str
        Outcome of shooting point evaluation. Options are 'accept',
        'reject', or 'inconclusive'.
    filename : str, optional
        Name of log file."""

    CVs = [str(n) for n in CVs]
    with open(filename, 'a') as f:
        line = "{} {} {}\n".format(sp_name, ' '.join(CVs), result)
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
    header[-1] = header[-1].rstrip()
    return header


def log_header(colvar_file, filename='log'):
    colvar_header = get_colvar_header(colvar_file)
    with open(filename, 'a') as f:
        f.write("{} {} {}\n".format('NAME', ' '.join(colvar_header), 'RESULT'))
    return

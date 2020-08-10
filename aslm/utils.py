import numpy as np


def log_run(sp_name, CVs, result, forward_commit_basin, filename='log'):
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
        line = "{} {} {} {}\n".format(sp_name, ' '.join(CVs), result,
                                      forward_commit_basin)
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
        f.write("{} {} {} {}\n".format('NAME', ' '.join(colvar_header), 'RESULT', 'FORWARD_COMMIT_BASIN'))
    return


def read_xyz_frame(ifile, store_atoms=False):
    """Reads a single frame from XYZ file.
    
    Parameters
    ----------
    ifile : TextIOWrapper
        opened file ready for reading
        
    Returns
    -------
    atoms : list
    xyz : np.ndarray
    eof : bool
        bool for if end of file has been reached."""
    
    
    _n_atoms = ifile.readline()
    if _n_atoms:
        _n_atoms = int(_n_atoms)
    else:
        atoms = None
        xyz = None
        eof = True
        return atoms, xyz, eof
    
    if store_atoms:
        atoms = []
    else:
        atoms = None
    
    xyz = np.zeros((_n_atoms, 3))
    # skip comment line
    next(ifile)
    for i in range(_n_atoms):
        line = ifile.readline().split()
        
        if store_atoms:
            atoms.append(line[0])
        else:
            pass
        
        xyz[i] = [float(x) for x in line[1:4]]
    eof = False
    return atoms, xyz, eof


def read_xyz_file(filename):
    """Reads xyz files into an atom list and xyz array
    
    Parameters
    ----------
    filename : str
        path to xyz file
        
    Returns
    -------
    atoms : list
        list containing atomic symbols
    xyz : np.ndarray"""
    
    with open(filename, 'r') as file:
        atoms, _xyz, eof = read_xyz_frame(file, store_atoms=True)
        if eof:
            raise Exception("File at '{}' is empty.".format(filename))
        else:
            pass
        
        n_atoms = len(atoms)
        # initially store each set of xyz coords in a list
        # and then join them all together at the end
        xyz_list = [np.array(_xyz)]
        
        eof = False
        while not eof:
            _, _xyz, eof = read_xyz_frame(file)
            if eof:
                pass
            else:
                xyz_list.append(_xyz)
    xyz = np.array(xyz_list)
    return atoms, xyz

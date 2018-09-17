import subprocess 


def run_MD(inputfile, jobname, topology, logfile, engine="AMBER"):
    """
    Function takes input files and runs MD from them.
    
    Parameters
    ----------
    inputfile : str 
        input file containing simulation details 
    jobname : str
        name for all files associated with run 
    topology : str 
        name of topology file (e.g. for amber it would be  "system.prmtop")
    logfile : str
        file to redirect run output too
    engine: stf
        engine being used for simulation, only takes AMBER inputs
    """ 
    commands = ["sander", "-O",
                "-i", inputfile,
                "-o", jobname + ".out",
                "-x", jobname + ".nc",
                "-p", topology,
                "-c", jobname + ".rst7",
                "-r", jobname + "_out.rst7",
                "-inf", jobname + ".info"]
    subprocess.run(commands, capture_output=True,
                   stdout=jobname + ".log", text=True)

import os
import re
import subprocess

import numpy as np
# import pandas as pd

from .utils import log_run


class AimlessShooting:
    """Runs the aimless shooting algorithm with a pool of starting points.

    Parameters
    ----------
    starting_points : str
        Location of a directory containing structures to start with.
    """
    def __init__(self, starting_points):
        self.starting_points = starting_points
        return
    pass


def run_MD(inputfile, jobname, logfile=None, topology="system.prmtop",
           engine="AMBER"):
    """
    Python wrapper to run MD from commandline

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
                   stdout=logfile, text=True)


def find_and_replace(line, substitutions):
    """
    Function performs find and replace on a string given a dictionary of
    substitution criteria. (Designed to be a sed-like command that can replace
    multiple strings in a single pass)

    Parameters
    ----------
    line : str
        line of text that contains phrases to be replaced
    substitutions: dict
        Dictionary in format {"words_to_find": "words_to_replace"}
    """
    substrings = sorted(substitutions, key=len, reverse=True)
    regex = re.compile('|'.join(map(re.escape, substrings)))
    replace = regex.sub(lambda match: substitutions[match.group(0)], line)
    return replace


class ShootingPoint:
    def __init__(self):
        self.name = None  # Name of shooting point
        self.forward_commit = None  # Status of forward simulation
        self.reverse_commit = None  # Status of reverse simulation
        self.result = None  # Result of shooting point
        self.path = None
        self.cv_values = None
        return

    def run_forward(self, inputfile="fwd.in"):
        """
        Function launches forward simulation based on provided MD input file
        and calls check_if_committed to evaluate if the run ends up in a basin

        Parameters
        ----------
        inputfile : str
            Input file for MD production run
        """
        jobname = self.name + "_f"
        logfile = self.name + "_f.log"
        run_MD(inputfile, jobname, logfile, topology="system.prmtop",
               engine="AMBER")
        self.forward_commit = self.check_if_committed(logfile)
        return

    def run_reverse(self, inputfile="rev.in"):
        """
        Function launches reverse simulation based on provided MD input file
        and calls check_if_committed to evaluate if the run ends up in a basin

        Parameters
        ----------
        inputfile: str
            Input file for MD production run
        """
        jobname = self.name + "_r"
        logfile = self.name + "_r.log"
        run_MD(inputfile, jobname, logfile, topology="system.prmtop",
               engine="AMBER")
        self.reverse_commit = self.check_if_committed(logfile)
        return

    def check_if_committed(self, logfile):
        """
        Function checks a logfile to see if a simulation commited to a basin

        Parameters
        ----------
        logfile : str
            Output of MD simulation (and MD_run) contains PLUMED output

         Returns
        ----------
        status: str/NoneType
            Returns str of basin commited to or none type if simulation did
            not commit
        """
        search_term = "COMMITED"
        status = None
        for line in open(logfile, 'r'):
            if re.search(search_term, line):
                status = int(line.split()[-1].rstrip())
        return status

    def generate_new_shooting_points(self, deltaT, tag,
                                     topology="system.prmtop",
                                     cpptrajskel="cpptraj_skel.in",
                                     cpptrajin="cpptraj.in"):
        """
        Wrapper for cpptraj (AmberTools14). Creates an input file from cpptraj
        skeleton, creates a new restart file shifted deltaT frames away from
        the original shooting point.

        Parameters
        ----------
        cpptrajin : str
            Input file for cpptraj
        deltaT : str
            Frames away from the original shooting point to create new
            restart file
        tag : str
            Identifier taged on to filename to track where the new shooting
            point is created from
        topology : str
            Topology file for system
        cpptrajskel : str
            Skeleton for cpptraj input file
        cpptrajin : str
            Will write to and treat this as the cpptraj input file
        """

        trajectory = self.name + ".nc"
        outfile = self.name + tag + ".rst7"

        topology_text = "PRMTOP"
        trajectory_text = "TRAJECTORY"
        deltaT_text = "FRAME"
        outfile_text = "OUTFILE"

        substitutions = {topology_text: topology, trajectory_text: trajectory,
                         deltaT_text: deltaT, outfile_text: outfile}

        # First create input file (cpptraj.in)
        with open(cpptrajskel, 'r') as file:
            lines = file.readlines()
        with open(cpptrajin, 'w') as file:
            for line in lines:
                file.write(find_and_replace(line, substitutions))

        # Run cpptraj to get new restart file
        commands = ["cpptraj", "-i", cpptrajin]
        subprocess.run(commands)

    def generate_velocity(self, initfile="init.in", logfile=None,
                          topology="system.prmtop", engine="AMBER",
                          solvated=False):
        """
        Function takes shooting point and generates velocities by running MD
        for 1 step with a very small timestep. It then generates the starting
        coordinte and velocity file for the forward and reverse simulation

        Parameters
        ----------
        initfile : str
            Input file, should contain params to run 1 MD step for 0.0001 fs
        logfile : str
            Name of output file
        topology : str
            Name of topology file
        engine : str
            Name of MD engine
        solvated : bool
            If system is solvated the restart file will save box coordinates
        """

        initname = self.name + "_init"
        fwd = self.name + "_f.rst7"
        rev = self.name + "_r.rst7"
        # Run short MD simulation to get self.name_init.rst7 file
        run_MD(initfile, initname, logfile, topology, engine)

        init = initname + ".rst7"
        with open(init, 'r') as init_file:
            init_lines = [line.rstrip('\n') for line in init_file]
        # Write forward restart file (same as init restart files)
        with open(fwd, 'w') as fwd_file:
            job_name = init_lines[0].split()[0]
            num_atoms = init_lines[1].split()[0]
            rst_time = init_lines[1].split()[1]
            fwd_file.write(job_name + os.linesep)
            fwd_file.write("   " + num_atoms + "  " + rst_time + os.linesep)
            for line in range(2, int(num_atoms) + 2):
                fwd_file.write(init_lines[line] + os.linesep)
            # If solvated system print out last line (box coordinates)
            if solvated:
                fwd_file.write(init_lines[-1])
            else:
                print("No box dimensions to print")
                pass
        # Write reverse restart file (with negative velocities)
        with open(rev, 'w') as rev_file:
            rev_file.write(job_name + os.linesep)
            rev_file.write("   " + num_atoms + "  " + rst_time + os.linesep)
            # Copy coordinates over exactly
            coord_lines = int(int(num_atoms) / 2)  # 2 xyz values per row
            for sameline in range(2, coord_lines + 2):
                rev_file.write(init_lines[sameline] + os.linesep)
            for revline in range(coord_lines + 2, coord_lines * 2 + 2):
                fline = map(float, init_lines[revline].split())
                rev_file.write("".join(" % 11.7f" % -num for num in fline))
                rev_file.write(os.linesep)
            # If solvated system print out last line (box coordinates)
            if solvated:
                rev_file.write(init_lines[-1])
            else:
                print("No box dimensions to print")
                pass

    def _read_cv_values(self, colvar_file='COLVAR'):
        """Grabs the values of all CVs in first line of COLVAR.

        Parameters
        ----------
        colvar_file : str, optional
            Name of colvar file.

        Returns
        -------
        data : numpy.ndarray
            Values of CVs in first line, excluding time value."""
        data = np.genfromtxt(colvar_file, max_rows=1)
        self.cv_values = data[1:]
        return

    def log(self, filename='log'):
        """Logging method."""
        if self.cv_values:
            pass
        else:
            self._read_cv_values()
        log_run(self.name, self.cv_values, self.result, filename=filename)
        return

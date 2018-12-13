import os
import os.path as op
import re
import shutil
import subprocess

import numpy as np

from .utils import log_run
from .utils import log_header


class AimlessShooting:
    """Runs the aimless shooting algorithm with a pool of starting points.

    Parameters
    ----------
    starting_points : str
        Location of a directory containing structures to start with.
    deltaT : int, optional
        Value of deltaT for creating new shooting points.
    logfile : str, optional
        Name of logfile.
    """
    def __init__(self, starting_points, accepts_goal=100, deltaT=10,
                 logfile='log', colvar_name='COLVAR'):
        self.starting_points = starting_points
        self.logfile = logfile
        self.queue = []
        self.guesses = []
        self.num_accepts = 0
        self.accepts_goal = accepts_goal
        self.deltaT = deltaT
        self.colvar_name = colvar_name
        # self.counter = 0
        return

    def generate_guesses(self):
        filenames = os.listdir('starting_structures')
        for f in filenames:
            self.guesses.append(op.splitext(f)[0])
        return

    def start(self):
        # TODO: if restarting, first check if guesses is already populated.
        self.generate_guesses()

        while self.num_accepts < self.accepts_goal:
            # Initialize a shooting point
            sp = self.initialize_shooting_point()

            # Generate velocities
            sp.generate_velocity()

            # Run forward simulation
            sp.run_forward()

            # Create log file and write header.
            if op.exists(colvar_file):
                log_header(colvar_file, filename=logfile)
            else:
                raise Exception("{} file does not exist ".format(colvar_file)
                                + "and header cannot be extracted.")

            # Check if forward simulation commits to basin
            if sp.forward_commit is None:
                # Log results
                sp.results = 'inconclusive'
                sp.log(logfile)
                continue

                # Check to see how many times sp has attempted to commit
                # if attempts < attempt_criteria:
                    # Restart the shooting point with new velocities
                    # restart = self.initialize_shooting_point()  # **this line might break
                    # self.queue.append(restart)
                # else:
                    # If there have been too many attempts on the same sp
                    # pass
            else:
                # Continue to run the reverse simulation
                sp.run_reverse()
                # If reverse simulation does not commit to a basin
                if sp.reverse_commit is None:
                    sp.result = 'inconclusive'
                    sp.log(logfile)  # Log run
                    continue
                    # initialize a new --- RES
                # If reverse simulation commits to the same basin as forward
                elif sp.reverse_commit == sp.forward_commit:
                    sp.result = 'reject'
                    sp.log(logfile)
                    # launch new trajectory from queue
                elif sp.reverse_commit != sp.forward_commit:
                    sp.result = 'accept'
                    sp.log(logfile)
                    self.num_accepts += 1
                    # Generate 3 new shooting points
                    job1 = sp.generate_new_shooting_points(deltaT, tag)
                    self.queue.append(job1)
                    job2 = sp.generate_new_shooting_points(deltaT, tag)
                    self.queue.append(job2)
                    job3 = sp.generate_new_shooting_points(deltaT, tag)
                    self.queue.append(job3)
                else:
                    raise Exception("Unknown shooting point result.")

            # Increment the job counter.
            # self.counter += 1
        return

    def initialize_shooting_point(self):
        """Creates new shooting point.

        This will always check the `self.queue` variable for new
        ShootingPoints. If `self.queue` is empty, it will instead move
        on to the next structure in `self.guesses`.
        
        Returns
        -------
        sp : ShootingPoint
            New ShootingPoint object.
        """
        # Check queue first.
        if self.queue:
            directory = './queue'
            print('current directory:', directory)
            name = self.queue[0]
            filename = name + '.rst7'
            shutil.copyfile(op.join(directory, filename), filename)
            del self.queue[0]
        else:
            directory = './guesses'
            print('current directory:', directory)
            name = self.guesses[0]
            filename = name + '.rst7'
            shutil.copyfile(op.join(directory, filename), filename)
            del self.guesses[0]
        sp = ShootingPoint(name, input_init=name+'.rst7',
                           topology_file="system.prmtop",
                           md_engine="AMBER")
        return sp


def run_MD(inputfile, jobname, output=None, topology="system.prmtop",
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
    output : str
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
                   stdout=output, text=True)


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
    def __init__(self, name, input_md="nvt.in", input_init="init.in", topology_file=None,
                 md_engine="AMBER"):
        self.name = name  # Name of shooting point
        self.input_md = input_md  # Name of AMBER input file for running MD.
        self.input_init = input_init  # Name of AMBER input file for velocity generation.
        self.topology_file = topology_file  # Path to topology file
        self.md_engine = md_engine
        # self._input_file_dir = op.dirname(op.abspath(self.input_file))
        # self._input_file_name = op.splitext(op.basename(self.input_file))[0]

        self.forward_commit = None  # Status of forward simulation
        self.reverse_commit = None  # Status of reverse simulation
        self.result = None  # Result of shooting point
        self.cv_values = None
        return

    def run_forward(self):
        """
        Function launches forward simulation based on provided MD input file
        and calls check_if_committed to evaluate if the run ends up in a basin

        Parameters
        ----------
        inputfile : str
            Input file for MD production run
        """
        jobname = self.name + "_f"
        output = self.name + "_f.log"
        run_MD(self.input_md, jobname, output, topology="system.prmtop",
               engine="AMBER")
        self.forward_commit = self.check_if_committed(output)

        # Copy plumed log file to new directory with unique name.
        if op.exists('plumed_log'):
            os.copy(output, op.join('plumed_log', output))
        else:
            raise Exception("There is no 'plumed_log' directory.")
        return

    def run_reverse(self):
        """
        Function launches reverse simulation based on provided MD input file
        and calls check_if_committed to evaluate if the run ends up in a basin
        """
        jobname = self.name + "_r"
        output = self.name + "_r.log"
        run_MD(self.input_md, jobname, output, topology="system.prmtop",
               engine="AMBER")
        self.reverse_commit = self.check_if_committed(output)

        # Copy plumed log file to new directory with unique name.
        if op.exists('plumed_log'):
            os.copy(output, op.join('plumed_log', output))
        else:
            raise Exception("There is no 'plumed_log' directory.")
        return

    def check_if_committed(self, output):
        """
        Function checks an output file to see if a simulation commited to a basin

        Parameters
        ----------
        output : str
            Output of MD simulation (and MD_run) contains PLUMED output

         Returns
        ----------
        status: str/NoneType
            Returns str of basin commited to or none type if simulation did
            not commit
        """
        search_term = "COMMITED"
        status = None
        for line in open(output, 'r'):
            if re.search(search_term, line):
                status = int(line.split()[-1].rstrip())
        return status

    def generate_new_shooting_points(self, deltaT, direction,
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
        direction : str
            Identifier taged on to filename to track where the new shooting
            point is created from
        topology : str
            Topology file for system
        cpptrajskel : str
            Skeleton for cpptraj input file
        cpptrajin : str
            Will write to and treat this as the cpptraj input file

         Returns
        ----------
        outfile: str
            Returns a strig of the restart file created from cpptraj
        """
        if direction == 'fwd':
            name = self.name + "_f"
            tag = "_f"
        elif direction == 'rev':
            name = self.name + "_r"
            tag = "_r"
        trajectory = name + ".nc"
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

        return outfile

    def generate_velocity(self, output=None, topology="system.prmtop",
                          engine="AMBER", solvated=False):
        """
        Function takes shooting point and generates velocities by running MD
        for 1 step with a very small timestep. It then generates the starting
        coordinte and velocity file for the forward and reverse simulation

        Parameters
        ----------
        output : str
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
        run_MD(self.input_init, initname, output, topology, engine)

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

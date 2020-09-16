import os
import os.path as op
import re
import shutil
import subprocess
import sys

# import mdtraj as md
import numpy as np

from .utils import log_run, log_header, read_xyz_file, write_xyz_frame

# Make sure to copy the original shooting point for the "_0" instance!!


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
                 retries=10, max_count=2000, logfile='log', cores=1,
                 colvar_name='COLVAR', md_engine='CP2K', exe=None,
                 mpiexec=None):
        self.starting_points = starting_points
        self.logfile = logfile
        self.queue = []
        self.guesses = []
        self.num_accepts = 0
        self.retries = retries
        self.max_count = max_count
        self.accepts_goal = accepts_goal
        self.deltaT = deltaT
        self.colvar_name = colvar_name
        self.cores = cores
        self.md_engine = md_engine
        self.exe = exe
        self.mpiexec = mpiexec
        # self.counter = 0
        return

    def generate_guesses(self):
        filenames = os.listdir(self.starting_points)
        for f in filenames:
            name = op.splitext(f)[0]
            self.guesses.append(name)
            if op.exists(op.join('guesses', f)):
                pass
            else:
                shutil.copyfile(op.join(self.starting_points, f),
                                op.join('guesses', f))
        return

    def start(self):
        # TODO: write a function to clean old files and COLVARs

        # Generate all necessary folders if they do not already exist.
        required_directories = ['guesses', 'queue', 'COLVARS', 'plumed_log']
        for directory in required_directories:
            if op.exists(directory):
                pass
            else:
                os.mkdir(directory)

        trial_num = 0
        counter = 0
        # max_count = 10000
        max_count_reached = False

        # while self.num_accepts < self.accepts_goal:
        # while trial_num < self.retries:
        for trial_num in range(self.retries):
            # TODO: if restarting, first check if guesses is already populated.
            self.generate_guesses()
            while self.queue or self.guesses:  # counter < max_count:
                # Initialize a shooting point
                sp = self.initialize_shooting_point(self.md_engine)

                # If `self.initialize_shooting_point()` cannot generate a
                # new shooting point, then it returns `1`.
                if sp == 1:
                    break
                else:
                    pass

                # Generate velocities
                # sp.generate_velocity(None, sp.atoms)

                # Run forward simulation
                sp.run_forward()

                # Create log file and write header.
                if not op.exists(self.logfile):
                    if op.exists(self.colvar_name):
                        log_header(self.colvar_name, filename=self.logfile)
                    else:
                        raise Exception("{} file does not exist "
                                        .format(self.colvar_name)
                                        + "and header cannot be extracted.")
                else:
                    pass

                # Check if forward simulation commits to basin
                if sp.forward_commit is None:
                    # Log results
                    sp.result = 'inconclusive'
                    sp.log(self.logfile)
                    save_colvar(self.colvar_name, sp.name + '_f')

                    # Check to see how many times sp has attempted to commit
                    # if attempts < attempt_criteria:
                        # Restart the shooting point with new velocities
                        # restart = self.initialize_shooting_point()  # **this line might break
                        # self.queue.append(restart)
                    # else:
                        # If there have been too many attempts on the same sp
                        # pass
                else:
                    # Save the COLVAR of the previous run before starting new run.
                    save_colvar(self.colvar_name, sp.name + '_f')
                    # Continue to run the reverse simulation
                    sp.run_reverse()
                    # If reverse simulation does not commit to a basin
                    if sp.reverse_commit is None:
                        sp.result = 'inconclusive'
                        sp.log(self.logfile)  # Log run
                        # initialize a new --- RES
                    # If reverse simulation commits to the same basin as forward
                    elif sp.reverse_commit == sp.forward_commit:
                        sp.result = 'reject'
                        sp.log(self.logfile)
                        # launch new trajectory from queue
                    elif sp.reverse_commit != sp.forward_commit:
                        sp.result = 'accept'
                        sp.log(self.logfile)
                        self.num_accepts += 1

                        # Generate 3 new shooting points
                        # TODO: create variable that contains extension based
                        #       on which MD engine so checking is not necessary
                        job1 = sp.generate_new_shooting_points(self.deltaT, 'fwd')
                        self.queue.append(job1)
                        if sp.md_engine == "AMBER":
                            job1name = job1 + '.rst7'
                        elif sp.md_engine == "CP2K":
                            job1name = job1 + '.xyz'
                        else:
                            raise Exception("MD engine not recognized.")
                        shutil.move(job1name, op.join('queue', job1name))

                        job2 = sp.generate_new_shooting_points(self.deltaT, 'rev')
                        self.queue.append(job2)
                        if sp.md_engine == "AMBER":
                            job2name = job2 + '.rst7'
                        elif sp.md_engine == "CP2K":
                            job2name = job2 + '.xyz'
                        else:
                            raise Exception("MD engine not recognized.")
                        shutil.move(job2name, op.join('queue', job2name))

                        job3 = sp.generate_new_shooting_points(self.deltaT, '0')
                        self.queue.append(job3)
                        if sp.md_engine == "AMBER":
                            job3name = job3 + '.rst7'
                        elif sp.md_engine == "CP2K":
                            job3name = job3 + '.xyz'
                        else:
                            raise Exception("MD engine not recognized.")
                        shutil.move(job3name, op.join('queue', job3name))
                    else:
                        raise Exception("Unknown shooting point result.")
                    save_colvar(self.colvar_name, sp.name + '_r')

                # Delete all run files after the shooting point is complete.
                self.remove_files(sp.name)
                # Increment the job counter.
                counter += 1
                if counter == self.max_count:
                    max_count_reached = True
                    with open(self.logfile, 'a') as f:
                        f.write("# Max count reached ({}). ".format(self.max_count) +
                                "No new shooting points will be generated.\n")
                else:
                    pass
                # End of while loop

            # If max_count is already reached, then do not restart the run.
            if max_count_reached:
                break
            else:
                pass
        return
                # self.counter += 1
            # if counter == max_count:
                # break
            # else:
                # pass
        # return

    def remove_files(self, basename):
        name_tags = ['_r', '_f', '_init', '_init']
        out_tags = ['', '_out']
        extensions = ['.rst7', '.out', '.nc', '.log', '.info', '.xyz']
        for tag in name_tags:
            for out_tag in out_tags:
                for ext in extensions:
                    if op.exists(basename + tag + out_tag + ext):
                        os.remove(basename + tag + out_tag + ext)
                    else:
                        pass
        return

    def initialize_shooting_point(self, engine):
        """Creates new shooting point.

        This will always check the `self.queue` variable for new
        ShootingPoints. If `self.queue` is empty, it will instead move
        on to the next structure in `self.guesses`.
        
        Returns
        -------
        sp : ShootingPoint
            New ShootingPoint object.
        """
        if engine == "AMBER":
            sp = self.initialize_AMBER_shooting_point()
        elif engine == "CP2K":
            sp = self.initialize_CP2K_shooting_point()
        else:
            raise Exception("Engine not recognized: {}".format(engine))
        return sp

    def initialize_AMBER_shooting_point(self):
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
            # print('current directory:', directory)
            name = self.queue[0]
            filename = name + '.rst7'
            target = name + '_init.rst7'
            shutil.copyfile(op.join(directory, filename), target)
            del self.queue[0]
        elif self.guesses:
            directory = './guesses'
            # print('current directory:', directory)
            name = self.guesses[0]
            filename = name + '.rst7'
            target = name + '_init.rst7'
            shutil.copyfile(op.join(directory, filename), target)
            del self.guesses[0]
        else:
            with open(self.logfile, 'a') as f:
                f.write("# No more guesses available. Exiting.\n")
            return 1
            # sys.exit()

        sp = ShootingPoint(name,
                           topology_file="system.prmtop",
                           md_engine="AMBER")
        return sp

    def initialize_CP2K_shooting_point(self):
        """Creates new shooting point in CP2K.

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
            name = self.queue[0]
            del self.queue[0]
        elif self.guesses:
            directory = './guesses'
            name = self.guesses[0]
            del self.guesses[0]
        else:
            with open(self.logfile, 'a') as f:
                f.write("# No more guesses available. Exiting.\n")
            return 1

        coordfilename = name + '.xyz'
        shutil.copyfile(op.join(directory, coordfilename), coordfilename)

        inputfilename = name + '.inp'
        if directory == './guesses':
            atoms, xyz = read_xyz_file(op.join(directory, name + '.xyz'))
            vel = get_maxwell_boltzmann_velocities(atoms, T=300,
                                                   dim=3, units='au')
            generate_cp2k_input_file(inputfilename,
                                     jobname=name,
                                     coord_file=coordfilename,
                                     wfn_restart_file=name + ".wfn",
                                     velocities=vel,
                                     template="template.inp")
        else:
            shutil.copyfile(op.join(directory, inputfilename), inputfilename)

        sp = ShootingPoint(name, input_md=inputfilename,
                           topology_file=coordfilename,
                           md_engine="CP2K", exe=self.exe, mpiexec=self.mpiexec)
        return sp


def run_MD(inputfile, jobname, output=None, cores=1, topology="system.prmtop",
           engine="AMBER", exe=None, mpiexec=None):
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
    engine: str
        engine being used for simulation; options are AMBER or CP2K;
        default is AMBER
    """
    if engine == "AMBER":
        run_AMBER(inputfile, jobname, output, cores, topology)
    elif engine == "CP2K":
        run_CP2K(inputfile, jobname, output, cores, exe, mpiexec)
    else:
        raise Exception("Engine '{}' is not supported, please".format(engine) +
                        " specify either 'AMBER' or 'CP2K'.")
    return


def run_AMBER(inputfile, jobname, output=None, cores=1,
              topology="system.prmtop"):
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
    engine: str
        engine being used for simulation, only takes AMBER inputs
    """
    commands = ["mpirun", "-np", str(cores),
                "sander.MPI", "-O",
                "-i", inputfile,
                "-o", jobname + ".out",
                "-x", jobname + ".nc",
                "-p", topology,
                "-c", jobname + ".rst7",
                "-r", jobname + "_out.rst7",
                "-inf", jobname + ".info"]
    if output:
        logf = open(output, 'w')
    else:
        logf = subprocess.DEVNULL

    subprocess.run(commands, stdout=logf)

    if output:
        logf.close()
    else:
        pass
    return


def run_CP2K(inputfile, jobname, output=None, cores=1, exe=None, mpiexec=None):
    """
    Python wrapper to run CP2K MD from commandline

    Parameters
    ----------
    inputfile : str
        input file containing simulation details
    jobname : str
        name for all files associated with run
    output : str
        file to redirect run output too
    exe : str
        path to the executable for CP2K
    mpiexec : str
        executable for MPI
    """
    assert op.exists(exe), "CP2K executable not found at: {}".format(exe)

    commands = [
                # mpiexec, "-np", str(cores),
                exe, 
                "-i", inputfile,
                "-o", jobname + ".out",
                ]
    if mpiexec:
        commands.insert(0, str(cores))
        commands.insert(0, "-np")
        commands.insert(0, mpiexec)
    else:
        pass

    if output:
        logf = open(output, 'w')
    else:
        logf = subprocess.DEVNULL

    subprocess.run(commands, stdout=logf)

    if output:
        logf.close()
    else:
        pass
    return


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
    def __init__(self, name, input_md="nvt.in", input_init="init.in", cores=1,
                 topology_file=None, md_engine="AMBER", exe=None, mpiexec=None):
        self.name = name  # Name of shooting point
        self.input_md = input_md  # Name of AMBER input file for running MD.
        self.input_init = input_init  # Name of AMBER input file for velocity generation.
        self.topology_file = topology_file  # Path to topology file
        self.md_engine = md_engine
        if md_engine == 'CP2K':
            self.atoms = self.get_atoms_list(topology_file)
        else:
            self.atoms = []
        self.cores = cores
        self.exe = exe
        self.mpiexec = mpiexec
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
        run_MD(self.input_md, jobname, output, cores=self.cores,
               topology=self.topology_file, engine=self.md_engine, exe=self.exe,
               mpiexec=self.mpiexec)
        self.forward_commit = self.check_if_committed(output)

        # Copy plumed log file to new directory with unique name.
        if op.exists('plumed_log'):
            shutil.move(output, op.join('plumed_log', output))
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
        run_MD(self.input_md, jobname, output, cores=self.cores,
               topology=self.topology, engine=self.md_engine, exe=self.exe,
               mpiexec=self.mpiexec)
        self.reverse_commit = self.check_if_committed(output)

        # Copy plumed log file to new directory with unique name.
        if op.exists('plumed_log'):
            shutil.move(output, op.join('plumed_log', output))
        else:
            raise Exception("There is no 'plumed_log' directory.")

        # Copy COLVAR file to new directory with unique name.
        # if op.exists('COLVARS'):
            # shutil.move('COLVAR', op.join('COLVARS', jobname + '_COLVAR'))
        # else:
            # raise Exception("There is no 'COLVARS' directory.")
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

    def generate_new_shooting_points(self, deltaT, direction):
        if self.md_engine == "AMBER":
            return self.generate_new_AMBER_shooting_points(self, deltaT,
                                                           direction,
                                                           self.topology)
        elif self.md_engine == "CP2K":
            return self.generate_new_CP2K_shooting_points(self, deltaT,
                                                          direction)
        else:
            raise Exception("Engine '{}' is not supported, ".format(engine) +
                            "please specify either 'AMBER' or 'CP2K'.")
        return

    def generate_new_AMBER_shooting_points(self, deltaT, direction,
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
            tag = "_f"
            name = self.name + tag
        elif direction == 'rev':
            tag = "_r"
            name = self.name + tag
        elif direction == '0':
            # Copies the original rst7 file and renames it.
            tag = "_0"
            name = self.name + tag
            shutil.copyfile(self.name + "_init.rst7", name + ".rst7")
            return name
        trajectory = name + ".nc"
        outfile = name + ".rst7"
        outfile_jobname = name

        topology_text = "PRMTOP"
        trajectory_text = "TRAJECTORY"
        deltaT_text = "FRAME"
        outfile_text = "OUTFILE"

        if isinstance(deltaT, str):
            pass
        else:
            deltaT = str(deltaT)

        substitutions = {topology_text: topology, trajectory_text: trajectory,
                         deltaT_text: deltaT, outfile_text: outfile}

        # First create input file (cpptraj.in)
        with open(cpptrajskel, 'r') as file:
            lines = file.readlines()
        with open(cpptrajin, 'w') as file:
            for line in lines:
                file.write(find_and_replace(line, substitutions))

        # Run cpptraj to get new restart file
        commands = ["cpptraj.MPI", "-i", cpptrajin]
        subprocess.run(commands, stdout=subprocess.DEVNULL)

        return outfile_jobname

    def generate_new_CP2K_shooting_points(self, deltaT, direction):
        """
        Generates new shooting points from original xyz trajectory

        Reads in CP2K xyz trajectory file and writes out frames that are
        deltaT away from original shooting point, or simply copies the original
        shooting point.

        Parameters
        ----------
        deltaT : str
            Frames away from the original shooting point to create new
            restart file
        direction : str
            Identifier taged on to filename to track where the new shooting
            point is created from

         Returns
        ----------
        outfile: str
            Returns a string of the new jobname
        """
        if direction == 'fwd':
            tag = "_f"
        elif direction == 'rev':
            tag = "_r"
        elif direction == '0':
            tag = "_0"
            shutil.copyfile(self.name + ".xyz", self.name + tag + ".xyz")
        name = self.name + tag
        trajectory = name + ".xyz"
        outfile_jobname = name

        atoms, xyz = read_xyz_file(self.name + '.xyz')
        if isinstance(deltaT, int):
            pass
        else:
            deltaT = int(deltaT)

        write_xyz_frame(trajectory, atoms, xyz[deltaT],
                        comment='dt = {}'.format(deltaT))
        vel = get_maxwell_boltzmann_velocities(atoms, T=300,
                                               dim=3, units='au')
        generate_cp2k_input_file(name + '.inp',
                                 jobname=name,
                                 coord_file=trajectory,
                                 wfn_restart_file=name + ".wfn",
                                 velocities=vel,
                                 template="template.inp")
        return outfile_jobname

    def generate_velocity(self, output, atoms, T=300):
        if self.md_engine == "AMBER":
            self.generate_AMBER_velocity(output=output)
            return
        else:
            pass

        vel = get_maxwell_boltzmann_velocities(atoms, T=T,
                                               dim=3, units='au')
        write_xyz_frame(self.name + "_vel_f.xyz", atoms, vel)
        write_xyz_frame(self.name + "_vel_r.xyz", atoms, -vel)
        return

    def generate_AMBER_velocity(self, output=None, topology="system.prmtop",
                                engine="AMBER", solvated=True):
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
        run_MD(self.input_init, initname, output=output, cores=self.cores,
               topology=topology, engine=engine)
        
        # Remove the COLVAR file because it is unnecessary.
        if op.exists('COLVAR'):
            os.remove('COLVAR')
        else:
            pass

        init = initname + "_out.rst7"
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
        log_run(self.name, self.cv_values, self.result,
                self.forward_commit, filename=filename)
        return

    def get_atoms_list(self, topology):
        atoms, xyz = read_xyz_file(topology)
        return atoms


def save_colvar(colvar_name, jobname):
    # Copy COLVAR file to new directory with unique name.
    if op.exists('COLVARS'):
        shutil.move(colvar_name, op.join('COLVARS', jobname + '_COLVAR'))
    else:
        raise Exception("There is no 'COLVARS' directory.")
    return


def get_maxwell_boltzmann_velocities(atoms, T=300, dim=3, units='au'):
    """Generates atomic velocities from a Maxwell-Boltzmann distribution.

    Parameters
    ----------
    atoms : array-like of str
        Array containing atomic symbols.
    T : float, optional
        Reference temperature for Maxwell-Boltzmann distribution.
        Default is 300 K.
    dim : int, optional
        Specifies the dimensionality of the velocities. Default is
        3 dimensions.
    units : str, optional
        Specifies units for velocity. Default is au (bohr / au_time).
        Possible choices are ['au', 'm/s'].

    Returns
    -------
    v : np.ndarray
        Atomic velocities in specified number of dimensions.
        Units of `v` are in au (bohr / au_time)."""

    kB = 1.380649e-23            # J / K
    au_time_factor = 0.0242e-15  # s / au_time
    bohr_factor = 5.29e-11       # m / bohr

    mass = atomic_symbols_to_mass(atoms)
    n_atoms = len(mass)

    # Number of degrees of freedom
    if dim == 3:
        dof = dim * n_atoms - 6
    elif dim == 2:
        dof = dim * n_atoms - 3
    elif dim == 1:
        dof = dim * n_atoms - 1
    else:
        raise Exception("Number of dimensions must be 3, 2, or 1")

    # Convert mass from amu to kg
    mass = np.asarray(mass).reshape(-1, 1) / 1000 / 6.022e23

    v_raw = np.sqrt(kB * T / mass) * np.random.normal(size=(n_atoms, dim))

    # Shift velocities by mean momentum such that total
    # box momentum is 0 in all dimensions.
    p_mean = np.sum(mass * v_raw, axis=0) / n_atoms
    v_mean = p_mean / mass.mean()
    v_shifted = v_raw - v_mean
    
    # Scale velocities to exact target temperature.
    # Prevents systems with few atoms from sampling far away from target T.
    T_shifted = np.sum(mass * v_shifted ** 2) / (dof * kB)
    scale = np.sqrt(T / T_shifted)
    v_scaled = v_shifted * scale

    # Unit conversion from m/s to final units
    if units == 'm/s':
        v = v_scaled
    elif units == 'au':
        v = v_scaled * au_time_factor / bohr_factor
    else:
        raise Exception("Invalid units for velocity: {}.".format(units)
                        + " Choose from ['au', 'm/s'].")
    return v


def atomic_symbols_to_mass(atoms):
    """Converts atomic symbols to their atomic masses in amu."""
    atomic_mass_dict = get_atomic_mass_dict()
    masses = []
    for atom in atoms:
        masses.append(atomic_mass_dict[atom])
    return masses


def get_atomic_mass_dict():
    atomic_mass_dict = {}
    dir_path = op.dirname(op.realpath(__file__))
    with open(op.join(dir_path, 'data', 'atomic_info.dat')) as f:
        for line in f:
            data = line.split()
            if data[-1].startswith('('):
                data[-1] = data[-1][1:-1]
            atomic_mass_dict[data[1]] = float(data[-1])
    return atomic_mass_dict


def generate_cp2k_input_file(input_file, jobname, coord_file,
                             wfn_restart_file, velocities, template):

    vel_subst = '\n'.join(["       {: .8e} {: .8e} {: .8e}".format(vx, vy, vz)
                           for vx, vy, vz in velocities])

    subst_dict = {"jobname": jobname,
                  "system_coords": coord_file,
                  "wfn_restart": wfn_restart_file,
                  "INSERT_VELOCITIES_HERE": vel_subst}

    with open(input_file, 'w') as new_file:
        with open(template) as old_file:
            for line in old_file:
                pattern = None
                subst = None

                for _pattern, _subst in subst_dict.items():
                    if _pattern in line:
                        pattern = _pattern
                        subst = _subst
                    else:
                        pass

                if pattern:
                    new_file.write(line.replace(pattern, subst))
                else:
                    new_file.write(line)
    return

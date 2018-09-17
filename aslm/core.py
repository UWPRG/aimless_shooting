import subprocess 


def run_MD(inputfile, jobname, logfile=None, topology="system.prmtop", engine="AMBER"):
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
                   stdout=logfile, text=True)


class ShootingPoint:
    def __init__(self):
        self.name = None
        self.forward_commit = None
        self.backward_commit = None
        self.run_status = None
        self.path = None
        return

    def run_forward(self):
        run_MD(params)
        return

    def run_backward(self):
        run_MD(params)
        return

    def check_if_committed(self):
         
        return

    def generate_new_shooting_points(self):
        return

    def generate_velocity(self, initfile="init.in", logfile=None,
                          topology="system.prmtop", engine="AMBER"):
        """
        Function takes shooting point and generates velocities by running MD 
        for 1 step with a very small timestep. It then generates the starting 
        coordinte and velocity file for the forward and reverse simulation 
        """

        initname = self.name + "_init"
        fwd = self.name + "_f.rst7"
        rev = self.name + "_r.rst7"
        # Run short MD simulation to get self.name_init.rst7 file
        run_MD(initfile, initname, logfile, topology, engine)) 

        init = initname + ".rst7"
        with open(init, 'r') as init_file:
            init_lines = [line.rstrip('\n') for line in init_file]
        #write forward restart file (same as init restart files)
        with open(fwd, 'w') as fwd_file:
            job_name = init_lines[0].split()[0]
            num_atoms = init_lines[1].split()[0]
            rst_time = init_lines[1].split()[1]
            fwd_file.write(job_name + os.linesep)
            fwd_file.write("   " + num_atoms + "  " + rst_time + os.linesep)
            for line in range(2, int(num_atoms)+2):
                fwd_file.write(init_lines[line] + os.linesep)
            # if solvated run print out last line (box coordinates)
            if solvated == True:
                fwd_file.write(init_lines[-1])
            else:
                print ("No box dimensions to print")
                pass
        #write reverse restart file (with negative velocities)
        with open(rev, 'w') as rev_file:
            rev_file.write(job_name + os.linesep)
            rev_file.write("   " + num_atoms + "  " + rst_time + os.linesep)
            #copy coordinates over exactly
            coord_lines = int(int(num_atoms)/2) #2 xyz values per row
            for sameline in range(2, coord_lines+2):
                rev_file.write(init_lines[sameline] + os.linesep)
            for revline in range(coord_lines+2,coord_lines*2+2):
                fline = map(float, init_lines[revline].split())
                rev_file.write("".join(" % 11.7f" % -num for num in fline))
                rev_file.write(os.linesep)
            # if solvated run print out last line (box coordinates)
            if solvated == True:
                rev_file.write(init_lines[-1])
            else:
                print ("No box dimensions to print")
                pass

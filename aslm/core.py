import os 

def run_MD(name,inputfile,prmtop,logfile):
    """
    Function takes a name and runs MD
    """
    sander -O -i inputfile.in -o name.out -x name.nc -p system.prmtop -c name.rst7 -r name_out.rst7 -inf name.info >> log.log



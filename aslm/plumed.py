from itertools import combinations

def gen_plumed_file(atoms,basin_CV,print_stride='1',plumed_file='plumed.dat',
                       plumed_outfile='COLVAR'):
    """
    Function will take a list of atoms, COMS, and groups and create a plumed
    file containing distances, differences of distances, angles, dihedrals,
    and other interesting CVs

    Parameters
    ----------
    atoms : list of integers
        Given a list of important atoms, function will generate a plumed entry
        for each set of distances
        Example: [1,2,3,4]

    basin_CVs : list of strings
        List of variables that correspond to the CVs to be monitored for each 
        basin. Only accepts 1-2 CVs as an input.
        Example: ['phi','psi']

    """
    # Create pairs of 2 to describe distance
    pairs = list(combinations(atoms,2))
    # Create lists of 3 to describe angles
    triples = list(combinations(atoms,3))
    # Create lists of 4 to describe dihedrals
    quads = list(combinations(atoms,4))

    f=open(plumed_file,"w+") # Open file for writing 
    f.write('# vim:ft=plumed\n')
    if len(atoms) >= 4:
        for k in quads:
            f.write('\n# TORSION \n')
            f.write('t{}{}{}{}: TORSION ATOMS={},{},{},{}\n'.format(k[0],k[1],
                    k[2],k[3],k[0],k[1],k[2],k[3]))
    else:
        print("Must include at least 4 atoms to calculate torsion")
    if len(atoms) >= 3:
        f.write('\n# ANGLES \n')
        for j in triples: 
            f.write('a{}{}{}: ANGLE ATOMS={},{},{}\n'.format(j[0],j[1],j[2],
                    j[0],j[1],j[2]))
    else:
        print("Must include at least 3 atoms to calculate angles")
    if len(atoms) >= 2:
        f.write('\n# DISTANCES \n')
        for i in pairs:
            f.write('d{}{}: DISTANCE ATOMS={},{}\n'.format(i[0],i[1],i[0],
                    i[1]))
    else:
        print("Must include at least two atoms in the atom namelist to generate distances")

    f.write('\nCOMMITTOR ...\n')
    if len(basin_CV) == 1:
        f.write('ARG={}\n'.format(basin_CV[0]))
    if len(basin_CV) == 2:
        f.write('ARG={},{}\n'.format(basin_CV[0],basin_CV[1]))
    if len(basin_CV) >= 3:
        print ("Use of more than 2 CVs is not supported in this script", 
                "you may manually edit the plumed file")
    if len(basin_CV) < 1:
        print ("You must provide a CV to monitor if the shooting point has committed to a basin")
    f.write('... COMMITTOR\n')

    f.write('\nPRINT STRIDE={} ARG=* FILE={}'.format(print_stride,plumed_outfile))
    f.close() # Close file for writing 


# More CVs to add 

# difference of distances between CVs 
# cremer pople
# Hill-Reilly?
# Berces et. al36 parameter 1

# Finish scripting COMMITTOR    
# define basin 

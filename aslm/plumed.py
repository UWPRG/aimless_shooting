from itertools import combinations

def gen_plumed_file(atoms,basins,basin_CV,limit,print_stride='1',plumed_file='plumed.dat',
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
    
    basins : int
        Number of basins used to monitor COMMITTOR

    basin_CVs : list of strings
        List of variables that correspond to the CVs to be monitored for each 
        basin. Only accepts 1-2 CVs as an input.
        Example: ['phi','psi']

    limit : array 
        Array of integers that correspond to the lower and upper limits for 
        each basin. Must be enetered in a specific format
        Example for 2 basins: [[LL1_x,LL1_y,UL1_x,UL1_y],[LL2_x,LL2_y,UL2_x,UL2_y]]
        * where x and y correspond to those values for each supplied argument
        * 1 and 2 correspond to each basin
        * LL = lower limit 
        * UL = upper limit

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
    if basins == 2:
        f.write('ARG={}\n'.format(basin_CV[0], basin_CV[1]))
        f.write('STRIDE={}\n'.format(print_stride))
        f.write('BASIN_LL1={},{}\n'.format(limit[0][0], limit[0][1]))
        f.write('BASIN_UL1={},{}\n'.format(limit[0][2], limit[0][3]))
        f.write('BASIN_LL2={},{}\n'.format(limit[1][0], limit[1][1]))
        f.write('BASIN_UL2={},{}\n'.format(limit[1][2], limit[1][3]))
    if basins == 3:
        f.write('ARG={},{}\n'.format(basin_CV[0],basin_CV[1],basin_CV[3]))
        f.write('STRIDE={}\n'.format(print_stride))
        f.write('BASIN_LL1={},{}\n'.format(limit[0][0], limit[0][1]))
        f.write('BASIN_UL1={},{}\n'.format(limit[0][2], limit[0][3]))
        f.write('BASIN_LL2={},{}\n'.format(limit[1][0], limit[1][1]))
        f.write('BASIN_UL2={},{}\n'.format(limit[1][2], limit[1][3]))
        f.write('BASIN_LL3={},{}\n'.format(limit[2][0], limit[2][1]))
        f.write('BASIN_UL3={},{}\n'.format(limit[2][2], limit[2][3]))
    if len(basin_CV) >= 4:
        print ("Use of more than 3 CVs to monitor COMMITTOR is not currently supported", 
                "you may manually edit the plumed file")
    if len(basin_CV) < 2:
        print ("You must provide at least two CVs to monitor if the shooting point has committed to a basin")
    f.write('... COMMITTOR\n')

    f.write('\nPRINT STRIDE={} ARG=* FILE={}'.format(print_stride,plumed_outfile))
    f.close() # Close file for writing 


# More CVs to add 

# difference of distances between CVs 
# cremer pople
# Hill-Reilly?
# Berces et. al36 parameter 1

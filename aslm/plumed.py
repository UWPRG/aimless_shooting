from itertools import combinations

# Change to class 

class GeneratePlumed:
    """
    Class instantiates a plumed file and has functionality to add candidate 
    CVs such as distances, differences of distances, angles, dihedrals, etc.

    Parameters
    ----------
    print_stride : str
        Print stride to COLVAR file
    plumed_file : str
        Name of file to output CV information too. (Input to PLUMED)
    colvar : str
        Name of file containing CV values over time (Output of PLUMED)

    """

    def __init__(self, print_stride='1', plumed_file='plumed.dat', 
                 colvar='COLVAR'):
        self.print_stride='1'
        self.plumed_file='plumed.dat'
        self.colvar='COLVAR'
        self.string=''
        return

    
    def add_torsion(self, atoms):
        """Adds every combination of torsion CVs given a set of atoms
        
        Parameters
        ----------
        atoms : list of integers
            Given a list of important atoms, function will generate a plumed entry
            for each set of atoms

        Example: x.add_torsion=([1,2,3,4])
        """
        # Create lists of 4 to describe dihedrals
        quads = list(combinations(atoms, 4))
        if len(atoms) >= 4:
            self.string+='\n# TORSION \n'
            for k_counter,k in enumerate(quads):
                self.string+='t{}: TORSION ATOMS={},{},{},{}\n'.format(
                        k_counter,k[0], k[1], k[2], k[3])
        else:
            print("Must include at least 4 atoms to calculate torsion")


    def add_angles(self, atoms):
        """Adds every combination of angle CVs given a set of atoms

        Parameters
        ----------
        atoms : list of integers
            Given a list of important atoms, function will generate a plumed 
            entry for each set of atoms

        Example: x.add_angles=([1,2,3,4])
        """
        # Create lists of 3 to describe angles
        triples = list(combinations(atoms, 3))
        if len(atoms) >= 3:
            self.string+='\n# ANGLES \n'
            for j_counter,j in enumerate(triples):
                self.string+='a{}: ANGLE ATOMS={},{},{}\n'.format(j_counter,
                        j[0], j[1], j[2])
        else:
            print("Must include at least 3 atoms to calculate angles")


    def add_distances(self, atoms,diff_of_dist=True, periodic='NO'):
        """Adds distances and difference of distances for every combination
        of atoms given a set of atoms

        Parameters
        ----------
        atoms : list of integer
            Given a list of important atoms, function will generate a plumed 
            entry for each set of distances
        diff_of_dist : bool 
            If set to True this function Will also calculate the difference 
            of distances, useful for likelihood maximization 
        periodic : str
            Are the CVs periodic for a difference of CVs, should always be no 
            for distnaces

        Example: x.add_distances=([1,2,3,4])
        """
        # Create pairs of 2 to describe distance
        pairs = list(combinations(atoms, 2))
        pair_list=[] # List of all the distance CVs
        if len(atoms) >= 2:
            self.string+='\n# DISTANCES \n'
            for i_counter,i in enumerate(pairs):
                self.string+='d{}: DISTANCE ATOMS={},{}\n'.format(i_counter,
                        i[0], i[1])
                pair_list.append(i_counter)
            
        if len(pair_list) >= 2 and diff_of_dist==True:  # Check list of pairs for at least two pairs 
            pair_combos = list(combinations(pair_list,2))
            self.string+='\n# DIFFERENCE OF DISTANCES \n'
            for pair_combo_counter, pair_combo in enumerate(pair_combos):    
                self.string+='dd{}: COMBINE ARG=d{},d{} COEFFICIENTS=1.0,-1.0 ' \
                        'PARAMETERS=0.0,0.0 POWERS=1.0,1.0 PERIODIC={}\n'.format(
                        pair_combo_counter, pair_combo[0], pair_combo[1], periodic)
        else:
            print("Must include at least two atoms in the atom namelist to " \
                  "generate distance CVs")


    def add_committor(self, basins, basin_CV, limit):
        """
        Parameters
        ----------
        basins : int
            Number of basins used to monitor COMMITTOR

        basin_CVs : list of strings
            List of variables that correspond to the CVs to be monitored for each
            basin. Only accepts 2-3 CVs as an input.
            Example: ['d1','d2']

        limit : array
            Array of integers that correspond to the lower and upper limits for
            each basin. Must be enetered in a specific format
        
        Example for 2 basins:
            [[LL1_x,LL1_y,UL1_x,UL1_y],[LL2_x,LL2_y,UL2_x,UL2_y]]
        * where x and y correspond to those values for each supplied argument
        * 1 and 2 correspond to each basin
        * LL = lower limit
        * UL = upper limit
        * For 3 basins, provide 3 basin constraints 
        """
        self.string+='\nCOMMITTOR ...\n'
        if basins == 2:
            self.string+='ARG={},{}\n'.format(basin_CV[0], basin_CV[1])
            self.string+='STRIDE={}\n'.format(self.print_stride)
            self.string+='BASIN_LL1={},{}\n'.format(limit[0][0], limit[0][1])
            self.string+='BASIN_UL1={},{}\n'.format(limit[0][2], limit[0][3])
            self.string+='BASIN_LL2={},{}\n'.format(limit[1][0], limit[1][1])
            self.string+='BASIN_UL2={},{}\n'.format(limit[1][2], limit[1][3])
        if basins == 3:
            self.string+='ARG={},{}\n'.format(basin_CV[0], basin_CV[1], basin_CV[3])
            self.string+='STRIDE={}\n'.format(self.print_stride)
            self.string+='BASIN_LL1={},{}\n'.format(limit[0][0], limit[0][1])
            self.string+='BASIN_UL1={},{}\n'.format(limit[0][2], limit[0][3])
            self.string+='BASIN_LL2={},{}\n'.format(limit[1][0], limit[1][1])
            self.string+='BASIN_UL2={},{}\n'.format(limit[1][2], limit[1][3])
            self.string+='BASIN_LL3={},{}\n'.format(limit[2][0], limit[2][1])
            self.string+='BASIN_UL3={},{}\n'.format(limit[2][2], limit[2][3])
        if len(basin_CV) >= 4:
            print("Use of more than 3 CVs to monitor COMMITTOR is not currently " \
                  "supported you may manually edit the plumed file")
        if len(basin_CV) < 2:
            print("You must provide at least two CVs to monitor if the shooting " \
                  "point has committed to a basin")
        self.string+='... COMMITTOR\n'

    #def add_coordination_numbers():

    def write_file(self):
        f = open(self.plumed_file, "w+")
        f.write('# vim:ft=plumed\n')
        f.write(self.string)
        f.write('\nPRINT STRIDE={} ARG=* FILE={}'.format(self.print_stride,
                self.colvar))
        f.close()


# More CVs to add

# coordination numbers -> require the input of a pdb? and grep or use newer version of plumed
# cremer pople
# Hill-Reilly?
# Berces et. al36 parameter 1

from itertools import combinations

class GeneratePlumed:
    """
    Class instantiates a plumed file and has functionality to add candidate 
    CVs such as distances, differences of distances, angles, dihedrals, etc.

    Parameters
    ----------
    print_stride : str
        Print stride to COLVAR file
    colvar : str
        Name of file containing CV values over time (Output of PLUMED)

    """

    def __init__(self, print_stride='1', plumed_file='plumed.dat', 
                 colvar='COLVAR'):
        self.print_stride='1'
        self.colvar='COLVAR'
        self.string=''
        return

    
    def add_torsion(self, atoms, periodic=False):
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
                self.string+='t{}: TORSION ATOMS={},{},{},{} '.format(
                        k_counter,k[0], k[1], k[2], k[3])
                if periodic == False:
                    self.string+='NOPBC\n'
                elif periodic == True:
                    self.string+='\n'
                    print("Periodic distances not supported in this release")
                else:
                    self.string+='\n'
                    print("must supply a bool for periodic")
        else:
            print("Must include at least 4 atoms to calculate torsion")


    def add_angles(self, atoms, periodic=False):
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
                self.string+='a{}: ANGLE ATOMS={},{},{} '.format(j_counter,
                        j[0], j[1], j[2])
                if periodic == False:
                    self.string+='NOPBC\n'
                elif periodic == True:
                    self.string+='\n'
                    print("Periodic distances not supported in this release")
                else:
                    self.string+='\n'
                    print("must supply a bool for periodic")
        else:
            print("Must include at least 3 atoms to calculate angles")


    def add_distances(self, atoms,diff_of_dist=True, periodic=False, components=False):
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
        periodic : bool
            Are the CVs periodic for a difference of CVs
        components: bool
            Option to calculate scaled compondents of the distance CV

        Example: x.add_distances=([1,2,3,4])
        """
        # Create pairs of 2 to describe distance
        pairs = list(combinations(atoms, 2))
        pair_list=[] # List of all the distance CVs
        if len(atoms) >= 2:
            self.string+='\n# DISTANCES \n'
            for i_counter,i in enumerate(pairs):
                self.string+='d{}: DISTANCE ATOMS={},{} '.format(i_counter,
                        i[0], i[1])
                if periodic == False:
                    self.string+='NOPBC'
                elif periodic == True: 
                    print("Periodic distances not supported in this release")
                else: 
                    print("must supply a bool for periodic")
                if components == False:
                    self.string+='\n'
                elif components == True:
                    self.string+='\n'
                    print("Components not supported in this release")
                else:
                    self.string+='\n'
                    print("must supply a bool for components")
                pair_list.append(i_counter)
        else:
            print("Must include at least two atoms in the atom namelist to \
                  generate distance CVs")
        
        # Difference of distances
        if len(pair_list) >= 2 and diff_of_dist==True:  # Check list of pairs for at least two pairs 
            pair_combos = list(combinations(pair_list,2))
            self.string+='\n# DIFFERENCE OF DISTANCES \n'
            for pair_combo_counter, pair_combo in enumerate(pair_combos):    
                self.string+='dd{}: COMBINE ARG=d{},d{} COEFFICIENTS=1.0,-1.0 ' \
                        'PARAMETERS=0.0,0.0 POWERS=1.0,1.0 '.format(
                        pair_combo_counter, pair_combo[0], pair_combo[1])
                if periodic == False:
                    self.string+='PERIODIC=NO\n'
                elif periodic == True:
                    self.string+='\n'
                    print("Periodic distances not supported in this release")
                else:
                    self.string+='\n'
                    print("must supply a bool for periodic")
        elif len(pair_list) >= 2 and diff_of_dist==False:
            print("Difference of distance CV not generated")
        else:
            print("Difference of distance CV not generated, requires at least 2 distance CVs")

    def add_group(self, group_name, atoms):
        """
        Parameters
        ----------
        Define groups of atoms
        group_name: str
            Name of the list of atoms you provide
            Ex: 'wat' 
        atoms: list of str
            Atoms you would like to include in the group 
            Ex: ['1,5,6,9','1-100']
        """
        if len(group_name) >= 1:
            self.string+='\n#GROUPS OF ATOMS \n'
        if len(group_name) == len(atoms):
            for group_counter, group in enumerate(group_name):
                self.string+='{}: GROUP ATOMS={}\n'.format(group, atoms[group_counter])
        else: 
            print ("Length of group names must equal length of atoms")
            
    def add_coordination(self, group_a, group_b, r_0, nl_cutoff=[1.0], nlist=True,
                         nn=6, mm=0, d_0=0.0):
        """
        Parameters
        ----------
        groupa: str
            List of atoms used for group A, easiest to define a group with
            add_group() then use the group_name as group_a
        groupb: str
        r_0 : int 
            Parameter of switching function 
            (see https://plumed.github.io/doc-v2.4/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html)
        nl_cutoff : int
            Cuttoff for the neighborlist
        nlist : bool
            Option to turn on neighbor lists (will speed up calculation, should be on) 
        Ex: 
            # first define a group of atoms 
            add_group(group1, atoms=[1,100], True)
            add_group(group2, atoms=[5,8,10], False)
            # next add coordination of that group
            add_coordination(group1, group2, 0.3, 0.6, 500)

        """

        nl_stride = self.print_stride

        if len(group_a) == len(group_b) == len(r_0) == len(nl_cutoff):
            self.string+='\n#COORDINATION \n'
            for a_counter, a in enumerate(group_a):
                self.string+='c{}: COORDINATION GROUPA={} GROUPB={} R_0={} NN={} MM={} D_0={} '.format(
                        a_counter, a, group_b[a_counter], r_0[a_counter],
                        nn, mm, d_0)
                if nlist == False:
                    self.string+='\n'
                elif nlist == True:
                    self.string+='NLIST NL_CUTOFF={} NL_STRIDE={}\n'.format(nl_cutoff[a_counter],
                            nl_stride)
                else:
                    print ('nlist must be a bool')
        else:
            print ("group_a, group_b, r_0, nl_cutoff lists must all be the same length")

    def add_puckering(self, atoms):
        """
        Given a list of 5 or 6 atoms, or a group_name, will create a cv for puckering ring

        Parameters
        ----------
        atoms: list of int
            For a 5-member ring, the five atoms should be provided as:
            [C4,O4,C1,C2,C3]
        """

        self.string+='\n#PUCKERING \n'
        for atom_list_counter, atom_list in enumerate(atoms):
            if len(atom_list) == 1:
                if isinstance(atom_list[0],str):
                    self.string+='puck{}: PUCKERING ATOMS={}\n'.format(
                        atom_list_counter, atom_list[0])
                else:
                    print("Provide a list of 5 or 6 atoms, or a group name")
            elif len(atom_list) == 5:
                self.string+='puck{}: PUCKERING ATOMS={},{},{},{},{}\n'.format(
                        atom_list_counter,atom_list[0],atom_list[1],atom_list[2],
                        atom_list[3],atom_list[4])
            elif len(atom_list) == 6:
                self.string+='puck{}: PUCKERING ATOMS={},{},{},{},{},{}\n'.format(
                        atom_list_counter,atom_list[0],atom_list[1],atom_list[2],
                        atom_list[3],atom_list[4], atom_list[5])
            else:
                print("Atom list must contain 5 or 6 atoms")

    def add_committor(self, basins, basin_CV, limit):
        """
        Parameters
        ----------
        cvs_monitored : int 
            Number of CVs used to monitor basins
        basin_CV : list of strings
            List of variables that correspond to the CVs to be monitored for each
            basin. Only accepts 2-3 CVs as an input.
            Example: ['d1','d2']

        limit : array
            Array of integers that correspond to the lower and upper limits for
            each basin. Must be enetered in a specific format
        
        Example for 2 basins and 1CV:
            [[LL1,UL1],[LL2,UL2]]
        Example for 2 basins and 2CVs:
            [[LL1_x,LL1_y,UL1_x,UL1_y],[LL2_x,LL2_y,UL2_x,UL2_y]]
        * where x and y correspond to those values for each supplied argument
        * 1 and 2 correspond to each basin
        * LL = lower limit
        * UL = upper limit
        * For 3 basins, provide 3 basin constraints 
        """
        self.string+='\nCOMMITTOR ...\n'
        if basins == 2:
            if len(basin_CV) == 1:
                self.string+='ARG={}\n'.format(basin_CV[0])
                self.string+='BASIN_LL1={}\n'.format(limit[0][0])
                self.string+='BASIN_UL1={}\n'.format(limit[0][1])
                self.string+='BASIN_LL2={}\n'.format(limit[1][0])
                self.string+='BASIN_UL2={}\n'.format(limit[1][1])
            elif len(basin_CV) == 2:
                self.string+='ARG={},{}\n'.format(basin_CV[0], basin_CV[1])
                self.string+='BASIN_LL1={},{}\n'.format(limit[0][0], limit[0][1])
                self.string+='BASIN_UL1={},{}\n'.format(limit[0][2], limit[0][3])
                self.string+='BASIN_LL2={},{}\n'.format(limit[1][0], limit[1][1])
                self.string+='BASIN_UL2={},{}\n'.format(limit[1][2], limit[1][3])
            elif len(basin_CV) >= 3:
                print("Current version not supported to monitor more than one CV")
            else:
                print("Must supply at least one CV to monitor")
        elif basins == 3:
            if len(basin_CV) == 1:
                self.string+='ARG={}\n'.format(basin_CV[0])
                self.string+='BASIN_LL1={}\n'.format(limit[0][0])
                self.string+='BASIN_UL1={}\n'.format(limit[0][1])
                self.string+='BASIN_LL2={}\n'.format(limit[1][0])
                self.string+='BASIN_UL2={}\n'.format(limit[1][1])
                self.string+='BASIN_LL3={}\n'.format(limit[2][0])
                self.string+='BASIN_UL3={}\n'.format(limit[2][1])
            elif len(basin_CV) == 2:
                self.string+='ARG={},{}\n'.format(basin_CV[0], basin_CV[1])
                self.string+='BASIN_LL1={},{}\n'.format(limit[0][0], limit[0][1])
                self.string+='BASIN_UL1={},{}\n'.format(limit[0][2], limit[0][3])
                self.string+='BASIN_LL2={},{}\n'.format(limit[1][0], limit[1][1])
                self.string+='BASIN_UL2={},{}\n'.format(limit[1][2], limit[1][3])
                self.string+='BASIN_LL3={},{}\n'.format(limit[2][0], limit[2][1])
                self.string+='BASIN_UL3={},{}\n'.format(limit[2][2], limit[2][3])
            elif len(basin_CV) >= 3:
                print("Current version not supported to monitor more than one CV")
            else:
                print("Must supply at least one CV to monitor")
        elif basins >= 4:
            print("Ability to monitor more than 3 basins is not currently supported")
        else:
            print("Must supply at least two basins to monitor")
        self.string+='STRIDE={}\n'.format(self.print_stride)
        self.string+='... COMMITTOR\n'


    def write_file(self, plumed_file):
        """
        Takes all the information built in self.string and outputs to file
        in plumed format 
        
        Parameters
        ----------
        plumed_file : str
            Name of file to output CV information too. (Input to PLUMED)
        """
        f = open(plumed_file, "w+")
        f.write('# vim:ft=plumed\n')
        f.write(self.string)
        f.write('\nPRINT STRIDE={} ARG=* FILE={}'.format(self.print_stride,
                self.colvar))
        f.close()


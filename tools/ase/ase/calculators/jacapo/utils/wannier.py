""" class Wannier provides the function get_wannier_localization_matrix and initial_wannier.
"""

import numpy as np
from Scientific.IO.NetCDF import NetCDFFile as netCDF

def dagger(matrix,copy=1):
        # First change the axis: (Does not allocate a new array)
        matrix_conj=np.swapaxes(matrix,0,1)
        if copy: # Allocate space for new array
                return np.conjugate(matrix_conj)
        else:    # The array of matrix is used for output
                np.dot(matrix_conj.imag,-1,matrix_conj.imag)
                return matrix_conj

def project(a,b):
    """ returns the projection of b onto a"""
    return a*(np.dot(np.conjugate(a),b)/np.dot(np.conjugate(a),a))

def translate(array,translation):
        """Method for translating an array"""
        newarray=array
        size=array.shape
        for dim in range(len(translation)):
            axis=dim-len(translation)
            newarray=np.concatenate((np.take(newarray,range(translation[dim],size[axis]),axis),np.take(newarray,range(translation[dim]),axis)),axis)
        return newarray.copy()

def cartesian2scaled(basis,cart):
    return np.linalg.solve(np.transpose(basis),cart)

def scaled2cartesian(basis,scaled):
    return np.dot(scaled,basis)

class Translation_Operator:
    def __init__(self,dimensions,basis):
        self.set_dimensions(dimensions)
        self.set_basis(basis)

    def set_dimensions(self,dim):
        self.dimensions = dim

    def get_dimensions(self):
        return self.dimensions
    
    def set_basis(self,basis):
        self.basis = basis

    def get_basis(self):
        return self.basis

    def set_cartesian_translation_coordinates(self,trans_cartes):
        self.cartesian_translation_coordinates = trans_cartes

    def get_cartesian_translation_coordinates(self):
        return self.cartesian_translation_coordinates

    def get_coordinates(self):
        trans_cartes = self.get_cartesian_translation_coordinates()
        #trans_coor = np.linalg.solve(np.transpose(self.get_basis()),trans_cartes) 
        trans_coor = cartesian2scaled(self.get_basis(),trans_cartes)
        return np.multiply(trans_coor,self.get_dimensions())

    def get_translational_diagonal(self,dim_index):
        # Using that G = 2*pi*n/L : epx(-i*G*x_{1}) =
        #                pow(exp(-i*2*pi/N*x_{1}),n)

        length = self.get_dimensions()[dim_index]
        coordinates = self.get_coordinates()[dim_index]
        basis_neg = np.arange(0,-length/2,-1)
        basis_pos = np.arange(length/2,0,-1)
        basis = np.concatenate((basis_neg,basis_pos),-1)
        prefactor=np.exp(-complex(0,1)*2*np.pi*coordinates/length)
        translation=map(lambda x,prefactor=prefactor:prefactor**x,basis)
        return np.array(translation)

    def operate(self,state):
        translation = np.multiply.outer(self.get_translational_diagonal(0),self.get_translational_diagonal(1))
        translation = np.multiply.outer(translation,self.get_translational_diagonal(2))
        return np.multiply(translation,state)       

def coordinate_array_from_unit_vectors(shape, gridunitvectors,
                                       origin=[0, 0, 0], indexfunction=None):
        """
        This method can be used to obtain an array representing the coordinates
        of a space defined by 'gridunitvecors'. 'gridunitvectors' is in turn a
        list containing the vectors defining the cells of the grid, i.e. the
        vectors between neighboring grid points. These vectors are spanned
        according to the specified shape. 

        'origin' -- specifies the origin of the returned coordinate array. 
        
        'indexfunction' -- is a lambda expression that defines the indices 
        with which each of the specified gridunitvectors are to be multiplied. 
        'indexfunction' must take two arguments, 'i' and 'length' - default
        is 'lambda i,length:i'. During exection the input index 'i' will run 
        over the interval 0,1,..., 'length' -1.

        **An Example**

        To obtain a coordinate array of shape (10,10) with 
        'gridunitvectors' =[[2,0],[0,1]] and the origin at [10,0] use:

        'CoordinateArrayFromUnitVectors((10,10),[[2,0],[0,1],[10,0])'

        Note that the output array will be of shape 
        (< *dimension* > ,  < *spatialcoordinates* >).
        """
        
        if indexfunction is None:
            indexfunction = lambda i, length: i
        
        coordinatelist=[]
        gridunitvectors=np.asarray(gridunitvectors)
        # Looping over the dimensionality of the vectors
        for dim in range(gridunitvectors.shape[1]):
                coordinates=origin[dim]
                # Contribution from each unitvector
                for nunitvector in range(gridunitvectors.shape[0]):
                        # Finding the indices from which the coordinate grid
                        # is spanned
                        indices=map(lambda i,f=indexfunction,l=shape[nunitvector]:f(i,l),range(shape[nunitvector]))
                        coordinatefunc=lambda i,v=gridunitvectors[nunitvector,dim]:i*v
                        coordinates=np.add.outer(coordinates,map(coordinatefunc,indices))
                coordinatelist.append(coordinates)
        return np.array(coordinatelist)

def coordinates_from_function_values(dimensions,functionvalues):
        # In general: Using basis functions of the form:
        # 1/sqrt(N1*N2*N3)*exp(iG_dot_r)
        normalization=np.sqrt(np.multiply.reduce(dimensions))
        # Functionvalues of the form (a,x_i,x_j,x_k), where
        # a is different functions and x_i,x_j,x_k the function values
        # Return array of the form (a,G_i,G_j,G_k)

        # The three last axes are transformed
        coordinates=np.fft.ifft(np.fft.ifft(np.fft.ifft(functionvalues,n=functionvalues.shape[-1],axis=-1),n=functionvalues.shape[-2],axis=-2),n=functionvalues.shape[-3],axis=-3)
        # Scaling:
        np.multiply(coordinates,normalization,coordinates)
        return coordinates


class Wannier:

   def __init__(self,calc):

      self.calc = calc
      self.set_bands(None)
      self.set_spin(None)
      self.has_changed = True

   def set_spin(self,spin): 
      self.spin = spin
      self.has_changed = True

   def get_spin(self): 
      return self.spin

   def get_fft_index(self):
      if not hasattr(self,'fftindex'):
          fftindex = []
          kpoints = self.calc.get_bz_k_points()
          for kpt in range(len(kpoints)):
              fftindex.append(self.calc.get_reciprocal_fft_index())
          self.fftindex = fftindex
      return self.fftindex

   def get_grid_dimensions(self):
      fftgrids = self.calc.get_fftgrid()
      return fftgrids['soft']
      
   def get_list_of_wave_functions(self):
      if self.get_spin() == None or self.get_bands() == None:
          raise "Bands and spin must be set before wave function list can be created"

      if self.has_changed:
          listofwavefct = []
          kpoints = self.calc.get_bz_k_points()
          for kpt in range(len(kpoints)):
              eigenstates = []
              for band in range(self.get_bands()):
                  eigenstates.append(self.calc.get_reciprocal_bloch_function(band=band,kpt=kpt,spin=self.get_spin()))

              listofwavefct.append(eigenstates)
          self.listofwavefct = listofwavefct
          self.has_changed = False
      return self.listofwavefct

   def set_bands(self,numberofbands):	
      self.numberofbands = numberofbands
      self.has_changed = True
      
   def get_bands(self): 
      return self.numberofbands

   def get_zi_bloch_matrix(self,dirG,kpoint,nextkpoint,G_I):
      """ calculate matrix of ZIi,j values
      This matrix consist of 3 matrices each of dimension MxM, i.e. corresponding to the full space.
      """

      #print 'DacapoWannier: Initialize ZIBlochMatrix ..'
      if self.get_bands() == None or self.get_spin() == None:
          raise "Bands and spin must be set before wannier localization matrix can be calculated"
     	
      phi=np.swapaxes(np.array(self.get_list_of_wave_functions()[kpoint]),0,1)
      # K1 and reciprocal lattice vector G_I  given kpoint K
      # that fulfills the criteria : K1-K-K0+G1=0
      list1,list2 = self.get_gg_list(kpoint,nextkpoint,G_I)
				
      a=np.take(phi,list1,axis=0)
      a=np.swapaxes(a,0,1)
      phi1 = np.swapaxes(np.array(self.get_list_of_wave_functions()[nextkpoint]),0,1)
      b=np.take(phi1,list2,axis=0)

      ziblochmatrix = np.dot(np.conjugate(a),b)
      usziblochmatrix = self.get_ultra_soft_non_loc_matrix(dirG,kpoint,nextkpoint)
      ziblochmatrix += usziblochmatrix

      return ziblochmatrix


   def get_gg_list(self,kpt1,kpt2,GI):
       """ define list of (G,G+G1) defining the product
       phi(kpt1,G)*phi(kpt2,G+G1),

       GI is one of
       [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]]
			 
       The layout of fourier components is 
       1   2   3   4   5   6   7   8   ngx = 8 
       0   1   2   3   4  -3  -2  -1	n*2pi/L	
       """

       numberplanewaves = len(self.get_list_of_wave_functions()[kpt1][0])
       reciprocalindex = self.get_fft_index()[kpt1]

       ngrids = self.get_grid_dimensions()
			
       # setup the mapping from the 3D FFT grid to the wavefuction list
       map2 = self.get_index_map(kpt2) 
		
       gglist = []
       # print "Generating plane wave index list for direction ",GI," kpt12 ",kpt1,kpt2
       list1 = []
       list2 = []
	   
       # find G,G+GI
       for n in range(numberplanewaves):
           index = reciprocalindex[:,n]
	   index = index - 1
	   # find G,G+GI
	   for dir in range(3): 
	       index[dir] += GI[dir]
	       if index[dir]>=ngrids[dir]:
	           # wrap around
		   index[dir] = 0
				
           # now find the corresponding index into phi(kpt2)
           n1 = map2[index[0],index[1],index[2]]

           if n1>=0:
              list1.append(n)
              list2.append(n1)
			
       # print '  Number of elements in GG list ',len(list1)
       return list1,list2

		

   def get_index_map(self,kpt):
       """ generate mapping from 3D FFT grid to the wavefunction list
		
       A negative number is returned from map(g1,g2,g3) is the
       grid point does not exists in the wavefunction list
       """
       ngrids = self.get_grid_dimensions()
       map_to_wflist = np.zeros(ngrids,np.int)
       map_to_wflist = map_to_wflist - 1

       numberplanewaves = len(self.get_list_of_wave_functions()[kpt][0])
       reciprocalindex = self.get_fft_index()[kpt]
			
       for n in range(numberplanewaves):
           i0 = reciprocalindex[0][n]-1
	   i1 = reciprocalindex[1][n]-1
	   i2 = reciprocalindex[2][n]-1
	   map_to_wflist[i0,i1,i2] = n

       return map_to_wflist

   def get_ultra_soft_non_loc_matrix(self,GI,kpt,kpt1):
        """ calculate
              a                            I        I              I
             W    = sum(n,m,I) <psi  | beta  > <beta  | psi   > * q
               i,j                  ik      m        n      jk1     mn
 
             n,m : projectors
             I   : atom no
             a (nbands,nbands) matrix is returned.
        """

        nc = netCDF(self.calc.get_nc(),'r')
        if not 'WannierAugFactor' in nc.variables:
            nc.close()
            return

        # Read NLProjectorPsi, StructureFactor, WannierAugFactor
        vnlprojpsi = nc.variables['NLProjectorPsi'][:]
        nlprojpsi = np.zeros((np.shape(vnlprojpsi[:,:,:,:,:,0])),np.complex)
        nlprojpsi.real = vnlprojpsi[:,:,:,:,:,0]
        nlprojpsi.imag = vnlprojpsi[:,:,:,:,:,1]
        nlprojpsi = np.swapaxes(nlprojpsi,3,4)

        vstrfactor = nc.variables['StructureFactor'][:]
        strfactor = np.zeros((np.shape(vstrfactor[:,:,0])),np.complex)
        strfactor.real = vstrfactor[:,:,0]
        strfactor.imag = vstrfactor[:,:,1]

        vaugfactor = nc.variables['WannierAugFactor'][:]
        augfactor = np.zeros((np.shape(vaugfactor[:,:,:,:,0])),np.complex)
        augfactor.real = vaugfactor[:,:,:,:,0]
        augfactor.imag = vaugfactor[:,:,:,:,1]
        augfactor = np.swapaxes(augfactor,0,1)

        nc.close()

        lst=[[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]]

        # find direction corresponding to GI
        dir = lst.index(GI.tolist()) 

        natoms = nlprojpsi.shape[2]
        nbands = self.get_bands()
        matrix = np.zeros([nbands,nbands],np.complex)
        for atom in range(natoms):

            if kpt==kpt1:
                q=np.conjugate(augfactor[:,:,atom,dir+1])
            else:
                q=augfactor[:,:,atom,0]

                if dir<3:
                    q=q*np.conjugate(strfactor[atom,dir])

            A=nlprojpsi[kpt,self.get_spin(),atom,:nbands,:]
            B=nlprojpsi[kpt1,self.get_spin(),atom,:nbands,:]
            matrix=matrix+np.dot(A,np.dot(q,dagger(B)))
        return matrix

##### code to calculate initial wannier functions ######

   def set_data(self,data):
        self.data = data

   def set_k_point_grid(self,kpointgrid):
        self.kpointgrid = np.array(kpointgrid)

   def get_k_point_grid(self):
        return self.kpointgrid

   def get_repeated_unit_cell(self):
        basis=self.calc.get_atoms().get_cell()
        return np.transpose(np.transpose(basis)*self.get_k_point_grid())

   def get_repeated_grid_dimensions(self):
        return self.get_grid_dimensions()*self.get_k_point_grid()

   def get_detailed_data(self):
        if not hasattr(self,'detaileddata'):
            """ Coverts data on the form
                [atom,l,(m),a], or [center,l,(m),a].
                atom is a number in a bracket, example [2] for atom number 2
                center is 3 numbers in a bracket, example [1,0.5,0.3] denoting SCALED coordinates of the center.
                Here m is optional. If m is left out all m values are used."""
            datalist=self.data
            detaileddata=[]
            atoms=self.calc.get_atoms()
            for data in datalist:
                if len(data[0])==1:
                    r_c=atoms[data[0][0]].get_position()
                elif len(data[0])==3:
                    r_c = scaled2cartesian(atoms.get_cell(),data[0])
                else:
                    print "First element in initial data must be of the form [atom] or [c1,c2,c3], where the latter is scaled coordinates of the center"
                if len(data)==4:
                    # m is specified        
                    detaileddata.append([r_c,data[1],data[2],data[3]])
                else:
                    # Orbitals with all allowed m values are produced
                    for m in range(-data[1],data[1]+1):
                        detaileddata.append([r_c,data[1],m,data[2]])
            self.detaileddata = detaileddata

        return self.detaileddata

   def get_origin_index(self):

        griddim = self.get_repeated_grid_dimensions()
        originindex = map(lambda coor:int(coor),list(np.multiply(np.array([0.5,0.5,0.5]),griddim)))
        return originindex

   def get_cartesian_coordinates(self):

        griddim = self.get_repeated_grid_dimensions()
        basis = self.calc.get_atoms().get_cell()
         
        # Set origin to center of grid
        originindex = self.get_origin_index()
        # origincoord is in scaled coordinates:
        origincoord=-np.array(np.array(originindex,np.float)/griddim)
        origincart = scaled2cartesian(basis,origincoord)
        gridunitvectors = np.array(map(lambda unitvector,shape:unitvector/shape,basis,griddim))
        c = coordinate_array_from_unit_vectors(shape=griddim,
                                           gridunitvectors=gridunitvectors,
                                           origin=origincart)
        return c

   def get_normalized_coordinates(self):
 
        if not hasattr(self,'normalized_coordinates'):
            originindex = tuple(self.get_origin_index())
            c = self.get_cartesian_coordinates()
            dist = np.sqrt(c[0]**2+c[1]**2+c[2]**2)

            # We define "normalized" coordinates. To avoid undeterminancy at origin we move
            # to the point (1,1,1)*1e-8
            c[0][originindex]=1.0e-8
            c[1][originindex]=1.0e-8
            c[2][originindex]=1.0e-8

            dist[originindex]=np.sqrt(3)*1.0e-8
            self.normalized_coordinates = c/dist

        return self.normalized_coordinates

   def get_distance_array_at_origin(self):

        # Set origin to center of grid
        originindex = self.get_origin_index()

        c = self.get_cartesian_coordinates()
        dist=np.sqrt(c[0]**2+c[1]**2+c[2]**2)
        # Translate back to origin
        dist=translate(dist,originindex)
        return dist

   def setup_m_matrix(self,listofeigenstates,bzkpoints):

        if self.data == None or self.kpointgrid == None:
            raise "Must set data, kpointgrid, spin before calculating M matrix"

        fftindex = self.get_fft_index()

        unitcell=self.get_repeated_unit_cell()
        griddim = self.get_repeated_grid_dimensions()    
        data = self.get_detailed_data()
        nkpoints = len(bzkpoints)
        nbands = len(listofeigenstates[0])
        M = np.zeros([nkpoints,nbands,len(data)],np.complex)
        orbital = np.zeros(griddim,np.complex)
        dist = self.get_distance_array_at_origin() 

        transop = Translation_Operator(griddim,unitcell)
        rec_basis = 2*np.pi*np.linalg.inv(np.transpose(self.calc.get_atoms().get_cell()))
        large_rec_basis = 2*np.pi*np.linalg.inv(np.transpose(unitcell))

        for i in range(len(data)):
            # Translate orbital
            r_c=data[i][0]
            l,m=data[i][1],data[i][2]
            a=data[i][3]
            orbital = self.get_cubic_harmonic_at_origin(l,m)*np.exp(-dist/a)
            orbital_fft = coordinates_from_function_values(griddim,orbital)
            transop.set_cartesian_translation_coordinates(r_c)
	    orbital_fft = transop.operate(orbital_fft)
            for kpt in range(nkpoints):
                kpoint = bzkpoints[kpt]
                kptnumber = cartesian2scaled(large_rec_basis,scaled2cartesian(rec_basis,kpoint))
                kptnumber[0]=round(kptnumber[0])
                kptnumber[1]=round(kptnumber[1])
                kptnumber[2]=round(kptnumber[2])
                kptnumber=kptnumber.astype(int)
                u_k = self.extract_periodic_part_of_small_cell(orbital_fft,kptnumber)
                compact_u_k = self.get_compact_fft_representation(u_k,fftindex[kpt],len(listofeigenstates[kpt][0]))
                M[kpt,:,i] = np.dot(np.conjugate(np.array(listofeigenstates[kpt])),compact_u_k)

        self.mmatrix = M

   def get_cubic_harmonic_at_origin(self,l,m):
       """ l=0,1,2. m=-l,...,l"""

       griddim = self.get_repeated_grid_dimensions()
       harmonic = np.zeros(griddim,np.complex)
       originindex = self.get_origin_index()
       nc = self.get_normalized_coordinates()

       # Constructing cubic harmonic
       if l==0 and m==0:
           harmonic=(1/np.sqrt(4*np.pi))*np.ones(nc[0].shape,np.Complex)
           harmonic=translate(harmonic,originindex)
       if l==1 and m==0:
           # p_x
           harmonic=np.sqrt(3/(4*np.pi))*nc[0]
           harmonic=translate(harmonic,originindex)
       if l==1 and m==-1:
           # p_z
           harmonic=np.sqrt(3/(4*np.pi))*nc[2]
           harmonic=translate(harmonic,originindex)
       if l==1 and m==1:
           # p_y
           harmonic=np.sqrt(3/(4*np.pi))*nc[1]
           harmonic=translate(harmonic,originindex)
       if l==2 and m==0:
           harmonic=0.5*np.sqrt(5/(4*np.pi))*(3*(nc[0]**2)-np.ones(nc[0].shape,np.Complex))
           harmonic=translate(harmonic,originindex)
       if l==2 and m==-1:
           harmonic=np.sqrt(15/(16*np.pi))*(nc[2]**2-nc[1]**2)
           harmonic=translate(harmonic,originindex)
       if l==2 and m==1:
           harmonic=np.sqrt(15/(4*np.pi))*nc[0]*nc[1]
           harmonic=translate(harmonic,originindex)
       if l==2 and m==-2:
           harmonic=np.sqrt(15/(4*np.pi))*nc[2]*nc[1]
           harmonic=translate(harmonic,originindex)
       if l==2 and m==2:
           harmonic=np.sqrt(15/(4*np.pi))*nc[2]*nc[0]
           harmonic=translate(harmonic,originindex)
       return harmonic

   def extract_periodic_part_of_small_cell(self,f,k):
       n1,n2,n3=self.get_k_point_grid()
       trans=[0,0,0]
       if k[0]<0:
           k[0]+=n1
           trans[0]=self.get_grid_dimensions()[0]-1
       if k[1]<0:
           k[1]+=n2
           trans[1]=self.get_grid_dimensions()[1]-1
       if k[2]<0:
           k[2]+=n3
           trans[2]=self.get_grid_dimensions()[2]-1

       u=f[k[0]::n1,k[1]::n2,k[2]::n3].copy()
       return translate(u,trans)

   def get_compact_fft_representation(self,freciprocal,fftindex,numberofpws):
       wflist=np.zeros([numberofpws],np.complex)
       for i in range(numberofpws):
           wflist[i]=freciprocal[int(fftindex[0,i]-1),int(fftindex[1,i]-1),int(fftindex[2,i]-1)]
       return wflist

   def get_orthonormality_factor(self,matrix):
       defect = abs(np.dot(dagger(matrix),matrix))-np.identity(matrix.shape[1],np.Float)
       return max(abs(defect.flat))

   def get_list_of_coefficients_and_rotation_matrices(self,matrixdimensions):
       from ase.dft.wannier import normalize,gram_schmidt
       import random, copy
       M,N,L=matrixdimensions
       nkpt=len(N)
       Ulist=[]
       clist=[]
       if not hasattr(self,'mmatrix'):
           raise "Must setup M Matrix first!"
       coeffmatrix=self.mmatrix
       for kpt in range(nkpt):
           #First normalize the columns of coeffmatrix
           normalize(coeffmatrix[kpt])
           T=coeffmatrix[kpt][N[kpt]:].copy()
           numberoforbitals=T.shape[1]
           c=np.zeros([M-N[kpt],L[kpt]],np.complex)
           U=np.zeros([N[kpt]+L[kpt],N[kpt]+L[kpt]],np.complex)
           # Initialize weights
           w=abs(np.sum(T*np.conjugate(T)))
           for i in range(min(L[kpt],numberoforbitals)):
               # Find index of maximal element in w
               t=w.tolist().index(max(w))
               c[:,i]=T[:,t]
               # Orthogonalize c[:,i] on previous vectors
               for j in range(i):
                   c[:,i]=c[:,i]-project(c[:,j],T[:,t])
               c[:,i]=c[:,i]/np.sqrt(np.dot(c[:,i],np.conjugate(c[:,i])))
               # Update weights
               w=w-abs(np.dot(np.conjugate(c[:,i]),T))**2
           if numberoforbitals<L[kpt]:
               # Supplement c by random vectors
               for i in range(numberoforbitals,L[kpt]):
                   for j in range(M-N[kpt]):
                       c[j,i]=random.random()
               gram_schmidt(c)
           # Test whether columns are orthonormal
           if L[kpt]>0:
               test = self.get_orthonormality_factor(c)
               if test>1.0e-3:
                   print "ERROR: Columns of c not orthogonal!"
                   
           U[:N[kpt],:numberoforbitals]=coeffmatrix[kpt][:N[kpt]]
           U[N[kpt]:,:numberoforbitals]=np.dot(dagger(c),coeffmatrix[kpt][N[kpt]:])
           # Perform democratic Lowdin orthogonalization on U[:,numberoforbitals]
           gram_schmidt(U[:,:numberoforbitals])
           if numberoforbitals<(N[kpt]+L[kpt]):
               #Supplement U by random vectors
               for i in range(numberoforbitals,N[kpt]+L[kpt]):
                   for j in range(N[kpt]+L[kpt]):
                       U[j,i]=random.random()
               # Finally orthogonalize everything
               # Note, only random vectors are affected
               gram_schmidt(U)
           # Test whether columns are orthonormal
           test = self.get_orthonormality_factor(U)
           if test>1.0e-3:
               print "ERROR: Columns of U not orthogonal for kpoint",kpt
           Ulist.append(U)
           clist.append(c)

       return clist,Ulist


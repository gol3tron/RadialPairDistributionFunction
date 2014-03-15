# RADIAL PAIR DISTRIBUTION FUNCTION MODULE
# Adam Goler, Created 14 March 2014, Last Update
# adamgoler@gmail.com
#
# This module uses tXYZ coordinate data (taken from a NAMD DCD trajectory)
# and uses that data to calculate the RPDF, using two distinct selections
# of atoms.
#
# First use makeSel(...) to create the partitioned trajectories corresponding
# to the desired selection (residue indices must be known a priori). Then use
# getPairRadii(...) to calculate the pairwise distances between atoms of the
# two selections. Finally, use plotRDF(...) to calculate the normalized
# radial distribution function and plot it.
#
#
#
#

def getPairRadii(sel1,sel2):

    import numpy as np              #import numpy for zeros() fun
    import math as m                #import math for pi, sqrt, ceil, floor, etc.

    nframes = len(sel1)
    natoms1 = len(sel1[0])          #number of atoms in sel1
    natoms2 = len(sel2[0])          #number of atoms in sel2

    natomstot = natoms1 + natoms2   #total number of atoms in sel1 + sel2
    
    rvals = []
    temp = np.zeros(3)

    xyz_CoM = [m.sqrt(temp[j])/natoms1 for j in range(3)]

    # loop over particles in sel2
    # calculate distance from sel1 CoM
    # and sort into correct histogram bin
    for t in range(nframes):
        for i in range(natoms1):
            for j in range(natoms2):
                # reset temp
                temp = [0.0,0.0,0.0]
                rptemp = 0.0
    
                for k in range(3):
                    temp[k] = sel2[t][j][k] - sel1[t][i][k]
                    rptemp += temp[k]**2
    
                rvals.append(m.sqrt(rptemp))

    return rvals,natoms1,natoms2

def plotRDF(rvals,dr,rmax,nframes,natoms1,natoms2):

    import numpy as np
    import math as m
    import matplotlib.pyplot as plt

    n = int(m.ceil(rmax/dr))
    r = range(1,n+1)
    rdf_hist = np.zeros(n)
    nvals = len(rvals)
    ntot = natoms1 + natoms2
    Nr = natoms1 * natoms2

    for i in range(n):
        for j in range(nvals):
            
            test = int(m.ceil(rvals[j]))

            if (test == r[i]):
                rdf_hist[i] += 1

        # Normalize your shit, bro!
        # See Eq. (5.1) for reference : http://homepage.univie.ac.at/franz.vesely/simsp/dx/node22.html
    vol = (4./3.)*m.pi*rmax**3
    natoms = ntot
    numberden = natoms/vol
    area = [4.*m.pi*r[i]**2 for i in range(n)]
    
    
    for i in range(n):
        # divide by spherical shell volume at r = r[i] (dr = 1)
        rdf_hist[i] /= area[i]
        # divide by total number of atoms
        rdf_hist[i] /= natoms
        # divide by number density N/V
        rdf_hist[i] /= numberden
        # dont forget to average count over all frames! only have so many particles!
        rdf_hist[i] /= nframes


    plt.plot(r,rdf_hist,label='radial distribution function')
    plt.xlabel('radius (Angstroms)')
    plt.ylabel('radial pair distribution function g(r)')

    return rdf_hist


def makeSel(allcoords,beg,end):
    
    # this method takes allcoords as input txNx3 array and returns, with selection params
    # and returns a subset of coordinates of a specific set of atoms
    
    import numpy as np
    
    nframes = len(allcoords)
    natoms = len(allcoords[0])
    dim = len(allcoords[0][0])
    section = end-beg
    
    
    # create a partitioned array for selection
    sel = np.zeros((nframes,section,dim))

    for i in range(nframes):
        for j in range(beg-1,end):
            for k in range(dim):
                sel[i][beg-j][k] = allcoords[i][j][k]

    return sel



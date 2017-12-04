################################################################################
# OCTAHEDRAL TILTING IN PEROVSKITES in Python3                                 #
#                                                                              #
# This script creates a 2x2x2 supercell from any cubic perovskite. The         #
# octahedra are tilted according to the entered Glazer tilt system.            #
#                                                                              #
# created by Michael Giger and Fabian Haake                                    #
################################################################################

# IMPORT

import re
import math

# PARAMETERS

# alat is the lattice constant for the simple cubic unit cell. Cation A is the
# name of the atom surrounded by octahedra. Cation M is the name of the atom in
# the center of octahedra. Anion X is the name of the atom which sits in the
# corner of the octahedra.
# direction is a string with three characters. Each one of them represents the
# type of tilting. Possible characters are '+', '-' and '0'. '+' stands for
# inphase tilting, '-' stands for antiphase tilting. With '0', no tilting will
# be applied. tilt gives the amount of tilt in degree for each direction.
# file is the name of the created POSCAR file.

alat  = 3.845
cationA = 'Ca'
cationM = 'Ti'
anionX  = 'O'

direction = '0+-'
tilt = [0,10,15]

file = 'POSCAR'

# FUNCTION

# Addition of two vectors.
def additionVector(v1,v2):
    v = [i+j for i,j in zip(v1,v2)]
    return v

# Product of a vector and a scalar.
def scalarVector(v,a):
    v = [a*i for i in v]
    return v

# Returns corresponding vectors for a given vector v in a 2x2x2 supercell.
def supercellVector(v):
    stv = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ta = scalarVector(a,k)
                tb = scalarVector(b,j)
                tc = scalarVector(c,i)
                tv = additionVector(v,ta)
                tv = additionVector(tv,tb)
                tv = additionVector(tv,tc)
                stv.append(scalarVector(tv,.5))
    return stv

# Tiltes group of vectors v of a supercell according to the entered tilt system.
# t determines which of the three anions X was meant. i is the number which
# represents which axis is handeld (0 => a, 1 => b, 2 => c). dN saves if it has
# to apply inphase or antiphase tilting. j is the number of which vector of the
# of the group v is handeld.
def tiltVector8(v,t):
    dirRot = [[[[-1,-1,1,1,1,1,-1,-1],[0,0,0,0,0,0,0,0],[1,1,-1,-1,-1,-1,1,1]],
               [[1,-1,1,-1,-1,1,-1,1],[-1,1,-1,1,1,-1,1,-1],[0,0,0,0,0,0,0,0]],
               [[0,0,0,0,0,0,0,0],[1,-1,-1,1,1,-1,-1,1],[-1,1,1,-1,-1,1,1,-1]]],
              [[[1,-1,-1,1,-1,1,1,-1],[0,0,0,0,0,0,0,0],[-1,1,1,-1,1,-1,-1,1]],
               [[-1,1,1,-1,1,-1,-1,1],[1,-1,-1,1,-1,1,1,-1],[0,0,0,0,0,0,0,0]],
               [[0,0,0,0,0,0,0,0],[-1,1,1,-1,1,-1,-1,1],[1,-1,-1,1,-1,1,1,-1]]]]
    
    for i in range(3):
        ra = [tilt[i],0,0]
        rb = [0,tilt[i],0]
        rc = [0,0,tilt[i]]
        absRot = [[rb,[0,0,0],rc],
                  [ra,rc,[0,0,0]],
                  [[0,0,0],rb,ra]]
            
        if direction[i] == '+':
            dN = 0
        elif direction[i] == '-':
            dN = 1
        else:
            dN = 2
        if dN != 2:
            for j in range(8):
                v[j] = additionVector(v[j],scalarVector(absRot[i][t],dirRot[dN][i][t][j]))
    return v
                
            

# Output for a vector.
def writeVector(v):
    for i in range(3):
        f.write(f'    {v[i]:18.15f}')
    f.write('\n')

# Output for 8 vectors.
def writeVector8(v):
    for i in range(8):
        writeVector(v[i])

# MAIN

# Creates basis of the crystal and corresponding atom positions in a fixed
# sequence (which is important for the tiltVector8 function. The entered tilts
# will also converted to absolute displacements.
a = [1,0,0]
b = [0,1,0]
c = [0,0,1]
 
posA  = [0,0,0]
posM  = [.5,.5,.5]
posX1 = [.5,.5,0]
posX2 = [0,.5,.5]
posX3 = [.5,0,.5]

tilt = [.25*math.tan(math.radians(i)) for i in tilt]

# Creates supercell and applies tilts.
posA  = supercellVector(posA)
posM  = supercellVector(posM)
posX1 = supercellVector(posX1)
posX2 = supercellVector(posX2)
posX3 = supercellVector(posX3)

posX1 = tiltVector8(posX1,0)
posX2 = tiltVector8(posX2,1)
posX3 = tiltVector8(posX3,2)

# Creates new file.
f = open(file,'w')
f.write('Tilted Supercell 2x2x2\n')
f.write(f'   {alat:.3f}\n')
writeVector(scalarVector(a,2))
writeVector(scalarVector(b,2))
writeVector(scalarVector(c,2))
f.write('   %2s   %2s   %2s\n' % (cationA,cationM,anionX))
f.write('    8    8   24\n')
f.write('Default\n')
writeVector8(posA)
writeVector8(posM)
writeVector8(posX1)
writeVector8(posX2)
writeVector8(posX3)
f.close()





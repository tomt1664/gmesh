#!/usr/bin/python

# atom2obj.py
# T.Trevethan
# 2015
#
# Script to create a triangle mesh from an sp2 bonded atomistic configuration
#
# The structure is supplied as an xyz format file: input.xyz
# 
# The mesh is output as output.obj a OBJ format file (Wavefront)
#
# Optional command line arguments:  
#                               
#       -mn float : minimum bond length for connection table (Default 1.0)
#       -mx float : maximum bond length for connection table (Default 1.7) 

import sys
import math 
import argparse
import numpy as np

minb = 1.0
maxb = 1.7

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-mn', type = float, help = 'min bond length')
parser.add_argument('-mx', type = float, help = 'max bond length')
variables = parser.parse_args()
if variables.mn:
    minb = variables.mn
if variables.mx:
    maxb = variables.mx

print "atom2obj: convert xyz format atomic structure to OBJ format triangle mesh"

# Open input files
input  = open("input.xyz", 'r')

#read input file and store values in array
coords = []

numat = input.readline()
nat = int(numat)
blank = input.readline()

for i in range(nat):
    line = input.readline()
    fdata = line.split()
    data = (float(fdata[1]),float(fdata[2]),float(fdata[3]))
    coords.append(data)

print "Read in "+str(nat)+" atoms"
print "Creating bond table ..."

#create bond table
bond = np.zeros((nat,4),dtype=np.int)

for i in range(nat):
    ib = 0
    for j in range(nat):
        if i != j:
            xdst = coords[i][0] - coords[j][0]
            ydst = coords[i][1] - coords[j][1]
            zdst = coords[i][2] - coords[j][2]
            dst = math.sqrt(xdst*xdst+ydst*ydst+zdst*zdst)
            if dst <= minb:
                print "Error: atoms "+str(i)+" and "+str(j)+" closer than min dist"
                sys.exit(1)
            if dst > minb and dst < maxb:
                if ib > 3:
                    print "Error: atom "+str(i)+" has more than 4 bonds"
                    sys.exit(1)
                bond[i][ib] = j
                ib += 1

#getting the bond list
nn2 = []
num2 = 0
for i in range(nat):
    for ib in range(4):
        if bond[i][ib] > i:
            bdat = (i,bond[i][ib])
            nn2.append(bdat)
            num2 += 1

print "No. bonds: "+str(num2)

#getting the 3-seq list
nn3 = []
num3 = 0
for i in range(num2):
    for j in range(i+1,num2):
        if nn2[i][0] == nn2[j][0]:
            adat = (nn2[i][1],nn2[i][0],nn2[j][1])
            nn3.append(adat)
            num3 += 1
        elif nn2[i][1] == nn2[j][0]:
            adat = (nn2[i][0],nn2[i][1],nn2[j][1])
            nn3.append(adat)
            num3 += 1
        elif nn2[i][1] == nn2[j][1]:
            adat = (nn2[i][0],nn2[i][1],nn2[j][0])
            nn3.append(adat)
            num3 += 1
        elif nn2[i][0] == nn2[j][1]:
            adat = (nn2[i][1],nn2[i][0],nn2[j][0])
            nn3.append(adat)
            num3 += 1

print "No. 3-body: "+str(num3)

#getting the 4-seq list
nn4 = []
num4 = 0
for i in range(num3):
    for j in range(i+1,num3):
        if nn3[i][1] == nn3[j][0] and nn3[i][2] == nn3[j][1]:
            tdat = (nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][0] == nn3[j][1] and nn3[i][1] == nn3[j][2]:
            tdat = (nn3[j][0],nn3[i][0],nn3[i][1],nn3[i][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][1] == nn3[j][2] and nn3[i][2] == nn3[j][1]:
            tdat = (nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][0])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][0] == nn3[j][1] and nn3[i][1] == nn3[j][0]:
            tdat = (nn3[j][2],nn3[j][1],nn3[j][0],nn3[i][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][1] == nn3[j][0] and nn3[i][0] == nn3[j][1]:
            tdat = (nn3[j][2],nn3[i][0],nn3[i][1],nn3[i][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][2] == nn3[j][1] and nn3[i][1] == nn3[j][2]:
            tdat = (nn3[j][0],nn3[i][2],nn3[i][1],nn3[i][0])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][1] == nn3[j][2] and nn3[i][0] == nn3[j][1]:
            tdat = (nn3[j][0],nn3[i][2],nn3[i][1],nn3[i][0])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][2] == nn3[j][1] and nn3[i][1] == nn3[j][0]:
            tdat = (nn3[i][2],nn3[i][1],nn3[i][0],nn3[j][2])
            nn4.append(tdat)
            num4 += 1

print "No. 4-body: "+str(num4)

#getting the 5-seq list
nn5 = []
num5 = 0
for i in range(num4):
    for j in range(i+1,num4):
        if nn4[i][1] == nn4[j][0] and nn4[i][2] == nn4[j][1] and nn4[i][3] == nn4[j][2]:
            pdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][3])
            nn5.append(pdat)
            num5 += 1
        elif nn4[i][0] == nn4[j][1] and nn4[i][1] == nn4[j][2] and nn4[i][2] == nn4[j][3]:
            pdat = (nn4[j][0],nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3])
            nn5.append(pdat)
            num5 += 1
        elif nn4[i][1] == nn4[j][3] and nn4[i][2] == nn4[j][2] and nn4[i][3] == nn4[j][1]:
            pdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][0])
            nn5.append(pdat)
            num5 += 1
        elif nn4[i][0] == nn4[j][2] and nn4[i][1] == nn4[j][1] and nn4[i][2] == nn4[j][0]:
            pdat = (nn4[j][3],nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3])
            nn5.append(pdat)
            num5 += 1

print "No. 5-body: "+str(num5)

#find all the triangles
r3 = []
numr3 = 0
for i in range(num4):
    if nn4[i][0] == nn4[i][3]:
        tdat = [nn4[i][0],nn4[i][1],nn4[i][2]]
        iapp = True
        for cmp in r3:
            if sorted(tdat) == sorted(cmp): iapp = False
        if iapp:
            r3.append(tdat)
            numr3 += 1

print "No. 3-fold rings: "+str(numr3)

#find all the squares
r4 = []
numr4 = 0
for i in range(num3):
    for j in range(i+1,num3):
        if nn3[i][0] == nn3[j][0] and nn3[i][2] == nn3[j][2]:
            sdat = [nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][1]]
            iapp = True
            for cmp in r4:
                if sorted(sdat) == sorted(cmp): iapp = False
            if iapp:
                r4.append(sdat)
                numr4 += 1
        elif nn3[i][0] == nn3[j][2] and nn3[i][2] == nn3[j][0]:
            sdat = [nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][1]]
            iapp = True
            for cmp in r4:
                if sorted(sdat) == sorted(cmp): iapp = False
            if iapp:
                r4.append(sdat)
                numr4 += 1

print "No. 4-fold rings: "+str(numr4)

#find all the pentagons
r5 = []
numr5 = 0
for i in range(num4):
    for j in range(num3):
        if nn4[i][0] == nn3[j][0] and nn4[i][3] == nn3[j][2]:
            if nn4[i][1] != nn3[j][1] and nn4[i][2] != nn3[j][1]:
                sdat = [nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn3[j][1]]
                iapp = True
                for cmp in r5:
                    if sorted(cmp) == sorted(sdat): iapp = False
                if iapp:
                    r5.append(sdat)
                    numr5 += 1
        elif nn4[i][3] == nn3[j][0] and nn4[i][0] == nn3[j][2]: 
            if nn4[i][1] != nn3[j][1] and nn4[i][2] != nn3[j][1]:
                sdat = [nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn3[j][1]]
                iapp = True
                for cmp in r5:
                    if sorted(cmp) == sorted(sdat): iapp = False
                if iapp:
                    r5.append(sdat)
                    numr5 += 1

print "No. 5-fold rings: "+str(numr5)

#find all the hexagons
r6 = []
numr6 = 0
for i in range(num4):
    for j in range(i+1,num4):
        if nn4[i][0] == nn4[j][0] and nn4[i][3] == nn4[j][3]:
            if nn4[i][1] != nn4[j][1] and nn4[i][2] != nn4[j][2]:
                hdat = [nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][2],nn4[j][1]]
                iapp = True
                for cmp in r6:
                    if sorted(cmp) == sorted(hdat): iapp = False
                if iapp:
                    r6.append(hdat)
                    numr6 += 1
        elif nn4[i][0] == nn4[j][3] and nn4[i][3] == nn4[j][0]:
            if nn4[i][2] != nn4[j][1] and nn4[i][1] != nn4[j][2]:
                hdat = [nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][1],nn4[j][2]]
                iapp = True
                for cmp in r6:
                    if sorted(cmp) == sorted(hdat): iapp = False
                if iapp:
                    r6.append(hdat)
                    numr6 += 1

print "No. 6-fold rings: "+str(numr6)

#find all the heptagons
r7 = []
numr7 = 0
for i in range(num5):
    for j in range(num4):
        if nn5[i][0] == nn4[j][0] and nn5[i][4] == nn4[j][3]:
            if nn5[i][1] != nn4[j][1] and nn5[i][3] != nn4[j][2]:
                hdat = [nn5[i][0],nn5[i][1],nn5[i][2],nn5[i][3],nn5[i][4],nn4[j][2],nn4[j][1]]
                iapp = True
                for cmp in r7:
                    if sorted(cmp) == sorted(hdat): iapp = False
                if iapp:
                    r7.append(hdat)
                    numr7 += 1
        elif nn5[i][0] == nn4[j][3] and nn5[i][4] == nn4[j][0]:
            if nn5[i][1] != nn4[j][2] and nn5[i][3] != nn4[j][1]:
                hdat = [nn5[i][0],nn5[i][1],nn5[i][2],nn5[i][3],nn5[i][4],nn4[j][1],nn4[j][2]]
                iapp = True
                for cmp in r7:
                    if sorted(cmp) == sorted(hdat): iapp = False
                if iapp:
                    r7.append(hdat)
                    numr7 += 1

print "No. 7-fold rings: "+str(numr7)
print "Triangulating polygons ..."

#get triangles from: squares
for sqr in r4:
    tdat = [sqr[0],sqr[1],sqr[2]]
    r3.append(tdat)
    tdat = [sqr[0],sqr[2],sqr[3]]
    r3.append(tdat)
    numr3 += 2

#get triangles from: pentagons
for pnt in r5:
    tdat = [pnt[0],pnt[1],pnt[2]]
    r3.append(tdat)
    tdat = [pnt[0],pnt[2],pnt[3]]
    r3.append(tdat)
    tdat = [pnt[0],pnt[3],pnt[4]]
    r3.append(tdat)
    numr3 += 3

#get triangles from: hexagons                     
for hex in r6:
    tdat = [hex[0],hex[1],hex[2]]
    r3.append(tdat)
    tdat = [hex[0],hex[2],hex[3]]
    r3.append(tdat)
    tdat = [hex[0],hex[3],hex[4]]
    r3.append(tdat)
    tdat = [hex[0],hex[4],hex[5]]
    r3.append(tdat)
    numr3 += 4

print "Total triangles: "+str(numr3)

#write vertices to obj file
objfile = open("output.obj",'w')
objfile.write("File created by atom2obj.py\n")
objfile.write("o gmesh\n")
for crd in coords:
    crdline = "v "+str(crd[0]*0.1)+" "+str(crd[1]*0.1)+" "+str(crd[2]*0.1)+"\n"
    objfile.write(crdline)

#write faces to obj file
for tri in r3:
    fline = "f "+str(tri[0]+1)+" "+str(tri[1]+1)+" "+str(tri[2]+1)+"\n"
    objfile.write(fline)

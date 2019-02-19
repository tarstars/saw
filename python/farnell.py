#!/usr/bin/python

from Scientific.Functions import Polynomial
import numpy as np
from math import *
from cmath import *

def rotz(phi):
    return np.array([[cos(phi),-sin(phi),0 ],[sin(phi),cos(phi),0],[0,0,1]])

def roty(phi):
    return np.array([[cos(phi),0,sin(phi) ],[0,1,0],[-sin(phi),0,cos(phi)]])

def rotx(phi):
    return np.array([[1, 0, 0],[0, cos(phi), -sin(phi)],[0, sin(phi), cos(phi)]])

def mul(a,b):
    ret = np.empty((3,3))
    for p in range(3):
        for q in range(3):
            ret[p,q] = 0
            for t in range(3):
                ret[p,q] = ret[p,q] + a[p,t] * b[t,q]
    return ret
            

def mrot(alpha, beta, gamma):
    return mul(mul(rotx(alpha) , roty(beta)), rotz(gamma))

def ten_rot(c,ma):
    ret = np.empty((3,3,3,3))
    for p in range(3):
        for q in range(3):
            for r in range(3):
                for s in range(3):
                    ret[p,q,r,s] = 0
                    for k in range(3):
                        for l in range(3):
                            for m in range(3):
                                for n in range(3):
                                    ret[p,q,r,s] = ret[p,q,r,s] + c[k,l,m,n]*ma[p,k]*ma[q,l]*ma[r,m]*ma[s,n]
    return ret

c11 = 5.6
c12 = 5.1
c13 = 2.2
c33 = 10.6
c44 = 2.65
c66 = 6.6
rho = 6

c = np.empty((3,3,3,3))
for p in range(3):
    for q in range(3):
        for r in range(3):
            for s in range(3):
                c[p,q,r,s] = 0


c[0,0,0,0] = c[1,1,1,1] = c11

c[0,0,1,1] = c[1,1,0,0] = c12

c[0,0,2,2] = c[2,2,0,0] = c[1,1,2,2] = c[2,2,1,1] = c13

c[2,2,2,2] = c33

c[1,2,1,2] = c[1,2,2,1] = c[2,1,1,2] = c[2,1,2,1] = \
c[0,2,0,2] = c[0,2,2,0] = c[2,0,0,2] = c[2,0,2,0] = c44

c[0,1,0,1] = c[0,1,1,0] = c[1,0,0,1] = c[1,0,1,0] = c66


gamm = 0.1

crist = {}
for i in range(3):
    for l in range(3):
        if i!=l:
            m = 0
        else:
            m = gamm
        crist[(i,l)] = (c[2,i,2,l], c[0,i,2,l] + c[2,i,0,l], c[0,i,0,l] - m)

pd = ()
program = (    (0,0,1,1,2,2,1),
               (0,1,1,2,2,0,1),
               (0,2,1,0,2,1,1),
               (0,2,1,1,2,0,-1),
               (0,0,1,2,2,1,-1),
               (0,1,1,0,2,2,-1) )

det = 0
for t in program:
    a = crist[(t[0],t[1])]
    b = crist[(t[2],t[3])]
    m1 = np.polymul(a , b)
    m2 = np.polymul(crist[(t[4],t[5])], t[6])
    det = np.polyadd(det, np.polymul(m1,m2) )

n3s = [t for t in np.roots(det) if t.imag < 0]

for n3 in n3s:
    print "n3 = ", n3
    gc = [[0,0,0],[0,0,0],[0,0,0]]
    for p in range(3):
        for q in range(3):
            gc[p][q] = np.polyval(crist[(p,q)], n3)

    print "crist = "
    print gc
    print
          
    ei = np.linalg.eig(gc)

    v = abs(ei[0][0])
    nv = 0

    if (abs(ei[0][1]) < v):
        v = abs(ei[0][1])
        nv = 1

    if (abs(ei[0][2]) < v):
        v = abs(ei[0][2])
        nv = 2

    print "minimum eigval:"
    print ei[0][nv]
    print
    print "eigvect:"
    print ei[1][:,nv]
    print

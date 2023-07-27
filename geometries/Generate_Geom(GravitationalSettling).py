# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 17:23:31 2020

@author: chris
"""

E = 2e8
v = 0.3
rho = 7e-6
go = -9800

a = 75
Nx = 20
Ny = 20
Nz = 20
h = a/(Nx-1)

#%% Outer Rectangular Boundary

xleft = []
yleft = []
zleft = []
nxleft = []
nyleft = []
nzleft = []
dispXleft = []
dispYleft = [] 
dispZleft = []
forceXleft = []
forceYleft = []
forceZleft = []
bxleft = []
byleft = []
bzleft = []

for i in range(Ny):
    for j in range(Nz):
        xleft.append(0.0)
        yleft.append(i*h)
        zleft.append(j*h)
        nxleft.append(-1.0)
        nyleft.append(0.0)
        nzleft.append(0.0)
        dispXleft.append(0.0)
        dispYleft.append('nan')
        dispZleft.append('nan')
        forceXleft.append('nan')
        forceYleft.append(0.0)
        forceZleft.append(0.0)
        bxleft.append(0.0)
        byleft.append(-rho*go)
        bzleft.append(0.0)
    
    
xright = []
yright = []
zright = []
nxright = []
nyright = []
nzright = []
dispXright = []
dispYright = [] 
dispZright = []
forceXright = []
forceYright = []
forceZright = []
bxright = []
byright = []
bzright = []

for i in range(Ny):
    for j in range(Nz):
        xright.append((Nx-1)*h)
        yright.append(i*h)
        zright.append(j*h)
        nxright.append(1.0)
        nyright.append(0.0)
        nzright.append(0.0)
        dispXright.append('nan')
        dispYright.append('nan')
        dispZright.append('nan')
        forceXright.append(0.0)
        forceYright.append(0.0)
        forceZright.append(0.0)
        bxright.append(0.0)
        byright.append(-rho*go)
        bzright.append(0.0)
     
xbottom = []
ybottom = []
zbottom = []
nxbottom = []
nybottom = []
nzbottom = []
dispXbottom = []
dispYbottom = []
dispZbottom = []
forceXbottom = []
forceYbottom = []
forceZbottom = []
bxbottom = []
bybottom = []
bzbottom = []

for i in range(1,Nx-1):
    for j in range(Nz):
        xbottom.append(i*h)
        ybottom.append(0.0)
        zbottom.append(j*h)
        nxbottom.append(0.0)
        nybottom.append(-1.0)
        nzbottom.append(0.0)
        dispXbottom.append('nan')
        dispYbottom.append(0.0)
        dispZbottom.append('nan')
        forceXbottom.append(0.0)
        forceYbottom.append('nan')
        forceZbottom.append(0.0)
        bxbottom.append(0.0)
        bybottom.append(-rho*go)
        bzbottom.append(0.0)
    
    
xtop = []
ytop = []
ztop = []
nxtop = []
nytop = []
nztop = []
dispXtop = []
dispYtop = []
dispZtop = []
forceXtop = []
forceYtop = []
forceZtop = []
bxtop = []
bytop = []
bztop = []

for i in range(1,Nx-1):
    for j in range(Nz):
        xtop.append(i*h)
        ytop.append((Ny-1)*h)
        ztop.append(j*h)
        nxtop.append(0.0)
        nytop.append(1.0)
        nztop.append(0.0)
        dispXtop.append('nan')
        dispYtop.append('nan')
        dispZtop.append('nan')
        forceXtop.append(0.0)
        forceYtop.append(0.0)
        forceZtop.append(0.0)
        bxtop.append(0.0)
        bytop.append(-rho*go)
        bztop.append(0.0)
        

xfront= []
yfront = []
zfront = []
nxfront = []
nyfront = []
nzfront = []
dispXfront = []
dispYfront = []
dispZfront = []
forceXfront = []
forceYfront = []
forceZfront = []
bxfront = []
byfront = []
bzfront = []

for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        xfront.append(i*h)
        yfront.append(j*h)
        zfront.append(0)
        nxfront.append(0.0)
        nyfront.append(0.0)
        nzfront.append(-1.0)
        dispXfront.append('nan')
        dispYfront.append('nan')
        dispZfront.append(0.0)
        forceXfront.append(0.0)
        forceYfront.append(0.0)
        forceZfront.append('nan')
        bxfront.append(0.0)
        byfront.append(-rho*go)
        bzfront.append(0.0)
        
        
xback= []
yback = []
zback = []
nxback = []
nyback = []
nzback = []
dispXback = []
dispYback = []
dispZback = []
forceXback = []
forceYback = []
forceZback = []
bxback = []
byback = []
bzback = []

for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        xback.append(i*h)
        yback.append(j*h)
        zback.append((Nz-1)*h)
        nxback.append(0.0)
        nyback.append(0.0)
        nzback.append(1.0)
        dispXback.append('nan')
        dispYback.append('nan')
        dispZback.append('nan')
        forceXback.append(0.0)
        forceYback.append(0.0)
        forceZback.append(0.0)
        bxback.append(0.0)
        byback.append(-rho*go)
        bzback.append(0.0)


#%% Inner nodes
    
xinner = []
yinner = []
zinner = []
nxinner = []
nyinner = []
nzinner = []
dispXinner = []
dispYinner = []
dispZinner = []
forceXinner = []
forceYinner = []
forceZinner = []
bxinner = []
byinner = []
bzinner = []

for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        for k in range(1,Nz-1):
            xtemp = i*h
            ytemp = j*h
            ztemp = k*h
            xinner.append(xtemp)
            yinner.append(ytemp)
            zinner.append(ztemp)
            nxinner.append('nan')
            nyinner.append('nan')
            nzinner.append('nan')
            dispXinner.append('nan')
            dispYinner.append('nan')
            dispZinner.append('nan')
            forceXinner.append('nan')
            forceYinner.append('nan')
            forceZinner.append('nan')
            bxinner.append(0.0)
            byinner.append(-rho*go)
            bzinner.append(0.0)
        
    
#%% Combining the boundary and inner nodes list
    
x = xtop + xbottom + xright + xleft + xfront + xback + xinner
y = ytop + ybottom + yright + yleft + yfront + yback + yinner
z = ztop + zbottom + zright + zleft + zfront + zback + zinner
nx = nxtop + nxbottom + nxright + nxleft + nxfront + nxback + nxinner
ny = nytop + nybottom + nyright + nyleft + nyfront + nyback + nyinner
nz = nztop + nzbottom + nzright + nzleft + nzfront + nzback + nzinner
dispX = dispXtop + dispXbottom + dispXright + dispXleft + dispXfront + dispXback + dispXinner
dispY = dispYtop + dispYbottom + dispYright + dispYleft + dispYfront + dispYback + dispYinner
dispZ = dispZtop + dispZbottom + dispZright + dispZleft + dispZfront + dispZback + dispZinner
forceX = forceXtop + forceXbottom + forceXright + forceXleft + forceXfront + forceXback + forceXinner
forceY = forceYtop + forceYbottom + forceYright + forceYleft + forceYfront + forceYback + forceYinner
forceZ = forceZtop + forceZbottom + forceZright + forceZleft + forceZfront + forceZback + forceZinner
bx = bxtop + bxbottom + bxright + bxleft + bxfront + bxback + bxinner
by = bytop + bybottom + byright + byleft + byfront + byback + byinner
bz = bztop + bzbottom + bzright + bzleft + bzfront + bzback + bzinner

#print(f"{np.sum(probe)}")

#%% Writing the txt data

file = open("geom_param.txt","w")

file.write(f"{E}")
file.write("\n")
file.write(f"{v}")
file.write("\n")
file.write(f"{h}")

file.close()

file = open("geom_data.txt","w")

for i in range(len(x)):
    file.write(f"{x[i]}" + "\t" + f"{y[i]}" + "\t" + f"{z[i]}" "\t" + f"{nx[i]}" + "\t" + f"{ny[i]}" + "\t" + f"{nz[i]}"+ "\t" + f"{dispX[i]}" + "\t" + f"{dispY[i]}" + "\t" + f"{dispZ[i]}" + "\t" + f"{forceX[i]}" + "\t" + f"{forceY[i]}" + "\t" + f"{forceZ[i]}" + "\t" + f"{bx[i]}" + "\t" + f"{by[i]}" + "\t" + f"{bz[i]}")
    file.write("\n")
    
file.close()


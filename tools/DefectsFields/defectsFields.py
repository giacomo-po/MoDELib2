# conda deactivate
# cd build
# cmake ..
# make -j6
# cd ..
# /usr/local/bin/python3.10 defectsFields.py
import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("./build")
import DefectsFields

# Define the simulation directory. This includes subdirectories inputFiles and evl
simulationDir="../../tutorials/DislocationDynamics/periodicDomains/dipoleNoise/"

# create the DefectsFieldsExtractor object
dfe=DefectsFields.DefectsFieldsExtractor(simulationDir)

# Read a pre-existing configuration file (e.g. readConfiguration(X) reads file simulationDir/evl/evl_X.txt)
#ad.readConfiguration(567)

# Alternatively, generate a new configuration using simulationDir/inputFiles/initialMicrostructure.txt
dfe.readMicrostructure()
#ad.writeConfiguration(0) # Optional. Write the generated configuration to file (writeConfiguration(X) writes file simulationDir/evl/evl_X.txt)

# grab the domain corners
ldc=dfe.lowerDomainCorner()
udc=dfe.upperDomainCorner()

# Extracting grids of values on a a y-z plane: e.g. solid angle and sigma_11
n=100
y=np.linspace(ldc[1], udc[1], num=n) # grid x-range
z=np.linspace(ldc[2], udc[2], num=n) # grid z-range
x=0.5*(ldc[0]+udc[0]) # grid position in y
solidAngle=np.empty([n, n]) # grid of solid angle values
s11=np.empty([n, n]) # grid of sigma_11 values
for i in range(0,y.size):
    for j in range(0,z.size):
        solidAngle[j,i]=dfe.solidAngle(x,y[i],z[j])
        stress=dfe.dislocationStress(x,y[i],z[j])
        s11[j,i]=stress[0,0]

fig=plt.figure()
plt.imshow(solidAngle, origin='lower',cmap='jet')
plt.colorbar()
plt.show()

fig=plt.figure()
plt.imshow(s11, origin='lower',cmap='jet')
plt.colorbar()
plt.show()


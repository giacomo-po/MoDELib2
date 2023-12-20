import sys
sys.path.append("../../../../../python/")
from modlibUtils import *


# Make a local copy of simulation parameters file and modify that copy if necessary
DDfile='DD.txt'
shutil.copy2('../../'+DDfile, '.')
#setInputVariable(DDfile,'Lmin','10')  # min segment length (in Burgers vector units)
#setInputVariable(DDfile,'Lmax','50')  # max segment length (in Burgers vector units)
#setInputVariable(DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces


# set up overall simulation size
m=2 # there are (2m)^3 precipitates in the box
materialFile='../../../MaterialsLibrary/Ni.txt'
b_SI=getValueInFile(materialFile,'b_SI')
targetCubeSide=250e-9/b_SI # Cazares, 2022, Acta Mater 241
f=0.447 # Cazares, 2022, Acta Mater 241
targetCubeSpacing=targetCubeSide/np.cbrt(f)
nBox=np.round(2*m*targetCubeSpacing/np.sqrt(2.)) # simiulation box length in number of Burgers vectors
cubeSpacing=nBox*np.sqrt(2.)/2./m
cubeSide=cubeSpacing*np.cbrt(f)
print('cubes side='+str(cubeSide*b_SI))
print('cubes side='+str(cubeSpacing*b_SI))
print('volume fraction='+str(np.power(cubeSide/cubeSpacing,3)))

pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile='../../../MeshLibrary/unitCube.msh'
#pf.grain1globalX1=np.array([0,1,-1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
#pf.grain1globalX3=np.array([1,1,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
#pf.boxEdges=np.array([[0,1,-1],[-1,1,0],[1,1,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([nBox,nBox,nBox]) # must be a vector of integers
pf.X0=np.array([0.5,0.5,0.5]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([0,1,2,3,4,5])
pf.write()

# Edit polyhedron file
polyhedronFile='polyhedronInclusionIndividual.txt'
#C*s*(x+0.5)+C*L*[i j k]=C*s*(x+0.5+L/s*[i j k])
X0=np.empty((0,3))
for i in range(-m,m):
    for j in range(-m,m):
        for k in range(-m,m):
            X0=np.append(X0,cubeSpacing/cubeSide*(np.array([[i,j,k]])+0.5)+0.5,0)
setInputMatrix(polyhedronFile,'x0',X0)
setInputMatrix(polyhedronFile,'A',pf.C2G*cubeSide)

# Edit planarLoop file
planarLoopFile1='planarLoop1.txt'
L1=200
L2=(cubeSpacing-cubeSide)*np.sqrt(2.0)
n1=np.array([1,1,1])@pf.C2G.transpose() # primary plane normal
n1=n1/np.linalg.norm(n1)
#b1=np.array([1,0,-1])@pf.C2G.transpose() # primary plane Burgers
b1=np.array([0,1,-1])@pf.C2G.transpose() # primary plane Burgers
b1=b1/np.linalg.norm(b1)
d1=np.array([0,1,-1])@pf.C2G.transpose() # direction along channel in primary plane
d2=np.cross(n1,d1) # direction across channel in primary plane
d1=d1/np.linalg.norm(d1)*L1
d2=d2/np.linalg.norm(d2)*L2
loopNodes=np.array([-d1/2-d2/2,d1/2-d2/2,d1/2+d2/2,-d1/2+d2/2])
setInputMatrix(planarLoopFile1,'loopNodes',loopNodes)
setInputVector(planarLoopFile1,'loopNormal',n1,'loop normal')
setInputVector(planarLoopFile1,'loopBurgers',b1,'loop Burgers')

# Edit stress
uniformExternalLoadControllerFile='uniformExternalLoadController.txt'
loadAxis=np.array([0,0,1])@pf.C2G.transpose()
loadAxis=loadAxis/np.linalg.norm(loadAxis)
stress=0.01
ExternalStress0=stress*np.outer(loadAxis,loadAxis)
setInputMatrix(uniformExternalLoadControllerFile,'ExternalStress0',ExternalStress0)

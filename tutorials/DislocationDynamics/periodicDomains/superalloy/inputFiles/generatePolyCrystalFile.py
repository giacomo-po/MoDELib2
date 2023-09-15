import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

pf=PolyCrystalFile('../../../MaterialsLibrary/Ni.txt');
pf.absoluteTemperature=300;
pf.enablePartials=0;
pf.dislocationMobilityType='default'
pf.meshFile='../../../MeshLibrary/unitCube.msh'
#pf.grain1globalX1=np.array([0,1,-1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
#pf.grain1globalX3=np.array([1,1,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
#pf.alignToSlipSystem0=0
pf.boxEdges=np.array([[0,1,-1],[-1,1,0],[1,1,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([10000,10000,4000]) # must be a vector of integers
pf.X0=np.array([0.5,0.5,0.5]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([0,1,2,3,4,5])
pf.write()

# Edit polyhedron file
polyhedronFile='polyhedronInclusionIndividual.txt'
#C*s*(x+0.5)+C*L*[i j k]=C*s*(x+0.5+L/s*[i j k])
# f=Vcube/Vbox --> f^1/3 = Lcube/Lbox  --> cubeSpacing=cubeSide/f^1/3
f=0.447 # Cazares, 2022, Acta Mater 241
cubeSide=1000 # Cazares, 2022, Acta Mater 241
cubeSpacing=cubeSide/np.cbrt(f)
X0=np.empty((0,3))
for i in range(-1,1):
    for j in range(-1,1):
        for k in range(-1,1):
            X0=np.append(X0,cubeSpacing/cubeSide*(np.array([[i,j,k]])+0.5)+0.5,0)
setInputMatrix(polyhedronFile,'x0',X0)
setInputMatrix(polyhedronFile,'A',pf.C2G*cubeSide)

# Edit planarLoop file
planarLoopFile1='planarLoop1.txt'
L1=2000
L2=(cubeSpacing-cubeSide)*np.sqrt(2.0)
n1=np.array([1,1,1])@pf.C2G.transpose() # primary plane normal
n1=n1/np.linalg.norm(n1)
#b1=np.array([1,0,-1])@pf.C2G.transpose() # primary plane normal
b1=np.array([0,1,-1])@pf.C2G.transpose() # primary plane normal
b1=b1/np.linalg.norm(b1)
d1=np.array([0,1,-1])@pf.C2G.transpose() # direction along channel in primary plane
d2=np.cross(n1,d1) # direction across channel in primary plane
d1=d1/np.linalg.norm(d1)*L1
d2=d2/np.linalg.norm(d2)*L2
loopNodes=np.array([-d1/2-d2/2,d1/2-d2/2,d1/2+d2/2,-d1/2+d2/2])
setInputMatrix(planarLoopFile1,'loopNodes',loopNodes)
setInputMatrix(planarLoopFile1,'loopNormal',np.array([n1]))
setInputMatrix(planarLoopFile1,'loopBurgers',np.array([b1]))

# Edit stress
uniformExternalLoadControllerFile='uniformExternalLoadController.txt'
loadAxis=np.array([0,0,1])@pf.C2G.transpose()
loadAxis=loadAxis/np.linalg.norm(loadAxis)
stress=0.01
ExternalStress0=stress*np.outer(loadAxis,loadAxis)
setInputMatrix(uniformExternalLoadControllerFile,'ExternalStress0',ExternalStress0)

import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

pf=PolyCrystalFile('../../../MaterialsLibrary/Cu.txt');
pf.absoluteTemperature=300;
pf.enablePartials=1;
pf.meshFile='../../../MeshLibrary/unitCube.msh'
pf.grain1globalX1=np.array([0,1,1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
pf.grain1globalX3=np.array([-1,1,-1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[0,1,1],[2,1,-1],[-1,1,-1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([200,200,200]) # must be a vector of integers
pf.X0=np.array([0.5,0.5,0.5]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([0,1,2,3,4,5])
pf.write()

# Edit periodicDipole file
periodicDipoleFile='periodicDipole.txt'
setInputVector(periodicDipoleFile,'periodicDipoleSlipSystemIDs',np.array([0,1]),'slip system IDs for each dipole')
setInputVector(periodicDipoleFile,'periodicDipoleExitFaceIDs',np.array([0,0]),'1 is for edge, 0 for screw')
setInputMatrix(periodicDipoleFile,'periodicDipolePoints',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]))
setInputVector(periodicDipoleFile,'periodicDipoleHeights',np.array([200,200]),'height of each dipole, in number of planes')
setInputVector(periodicDipoleFile,'periodicDipoleGlideSteps',np.array([10.0,30.0]),'step of each dipole in the glide plane')

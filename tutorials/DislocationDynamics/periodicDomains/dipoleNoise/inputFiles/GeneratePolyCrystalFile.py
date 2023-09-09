import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

pf=PolyCrystalFile('../../../MaterialsLibrary/AlMg5.txt');
pf.absoluteTemperature=300;
pf.enablePartials=0;
pf.dislocationMobilityType='default'
pf.meshFile='../../../MeshLibrary/unitCube.msh'
pf.alignToSlipSystem0=1
pf.boxScaling=np.array([200,200,200]) # must be a vector of integers
pf.X0=np.array([0.5,0.5,0.5]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([0,1,2,3,4,5])

pf.write()

#print(pf.A)

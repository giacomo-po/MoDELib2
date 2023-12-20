import sys
sys.path.append("../../../../../python/")
from modlibUtils import *

# Make a local copy of simulation parameters file and modify that copy if necessary
DDfile='DD.txt'
shutil.copy2('../../'+DDfile, '.') 

pf=PolyCrystalFile('../../../MaterialsLibrary/Cu.txt');
pf.absoluteTemperature=300;
pf.dislocationMobilityType='default'
pf.meshFile='../../../MeshLibrary/unitCube.msh'
#pf.grain1globalX1=np.array([1,2,1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
#pf.grain1globalX3=np.array([1,1,-3])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.alignToSlipSystem0=0
#pf.boxEdges=np.array([[1,2,1],[2,1,1],[1,1,-3]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([1000,1000,1000]) # must be a vector of integers
pf.X0=np.array([0.5,0.5,0.5]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([0,1,2,3,4,5])

pf.write()

#print(pf.A)

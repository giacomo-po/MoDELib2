import sys, string, os, fileinput
from fractions import Fraction
import numpy as np
# evl file

class EVL:
    nodes=np.empty([0,0])
        
class PolyCrystalFile(dict):
    materialFile=''
    crystalStructure=''
    absoluteTemperature=300.0
    enablePartials=0
    dislocationMobilityType='default'
    meshFile='../../../MeshLibrary/unitCube.msh'
    grain1globalX1=np.array([1,0,0]) # overwritten if alignToSlipSystem0=true
    grain1globalX3=np.array([0,0,1]) # overwritten if alignToSlipSystem0=true
    A=np.zeros((3,3))
    invA=np.zeros((3,3))
    alignToSlipSystem0=0
    boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge. Overwritten if alignToSlipSystem0=true
    boxScaling=np.array([1000,1000,1000]) # must be a vector of integers
    C2G=np.zeros((3,3))
    F=np.zeros((3,3))
    X0=np.array([0,0,0])
    periodicFaceIDs=np.array([0,1,2,3,4,5])
    
    def __init__(self, materialFile):
        self.materialFile = materialFile
        self.crystalStructure = getStringInFile(materialFile,'crystalStructure')
        if self.crystalStructure == 'FCC':
            self.A=np.array([[0.,1.,1.],[1.,0.,1.],[1.,1.,0.]])/np.sqrt(2.0)
        elif self.crystalStructure == 'BCC':
            self.A=np.array([[-1.,1.,1.],[1.,-1.,1.],[1.,1.,-1.]])/np.sqrt(3.0)
        else:
            raise Exception("Unkonwn crystalStructure "+self.crystalStructure)
        self.invA=np.linalg.inv(self.A)
        np.set_printoptions(precision=15)
        
    def write(self):
        if self.alignToSlipSystem0:
            if self.crystalStructure == 'BCC':
                self.grain1globalX1=np.array([1,1,-1]) # overwrite
                self.grain1globalX3=np.array([1,0,1])  # overwrite
            elif self.crystalStructure == 'FCC':
                self.grain1globalX1=np.array([0,1,1]) # overwrite
                self.grain1globalX3=np.array([-1,1,-1])  # overwrite
            else:
                raise Exception("Unkonwn crystalStructure "+self.crystalStructure)
        
        x1=self.grain1globalX1/np.linalg.norm(self.grain1globalX1);
        x3=self.grain1globalX3/np.linalg.norm(self.grain1globalX3);
        self.C2G=np.array([x1,np.cross(x3,x1),x3]);
        
        if self.alignToSlipSystem0:
            self.boxEdges=self.C2G
        
        # Find lattice vectors (columns of L) aligned to boxEdges
        L=np.zeros((3,3))
        B=self.invA@self.boxEdges.transpose()
        for j in range(0, 3):
            b=B[:,j]/np.max(np.abs(B[:,j]))
            n=np.array([0,0,0], dtype=int)
            d=np.array([1,1,1], dtype=int)
            for i in range(0, 3):
                f=Fraction(b[i]).limit_denominator(100)
                n[i]=f.numerator
                d[i]=f.denominator
            dp=np.prod(d);
            nr=np.array([1,1,1], dtype=int)
            for i in range(0, 3):
                nr[i]=n[i]*dp/d[i]
            nr=nr/np.gcd.reduce(nr)
            L[:,j]=self.A@nr.transpose()
            self.F[:,j]=self.C2G@L[:,j]*self.boxScaling[j]
        #print(self.F)

        polyFile = open("polycrystal.txt", "w")
        polyFile.write('materialFile='+self.materialFile+';\n')
        polyFile.write('absoluteTemperature='+str(self.absoluteTemperature)+'; # [K] simulation temperature \n')
        polyFile.write('enablePartials='+str(self.enablePartials)+'; # enables partial slip systems \n')
        polyFile.write('dislocationMobilityType='+self.dislocationMobilityType+'; # default or FCC,BCC,HEXbasal,HEXprismatic,HEXpyramidal \n')
        polyFile.write('meshFile='+self.meshFile+'; # mesh file \n')
        polyFile.write('C2G1='+' '.join(map(str, self.C2G[0,:]))+'\n')
        polyFile.write('     '+' '.join(map(str, self.C2G[1,:]))+'\n')
        polyFile.write('     '+' '.join(map(str, self.C2G[2,:]))+'; # crystal rotation matrix \n')
        polyFile.write('F='+' '.join(map(str, self.F[0,:]))+'\n')
        polyFile.write('  '+' '.join(map(str, self.F[1,:]))+'\n')
        polyFile.write('  '+' '.join(map(str, self.F[2,:]))+'; # mesh deformation gradient. Mesh nodes X are mapped to x=F*(X-X0) \n')
        polyFile.write('X0='+' '.join(map(str, self.X0))+'; # mesh shift. Mesh nodes X are mapped to x=F*(X-X0) \n')
        polyFile.write('periodicFaceIDs='+' '.join(map(str, self.periodicFaceIDs))+'; # IDs of faces labelled as periodic \n')

        polyFile.write('gridSize=256 256; # size of grid on the glide plane\n');
        polyFile.write('gridSpacing_SI=1e-10 1e-10; # [m] spacing of grid on the glide plane\n');
        polyFile.write('solidSolutionNoiseMode=0; # 0=no noise, 1= read noise, 2=compute noise\n')
        polyFile.write('solidSolutionNoiseFile_xz=../../../NoiseLibrary/noise_xz.vtk;\n');
        polyFile.write('solidSolutionNoiseFile_yz=../../../NoiseLibrary/noise_yz.vtk;\n');
        polyFile.write('stackingFaultNoiseMode=0; # 0=no noise\n');
        polyFile.write('spreadLstress_A=1; # add comment\n');
        polyFile.write('a_cai_A=1; # add comment\n');
        polyFile.write('seed=0; # add comment\n');

        polyFile.close()

def readEVLtxt(filename):
    evlFile = open(filename+'.txt', "r")
    numNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    evl=EVL();
    evl.nodes=np.empty([numNodes, 3])
    for k in np.arange(numNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.nodes[k,:]=data[1:4]
    return evl

def getStringInFile(fileName,variable):
    with open(fileName) as f:
        datafile = f.readlines()
    found = False  # This isn't really necessary
    for line in datafile:
#        print(line)
        if variable in line:
            # found = True # Not necessary
            foundEqual=line.find('=');
            foundSemiCol=line.find(';');
            if line[0:foundEqual].strip()==variable:
                return line[foundEqual+1:foundSemiCol].strip()
    return 'Not Found'  # Because you finished the search without finding

def getValueInFile(fileName,variable):
    return float(getStringInFile(fileName,variable))

def readFfile(folder):
    F=np.loadtxt(folder +'/F_0.txt');
    with open('./F/F_labels.txt') as f:
        lines = f.readlines()
        for idx in range(len(lines)):
            lines[idx] = lines[idx].rstrip()
    return F,lines
    
def getFarray(F,Flabels,label):
    k=0;
    for line in Flabels:
        if line==label:
            return F[:,k]
        k=k+1
    return np.zeros(shape=(0,0))

def setInputVariable(fileName,variable,newVal):
    with fileinput.FileInput(fileName, inplace=True) as file:
        for line in file:
            if variable in line:
                foundPound=line.find('#');
                foundEqual=line.find('=');
                foundSemiCol=line.find(';');
#                print(foundPound, file=sys.stderr)
#                print(foundEqual, file=sys.stderr)
#                print(foundSemiCol, file=sys.stderr)
                if (foundPound==-1 or foundPound > foundSemiCol) and (foundSemiCol>foundEqual):
                    oldVal=line[foundEqual+1:foundSemiCol]
                    line = line.replace(oldVal,newVal)
                    print(line, file=sys.stderr)
            sys.stdout.write(line)
            
def setInputMatrix(fileName,variable,newVal):
    with fileinput.FileInput(fileName, inplace=True) as file:
        num_rows, num_cols = newVal.shape
        writeLine=True
        for line in file:
            foundPound=line.find('#');
            foundEqual=line.find('=');
            foundSemiCol=line.find(';');
            if (foundPound==-1 or foundPound > foundSemiCol) and foundEqual>0:
                writeLine=True
            if variable in line:
                if foundPound==-1 and foundEqual>=len(variable):
                    line=variable+'='
                    for i in range(0,num_rows-1):
                        line=line+' '.join(map(str, newVal[i,:]))+'\n'
                    line=line+' '.join(map(str, newVal[num_rows-1,:]))+';\n'
                    sys.stdout.write(line)
                    writeLine=False # need to skip lines until we find another variable declaration to erase old matrix
            if writeLine:
                sys.stdout.write(line)


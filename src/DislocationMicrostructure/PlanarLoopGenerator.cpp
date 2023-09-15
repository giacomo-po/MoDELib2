/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarLoopGenerator_cpp_
#define model_PlanarLoopGenerator_cpp_

#include <numbers>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PlanarLoopGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    PlanarLoopGenerator::PlanarLoopGenerator(const std::string& fileName) :
    /* init */ MicrostructureGeneratorBase(fileName)
    {
        
    }

    void PlanarLoopGenerator::generateSingle(MicrostructureGenerator& mg,const VectorDimD& b,const VectorDimD& n,const std::vector<VectorDimD>& loopPoints)
    {
        if(loopPoints.size()>=3)
        {
            // Grain is defined by first point
            std::pair<bool,const Simplex<dim,dim>*> found(mg.ddBase.mesh.search(loopPoints[0]));
            if(!found.first)
            {
                std::cout<<"Point "<<loopPoints[0].transpose()<<" is outside mesh. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            const int grainID(found.second->region->regionID);
            //assert(mg.ddBase.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
            const auto& grain(mg.ddBase.poly.grain(grainID));
            
            const ReciprocalLatticeDirection<dim> r(grain.singleCrystal->reciprocalLatticeDirection(n));
            GlidePlaneKey<3> glidePlaneKey(loopPoints[0], r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
//            const GlidePlane<dim> glidePlane(loopPoints[0],r);
            for(const auto& point : loopPoints)
            {
                if(!glidePlane->referencePlane->contains(point))
                {
                    throw std::runtime_error("Plane does not contain point.");
                }
            }
            

                mg.insertJunctionLoop(loopPoints,glidePlane,
                                      b,n,
                                      loopPoints[0],grainID,std::fabs(b.normalized().dot(n.normalized()))<FLT_EPSILON? DislocationLoopIO<dim>::GLISSILELOOP : DislocationLoopIO<dim>::SESSILELOOP);

        }
        

    }

    void PlanarLoopGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        throw std::runtime_error("PlanarLoopGenerator::generateDensity not implemented");
    }

    void PlanarLoopGenerator::generateIndividual(MicrostructureGenerator& mg)
    {
        
            std::cout<<magentaBoldColor<<"Generating individual planar loop"<<defaultColor<<std::endl;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> loopPointsMatrix(this->parser.readMatrixCols<double>("loopNodes",dim,true));
        std::vector<VectorDimD> loopPoints;
        for(int i=0;i<loopPointsMatrix.rows();++i)
        {
            loopPoints.push_back(loopPointsMatrix.row(i));
        }
        generateSingle(mg,
                       this->parser.readMatrix<double>("loopBurgers",1,dim,true),
                       this->parser.readMatrix<double>("loopNormal",1,dim,true),
                       loopPoints);
        
    }

}
#endif

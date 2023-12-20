/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * Additional contributors:
 *  Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationDynamicsBase_CPP_
#define model_DislocationDynamicsBase_CPP_

#include <DislocationDynamicsBase.h>
#include <CTM.h>
namespace model
{

template <int dim>
double getEwaldLength(const std::vector<Eigen::Matrix<double,dim,1>>& periodicBasis,const double& EwaldLengthFactor)
{
    
    if(periodicBasis.size())
    {
        // Compute generalized volume of periodic lattice cell
        Eigen::MatrixXd B(Eigen::MatrixXd::Zero(3,periodicBasis.size()));
        for(size_t k=0;k<periodicBasis.size();++k)
        {
            B.col(k)=periodicBasis[k];
        }
        const double vol(sqrt((B.transpose()*B).determinant())/CTM::factorial(periodicBasis.size()));
//        std::cout<<"vol="<<vol<<std::endl;
//        std::cout<<"edge="<<std::pow(vol,1.0/periodicBasis.size())<<std::endl;
//        std::cout<<"elength="<<EwaldLengthFactor*std::pow(vol,1.0/periodicBasis.size())<<std::endl;
        return EwaldLengthFactor*std::pow(vol,1.0/periodicBasis.size());
    }
    else
    {
        return 0.0;
    }
}

template <int dim>
Eigen::Matrix<double,dim,Eigen::Dynamic> getPeriodicLatticeBasis(const std::vector<Eigen::Matrix<double,dim,1>>& basisVector)
{
    Eigen::Matrix<double,dim,Eigen::Dynamic> temp(Eigen::MatrixXd::Zero(dim,basisVector.size()));
    for(size_t k=0;k<basisVector.size();++k)
    {
        temp.col(k)=basisVector[k];
    }
    return temp;
}

template <int _dim>
DislocationDynamicsBase<_dim>::DislocationDynamicsBase( const std::string& folderName) :
/* init */ simulationParameters( folderName)
/* init */,mesh(simulationParameters.traitsIO.meshFile,
                TextFileParser(
                               simulationParameters.traitsIO.polyFile
                               ).readMatrix<double>("F",_dim,_dim,true),
                TextFileParser(
                               simulationParameters.traitsIO.polyFile
                               ).readMatrix<double>("X0",1,_dim,true).transpose(),
                simulationParameters.periodicFaceIDs
                )
/* init */,poly(simulationParameters.traitsIO.polyFile,mesh)
/* init */,glidePlaneFactory(poly)
/* init */,periodicGlidePlaneFactory(poly,glidePlaneFactory)
/* init */,periodicBasis(mesh.periodicBasis())
/* init */,periodicLatticeBasis(getPeriodicLatticeBasis(periodicBasis))
/* init */,periodicLatticeReciprocalBasis(periodicLatticeBasis*(periodicLatticeBasis.transpose()*periodicLatticeBasis).inverse())
/* init */,periodicShifts(mesh.periodicShifts(simulationParameters.periodicImageSize))
/* init */,EwaldLength(getEwaldLength(periodicBasis,TextFileParser(simulationParameters.traitsIO.ddFile).readScalar<double>("EwaldLengthFactor",true)))
{
    if(!mesh.simplices().size())
    {
        throw std::runtime_error("Mesh is empty");
    }
    
    std::cout<<"EwaldLength="<<EwaldLength<<std::endl;
}


template <int _dim>
Eigen::VectorXd DislocationDynamicsBase<_dim>::periodicCoordinates(const VectorDim& x) const
{
    const VectorDim dx(x-mesh.xCenter());
    return periodicLatticeReciprocalBasis.transpose()*dx;
}


template struct DislocationDynamicsBase<3>;

} // namespace model
#endif

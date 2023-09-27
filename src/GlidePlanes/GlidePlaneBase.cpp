/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneBase_cpp_
#define model_GlidePlaneBase_cpp_

#include <map>
#include <LatticeModule.h>
#include <GlidePlaneBase.h>

namespace model
{
    
    GlidePlaneBase::GlidePlaneBase(const LatticeVectorType& v1_in,
                                   const LatticeVectorType& v2_in,
                                   const std::shared_ptr<GammaSurface>& gammaSurface_in) :
    /* init */ ReciprocalLatticeDirectionType(v1_in,v2_in),
    /* init */ primitiveVectors(std::make_pair(v1_in,v2_in))
    /* init */,G2L(getG2L(v1_in.cartesian(),this->cartesian()))
    /* init */,gammaSurface(gammaSurface_in)
    {
    }
        
    typename GlidePlaneBase::LatticeVectorType GlidePlaneBase::snapToLattice(const VectorDimD& P) const
    {
        const Eigen::Matrix<double, 3, 2> B((Eigen::Matrix<double, 2, 3>() << primitiveVectors.first.cartesian().transpose(), primitiveVectors.second.cartesian().transpose()).finished().transpose());
        const Eigen::Matrix<double, 2, 1> nd = (B.transpose() * B).inverse() * B.transpose() * P;
        //            const Eigen::Matrix<double,3,1>  p=B*RoundEigen<double,2>::round(nd);
        const Eigen::Matrix<double, 3, 1> p = B * nd.array().round().matrix();
        //            return LatticeVectorType(p,primitiveVectors.first.covBasis,primitiveVectors.second.contraBasis);
        return LatticeVectorType(p, primitiveVectors.first.lattice);
    }

    typename GlidePlaneBase::VectorDimD GlidePlaneBase::snapToPlane(const VectorDimD& P) const
    {
        //            return LatticeVectorType(p,primitiveVectors.first.covBasis,primitiveVectors.second.contraBasis);
        VectorDimD cn(this->cartesian());
        const double cnNorm(cn.norm());
        assert(cnNorm > FLT_EPSILON);
        cn /= cnNorm;
        return P - P.dot(cn) * cn;
    }

    typename GlidePlaneBase::MatrixDim GlidePlaneBase::getG2L(const VectorDimD& x,const VectorDimD& z)
    {
    //    const double xNorm(x.norm());
    //    const double zNorm(z.norm());
    //    if(xNorm<FLT_EPSILON || zNorm<FLT_EPSILON)
    //    {
    //        throw std::runtime_error("GlidePlaneBase: cannot cumpute G2L from zero vectors.");
    //    }
    //    assert(xNorm>FLT_EPSILON);
    //    assert(zNorm>FLT_EPSILON);
    //    assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
    //    Eigen::Matrix3d temp(Eigen::Matrix3d::Identity());
    //    temp.col(2)=z/zNorm;
    //    temp.col(0)=x/xNorm;
    //    temp.col(1)=temp.col(2).cross(temp.col(0));

        Eigen::Matrix3d temp(Eigen::Matrix3d::Identity());
        temp.row(2)=z.normalized();
        temp.row(0)=x.normalized();
        temp.row(1)=temp.row(2).cross(temp.row(0));
        if((temp*temp.transpose()-Eigen::Matrix3d::Identity()).norm()>FLT_EPSILON)
        {
            std::cout<<"x1="<<temp.row(0)<<std::endl;
            std::cout<<"x2="<<temp.row(1)<<std::endl;
            std::cout<<"x3="<<temp.row(2)<<std::endl;
            std::cout<<"temp="<<temp<<std::endl;
            std::cout<<"temp*temp^T="<<temp*temp.transpose()<<std::endl;
            throw std::runtime_error("GlidePlaneBase: G2L is not orthogonal.");
        }
        if(std::fabs(temp.determinant()-1.0)>FLT_EPSILON)
        {
            throw std::runtime_error("GlidePlaneBase: G2L is not proper.");
        }
        return temp;
        
    //    return temp.transpose();
    }

double GlidePlaneBase::misfitEnergy(const VectorDimD& b) const
{
    if(gammaSurface)
    {
        const VectorDimD bL(G2L*b);
        if(std::fabs(bL(2))>FLT_EPSILON)
        {
            throw std::runtime_error("SLIP VECTOR NOT ON GlidePlaneBase");
        }
        return gammaSurface->operator()(bL.segment<2>(0));
    }
    else
    {
        return 0.0;
    }
}

} // end namespace
#endif

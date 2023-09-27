/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneBase_h_
#define model_GlidePlaneBase_h_

#include <map>
#include <LatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>
#include <GammaSurface.h>

namespace model
{
    struct GlidePlaneBase : public ReciprocalLatticeDirection<3>
    {
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<long int,3,1> VectorDimI;
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3>    LatticeDirectionType;
        typedef ReciprocalLatticeDirection<3> ReciprocalLatticeDirectionType;
        
        const std::pair<LatticeDirectionType,LatticeDirectionType> primitiveVectors;
        const Eigen::Matrix3d G2L;
        const std::shared_ptr<GammaSurface> gammaSurface;

        GlidePlaneBase(const LatticeVectorType& v1_in,
                       const LatticeVectorType& v2_in,
                       const std::shared_ptr<GammaSurface>& gammaSurface_in);


        LatticeVectorType snapToLattice(const VectorDimD& P) const;
        VectorDimD snapToPlane(const VectorDimD& P) const;
        double misfitEnergy(const VectorDimD& b) const;

        static MatrixDim getG2L(const VectorDimD& x,const VectorDimD& z);

        
    };
    
} // end namespace
#endif

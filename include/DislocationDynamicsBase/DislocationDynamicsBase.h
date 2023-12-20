/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * Additional contributors:
 *   Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationDynamicsBase_H_
#define model_DislocationDynamicsBase_H_

#include <string>
#include <Eigen/Dense>

#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <GlidePlaneFactory.h>
#include <PeriodicGlidePlaneFactory.h>


namespace model
{
    template <int _dim>
    struct DislocationDynamicsBase
    {
        // desired members:
        //  those common to MicrostructureGenerator and DefectiveCrystal
        //  so that their members can be references to those in this class
        static constexpr int dim=_dim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef DislocationDynamicsBase<dim> DislocationDynamicsBaseType;

        
        DislocationDynamicsBase( const std::string& folderName);
        DefectiveCrystalParameters simulationParameters;
        const SimplicialMesh<dim> mesh;
        const Polycrystal<dim> poly;
        GlidePlaneFactory<3> glidePlaneFactory;
        PeriodicGlidePlaneFactory<3> periodicGlidePlaneFactory;
        const std::vector<VectorDim> periodicBasis;
        const Eigen::Matrix<double,dim,Eigen::Dynamic> periodicLatticeBasis;
        const Eigen::Matrix<double,dim,Eigen::Dynamic> periodicLatticeReciprocalBasis;
        const std::vector<VectorDim> periodicShifts;
        const double EwaldLength;

        
        Eigen::VectorXd periodicCoordinates(const VectorDim& x) const;

    }; // class DislocationDynamicsBase
} // namespace model
#endif

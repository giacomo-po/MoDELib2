/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_H_
#define model_DefectiveCrystal_H_

#include <iostream>
#include <vector>
#include <memory>

#include <DislocationDynamicsBase.h>
#include <DDtimeIntegrator.h>
#include <DefectiveCrystalParameters.h>
#include <DislocationDynamicsModule.h>
#include <CrackSystem.h>
#include <UniformExternalLoadController.h>


namespace model
{
    
    template <int _dim, short unsigned int corder>
    class DefectiveCrystal
    {
        
    public:
        static constexpr int dim=_dim; // make dim available outside class
        typedef DefectiveCrystal<dim,corder> DefectiveCrystalType;
        typedef DislocationNetwork<dim,corder> DislocationNetworkType;
        typedef CrackSystem<dim> CrackSystemType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef BVPsolver<dim,2> BVPsolverType;
        typedef typename BVPsolverType::ElementType ElementType;
        
        
        DislocationDynamicsBase<dim>& ddBase;
        const std::unique_ptr<DislocationNetworkType> DN;
        const std::unique_ptr<CrackSystemType> CS;
        const std::unique_ptr<BVPsolverType> bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>> externalLoadController;
        
        
        
        static std::unique_ptr<ExternalLoadControllerBase<_dim>> getExternalLoadController(const DislocationDynamicsBase<dim>& ddBase_in, const MatrixDim& plasticStrain_in);
                
//        static std::vector<VectorDim> getPeriodicShifts(const SimplicialMesh<dim>& m,const DefectiveCrystalParameters& params);
        void updateLoadControllers(const long int& runID, const bool& isClimbStep);
        
    public:
        
        DefectiveCrystal(DislocationDynamicsBase<_dim>& ddBase_in) ;
//        DefectiveCrystal(const std::string& folderName) ;
        void singleGlideStep();;
        void runGlideSteps();
        MatrixDim plasticDistortion() const;
        MatrixDim plasticDistortionRate() const;
        MatrixDim plasticStrainRate() const;
        MatrixDim plasticStrain() const;
        double getMaxVelocity() const;
    };
}
#endif

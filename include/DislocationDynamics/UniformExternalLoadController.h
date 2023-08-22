/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UniformExternalLoadController_H_
#define model_UniformExternalLoadController_H_

#include <iostream>
#include <sstream>      // std::stringstream
#include <cmath>
#include <cfloat>

#include <ExternalLoadControllerBase.h>
#include <VoigtTraits.h>

#include <IDreader.h>



namespace model
{
    
    template <int dim>
    class UniformExternalLoadController : public ExternalLoadControllerBase<dim>
    {
        
        static constexpr int voigtSize=dim*(dim+1)/2;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        const SymmetricVoigtTraits<dim> voigtTraits;
        const bool enable;
        const int relaxSteps;
                
    public:
        
        UniformExternalLoadController(const DislocationDynamicsBase<dim>& ddBase,const MatrixDim& plasticStrain_in);
        MatrixDim stressconsidermachinestiffness(const MatrixDim& S_strain,const MatrixDim& S_stress) const override;
        MatrixDim straininducedStress(const MatrixDim& dstrain, const double& lambda) const;
        MatrixDim elasticstrain(const MatrixDim& _externalStress,const double& nu_use) const override;
        MatrixDim stress(const VectorDim&) const override;
        MatrixDim strain(const VectorDim&) const override;
        void update(const MatrixDim& plasticStrain_in) override;
        void output(std::ofstream& f_file,std::ofstream& F_labels) const override;
        static SymmetricVoigtTraits<dim> getVoigtTraits(const std::string& fileName);

    };
}
#endif

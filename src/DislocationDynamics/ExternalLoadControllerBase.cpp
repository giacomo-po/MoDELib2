/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalLoadControllerBase_cpp_
#define model_ExternalLoadControllerBase_cpp_


#include <ExternalLoadControllerBase.h>


namespace model
{
    
 
        /**************************************************************************/
        template <int dim>
        ExternalLoadControllerBase<dim>::ExternalLoadControllerBase(const DislocationDynamicsBase<dim>& ddBase_in,const std::string& _inputFileName) :
        /* init list */ ddBase(ddBase_in)
        /* init list */,inputFileName(_inputFileName)
        /* init list */,ExternalStress0(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStress0",dim,dim,true))
        /* init list */,ExternalStressRate(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStressRate",dim,dim,true))
        /* init list */,ExternalStrain0(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStrain0",dim,dim,true))
        /* init list */,ExternalStrainRate(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStrainRate",dim,dim,true))
        /* init list */,ExternalStress(MatrixDim::Zero())
        /* init list */,ExternalStrain(MatrixDim::Zero())
        /* init list */,plasticStrain(MatrixDim::Zero())
        /* init list */,MachineStiffnessRatio(TextFileParser(this->inputFileName).readMatrix<double>("MachineStiffnessRatio",1,voigtSize,true))
        /* init list */,nu(ddBase.poly.nu)
        /* init list */,lambda(2.0*nu/(1.0-2.0*nu))
        /* init list */,nu_use(nu/(1.0+nu)/2.0)
        /* init list */,stressmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero())
        /* init list */,strainmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero())
        {
            std::cout<<greenBoldColor<<"Reading ExternalLoadController file: "<<inputFileName<<defaultColor<<std::endl;
            if((ExternalStress0-ExternalStress0.transpose()).norm()>DBL_EPSILON)
            {
                throw std::runtime_error("ExternalStress0 is not symmetric.");
            }
            if((ExternalStressRate-ExternalStressRate.transpose()).norm()>DBL_EPSILON)
            {
                throw std::runtime_error("ExternalStressRate is not symmetric.");
            }
            if((ExternalStrain0-ExternalStrain0.transpose()).norm()>DBL_EPSILON)
            {
                throw std::runtime_error("ExternalStrain0 is not symmetric.");
            }
            if((ExternalStrainRate-ExternalStrainRate.transpose()).norm()>DBL_EPSILON)
            {
                throw std::runtime_error("ExternalStrainRate is not symmetric.");
            }
        }
        
        /**************************************************************************/
        template <int dim>
        ExternalLoadControllerBase<dim>::~ExternalLoadControllerBase()
        {
        }
        
        template class ExternalLoadControllerBase<3>;
}
#endif

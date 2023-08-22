/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalLoadControllerBase_H_
#define model_ExternalLoadControllerBase_H_

#include <cfloat>
#include <iostream>
#include <fstream>      // std::stringstream
#include <Eigen/Dense>
#include <TerminalColors.h>
#include <TextFileParser.h>
#include <DislocationDynamicsBase.h>
//#include <UniqueOutputFile.h>
//#include <cmath>
//#include <cfloat>
//#include <Material.h>
//#include <EigenDataReader.h>
//#include <IDreader.h>


namespace model
{

template <int dim>
class ExternalLoadControllerBase
{
    static constexpr int voigtSize=dim*(dim+1)/2;
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;
    typedef Eigen::Matrix<double,dim,1>   VectorDim;
    
    
public:
    
    const DislocationDynamicsBase<dim>& ddBase;
    const std::string inputFileName;
    
    //External stress control parameter
    
    const MatrixDim ExternalStress0;
    const MatrixDim ExternalStressRate;
    const MatrixDim ExternalStrain0;
    const MatrixDim ExternalStrainRate;
    
    
    MatrixDim ExternalStress;
    MatrixDim ExternalStrain;
    //External strain control parameter
    MatrixDim plasticStrain;
    //finite machine stiffness effect
    Eigen::Matrix<double,1,voigtSize> MachineStiffnessRatio;    //0 is stress control; infinity is pure strain control.
    
    const double nu;
    const double lambda;  //unit is mu, lambda=2*v/(1-2*v)
    const double nu_use;
    Eigen::Matrix<double,voigtSize,voigtSize> stressmultimachinestiffness;
    Eigen::Matrix<double,voigtSize,voigtSize> strainmultimachinestiffness;
    

    ExternalLoadControllerBase(const DislocationDynamicsBase<dim>& ddBase_in,const std::string& _inputFileName) ;
    virtual ~ExternalLoadControllerBase();
    virtual MatrixDim elasticstrain(const MatrixDim& _externalStress,const double& nu_use) const = 0;
    virtual MatrixDim stressconsidermachinestiffness(const MatrixDim& S_strain,const MatrixDim& S_stress) const = 0;
    virtual MatrixDim stress(const VectorDim&) const = 0;
    virtual MatrixDim strain(const VectorDim&) const = 0;
    virtual void update(const MatrixDim& plasticStrain_in) = 0;
    virtual void output(std::ofstream& f_file,std::ofstream& F_labels) const = 0;
    
};
}
#endif

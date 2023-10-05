/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectsFieldsExtractor_H_
#define model_DefectsFieldsExtractor_H_

#include <string>


#include <Eigen/Dense>

#include <DislocationDynamicsBase.h>
#include <MicrostructureGenerator.h>
#include <ConfigurationFields.h>

namespace model
{

    

    struct DefectsFieldsExtractor
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDim;
        typedef Eigen::Matrix<double,3,3> MatrixDim;

    
        DefectsFieldsExtractor(const std::string& folderName);
        
        DislocationDynamicsBase<3> ddBase;
//        DDconfigIO<3> configIO;
        MicrostructureGenerator microstructureGenerator;
        ConfigurationFields<3> configFields;
        
        const DDconfigIO<3>& config() const;
        DDconfigIO<3>& config();
        void readMicrostructure();
        void readConfiguration(const size_t& runID);
        void writeConfiguration(const size_t& runID);
        VectorDim lowerDomainCorner() const;
        VectorDim upperDomainCorner() const;
        double solidAngle(const VectorDim& x) const;
        double solidAngle(const double& x,const double& y,const double& z) const;
        VectorDim dislocationPlasticDisplacement(const VectorDim& x) const;
        VectorDim dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const;
        MatrixDim dislocationStress(const VectorDim& x) const;
        MatrixDim dislocationStress(const double& x,const double& y,const double& z) const;
        MatrixDim inclusionStress(const VectorDim& x) const;
        MatrixDim inclusionStress(const double& x,const double& y,const double& z) const;

//        double latticeParameter() const;

    };


}
#endif

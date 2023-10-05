/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityPy_cpp_
#define _model_DislocationMobilityPy_cpp_

#include <DislocationMobilityPy.h>

namespace model
{
#ifdef _MODEL_PYBIND11_ // COMPLIED WITH PYBIND11
    DislocationMobilityPy::DislocationMobilityPy(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in) :
    /* init */ DislocationMobilityBase("Py mobility for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,pyModuleName(pyModuleName_in)
    {// Set up pyModule
        //pyModule
    }

    double DislocationMobilityPy::velocity(const MatrixDim& ,
                    const VectorDim& ,
                    const VectorDim& ,
                    const VectorDim& ,
                    const double& ,
                    const double& ,
                    const double& ,
                    const std::shared_ptr<StochasticForceGenerator>& )
    {
        // use pyModule to compute velocity
        throw std::runtime_error("DislocationMobilityPy FINISH THIS");
        return 0.0;
    }

#else // COMPLIED WITHOUT PYBIND11
    DislocationMobilityPy::DislocationMobilityPy(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in) :
    /* init */ DislocationMobilityBase("Py mobility for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,pyModuleName(pyModuleName_in)
    {
    }

    double DislocationMobilityPy::velocity(const MatrixDim& ,
                    const VectorDim& ,
                    const VectorDim& ,
                    const VectorDim& ,
                    const double& ,
                    const double& ,
                    const double& ,
                    const std::shared_ptr<StochasticForceGenerator>& )
    {
        throw std::runtime_error("DislocationMobilityPy used without pybind11");
        return 0.0;
    }
#endif

   
    
}
#endif

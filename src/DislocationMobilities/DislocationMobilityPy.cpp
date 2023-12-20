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
#include <filesystem>

namespace model
{
#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    DislocationMobilityPy::DislocationMobilityPy(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in) :
    /* init */ DislocationMobilityBase("Py mobility for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,pyModuleName(pyModuleName_in)
    {// Set up pyModule
        std::filesystem::path modulePath(pyModuleName);
        const std::string moduleStem(modulePath.stem().string());
        const std::string moduleDir(modulePath.parent_path().string());
        std::cout<<"moduleStem="<<moduleStem<<std::endl;
        std::cout<<"moduleDir="<<moduleDir<<std::endl;
        pybind11::module sys = pybind11::module::import("sys");
        pybind11::list path = sys.attr("path");
        path.append(moduleDir);
        pyModule = pybind11::module::import(moduleStem.c_str());
    }

    double DislocationMobilityPy::velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& )
    {
        return pyModule.attr("velocityPy")(S,b,xi,n,T,dL,dt).cast<double>();
    }

#else // COMPILED WITHOUT PYBIND11
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

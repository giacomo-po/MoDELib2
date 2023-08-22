/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConfigurationFields_cpp_
#define model_ConfigurationFields_cpp_

#include <ConfigurationFields.h>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace model
{

ConfigurationFields::ConfigurationFields(const std::string& folderName):
/* init */ ddBase(folderName)
/* init */,microstructureGenerator(ddBase)
/* init */,configFields(ddBase,microstructureGenerator.config())
{
    
}

const DDconfigIO<3>& ConfigurationFields::config() const
{
    return microstructureGenerator.config();
}

DDconfigIO<3>& ConfigurationFields::config()
{
    return microstructureGenerator.config();
}

void ConfigurationFields::readConfiguration(const size_t& runID)
{
    config().read(runID);
    configFields.updateConfiguration();
}

void ConfigurationFields::readMicrostructure()
{
    microstructureGenerator.readMicrostructureFile();
    config().print();
    configFields.updateConfiguration();
}

void ConfigurationFields::writeConfiguration(const size_t& runID)
{
    microstructureGenerator.writeConfigFiles(runID);
}

double ConfigurationFields::solidAngle(const VectorDim& x) const
{
    return configFields.solidAngle(x);
}

double ConfigurationFields::solidAngle(const double& x,const double& y,const double& z) const
{
    return solidAngle((VectorDim()<<x,y,z).finished());
}

typename ConfigurationFields::VectorDim ConfigurationFields::dislocationPlasticDisplacement(const VectorDim& x) const
{
    return configFields.dislocationPlasticDisplacement(x);
}

typename ConfigurationFields::VectorDim ConfigurationFields::dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const
{
    return dislocationPlasticDisplacement((VectorDim()<<x,y,z).finished());
}

typename ConfigurationFields::MatrixDim ConfigurationFields::dislocationStress(const VectorDim& x) const
{
    return configFields.dislocationStress(x);
}

typename ConfigurationFields::MatrixDim ConfigurationFields::dislocationStress(const double& x,const double& y,const double& z) const
{
    return dislocationStress((VectorDim()<<x,y,z).finished());
}

typename ConfigurationFields::MatrixDim ConfigurationFields::inclusionStress(const VectorDim& x) const
{
    return configFields.inclusionStress(x);
}

typename ConfigurationFields::MatrixDim ConfigurationFields::inclusionStress(const double& x,const double& y,const double& z) const
{
    return inclusionStress((VectorDim()<<x,y,z).finished());
}

PYBIND11_MODULE(ConfigFields,m)
{
    namespace py=pybind11;
    py::class_<model::ConfigurationFields>(m,"ConfigurationFields")
        .def(py::init<const std::string&>())
        .def("readMicrostructure", &model::ConfigurationFields::readMicrostructure)
        .def("readConfiguration", &model::ConfigurationFields::readConfiguration)
        .def("writeConfiguration", &model::ConfigurationFields::writeConfiguration)
        .def("solidAngle", static_cast<double (model::ConfigurationFields::*)(const double&,const double&,const double&) const>(&model::ConfigurationFields::solidAngle))
        .def("dislocationPlasticDisplacement", static_cast<typename model::ConfigurationFields::VectorDim (model::ConfigurationFields::*)(const double&,const double&,const double&) const>(&model::ConfigurationFields::dislocationPlasticDisplacement))
        .def("dislocationStress", static_cast<typename model::ConfigurationFields::MatrixDim (model::ConfigurationFields::*)(const double&,const double&,const double&) const>(&model::ConfigurationFields::dislocationStress))
        .def("inclusionStress", static_cast<typename model::ConfigurationFields::MatrixDim (model::ConfigurationFields::*)(const double&,const double&,const double&) const>(&model::ConfigurationFields::inclusionStress))
    ;
}

}
#endif

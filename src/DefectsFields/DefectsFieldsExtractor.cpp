/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectsFieldsExtractor_cpp_
#define model_DefectsFieldsExtractor_cpp_

#include <DefectsFieldsExtractor.h>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace model
{

DefectsFieldsExtractor::DefectsFieldsExtractor(const std::string& folderName):
/* init */ ddBase(folderName)
/* init */,microstructureGenerator(ddBase)
/* init */,configFields(ddBase,microstructureGenerator.config())
{
    
}

const DDconfigIO<3>& DefectsFieldsExtractor::config() const
{
    return microstructureGenerator.config();
}

DDconfigIO<3>& DefectsFieldsExtractor::config()
{
    return microstructureGenerator.config();
}

void DefectsFieldsExtractor::readConfiguration(const size_t& runID)
{
    config().read(runID);
    configFields.updateConfiguration();
}

void DefectsFieldsExtractor::readMicrostructure()
{
    microstructureGenerator.readMicrostructureFile();
    config().print();
    configFields.updateConfiguration();
}

void DefectsFieldsExtractor::writeConfiguration(const size_t& runID)
{
    microstructureGenerator.writeConfigFiles(runID);
}

double DefectsFieldsExtractor::solidAngle(const VectorDim& x) const
{
    return configFields.solidAngle(x);
}

double DefectsFieldsExtractor::solidAngle(const double& x,const double& y,const double& z) const
{
    return solidAngle((VectorDim()<<x,y,z).finished());
}

typename DefectsFieldsExtractor::VectorDim DefectsFieldsExtractor::dislocationPlasticDisplacement(const VectorDim& x) const
{
    return configFields.dislocationPlasticDisplacement(x);
}

typename DefectsFieldsExtractor::VectorDim DefectsFieldsExtractor::dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const
{
    return dislocationPlasticDisplacement((VectorDim()<<x,y,z).finished());
}

typename DefectsFieldsExtractor::MatrixDim DefectsFieldsExtractor::dislocationStress(const VectorDim& x) const
{
    return configFields.dislocationStress(x);
}

typename DefectsFieldsExtractor::MatrixDim DefectsFieldsExtractor::dislocationStress(const double& x,const double& y,const double& z) const
{
    return dislocationStress((VectorDim()<<x,y,z).finished());
}

typename DefectsFieldsExtractor::MatrixDim DefectsFieldsExtractor::inclusionStress(const VectorDim& x) const
{
    return configFields.inclusionStress(x);
}

typename DefectsFieldsExtractor::MatrixDim DefectsFieldsExtractor::inclusionStress(const double& x,const double& y,const double& z) const
{
    return inclusionStress((VectorDim()<<x,y,z).finished());
}

typename DefectsFieldsExtractor::VectorDim DefectsFieldsExtractor::lowerDomainCorner() const
{
    return ddBase.mesh.xMin();
}

typename DefectsFieldsExtractor::VectorDim DefectsFieldsExtractor::upperDomainCorner() const
{
    return ddBase.mesh.xMax();
}


PYBIND11_MODULE(DefectsFields,m)
{
    namespace py=pybind11;
    py::class_<model::DefectsFieldsExtractor>(m,"DefectsFieldsExtractor")
        .def(py::init<const std::string&>())
        .def("readMicrostructure", &model::DefectsFieldsExtractor::readMicrostructure)
        .def("readConfiguration", &model::DefectsFieldsExtractor::readConfiguration)
        .def("writeConfiguration", &model::DefectsFieldsExtractor::writeConfiguration)
        .def("lowerDomainCorner", &model::DefectsFieldsExtractor::lowerDomainCorner)
        .def("upperDomainCorner", &model::DefectsFieldsExtractor::upperDomainCorner)
        .def("solidAngle", static_cast<double (model::DefectsFieldsExtractor::*)(const double&,const double&,const double&) const>(&model::DefectsFieldsExtractor::solidAngle))
        .def("dislocationPlasticDisplacement", static_cast<typename model::DefectsFieldsExtractor::VectorDim (model::DefectsFieldsExtractor::*)(const double&,const double&,const double&) const>(&model::DefectsFieldsExtractor::dislocationPlasticDisplacement))
        .def("dislocationStress", static_cast<typename model::DefectsFieldsExtractor::MatrixDim (model::DefectsFieldsExtractor::*)(const double&,const double&,const double&) const>(&model::DefectsFieldsExtractor::dislocationStress))
        .def("inclusionStress", static_cast<typename model::DefectsFieldsExtractor::MatrixDim (model::DefectsFieldsExtractor::*)(const double&,const double&,const double&) const>(&model::DefectsFieldsExtractor::inclusionStress))
    ;
}

}
#endif

/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PyGlideSolver_cpp_
#define model_PyGlideSolver_cpp_

#include <filesystem>
#include <PyGlideSolver.h>
#include <TerminalColors.h>


namespace model
{
    

        
    #ifdef _MODEL_PYBIND11_ // COMPLIED WITH PYBIND11

    template <typename DislocationNetworkType>
    PyGlideSolver<DislocationNetworkType>::PyGlideSolver(const DislocationNetworkType& DN_in) :
    /* init */ DislocationGlideSolverBase<DislocationNetworkType>(DN_in)
    /* init */,pyModuleName(this->DN.simulationParameters.traitsIO.inputFilesFolder+"/"+TextFileParser(this->DN.simulationParameters.traitsIO.ddFile).readString("pyModuleName"))
    {
        std::filesystem::path pyModulePath(pyModuleName);
        std::string pyModuleDir(pyModulePath.parent_path());
        std::string pyModuleFile(pyModulePath.filename());
        std::cout<<greenBoldColor<<"Creating PyGlideSolver: dir= "<<pyModuleDir<<", file="<<pyModuleFile<<defaultColor<<std::endl;
        pybind11::scoped_interpreter guard{};
        pybind11::module_ sys = pybind11::module_::import("sys");
        pybind11::list path = sys.attr("path");
        path.append(pyModuleDir.c_str());

        pyModule=pybind11::module_::import(pyModuleFile.c_str());
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd PyGlideSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"PyGlideSolver solving "<<std::flush;
        Eigen::VectorXd nodeVelocities;
        // compute nodeVelocities

        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        return nodeVelocities;
    }

    #else // COMPLIED WITHOUT PYBIND11

    template <typename DislocationNetworkType>
    PyGlideSolver<DislocationNetworkType>::PyGlideSolver(const DislocationNetworkType& DN_in) :
    /* init */ DislocationGlideSolverBase<DislocationNetworkType>(DN_in)
    {
        std::cout<<greenBoldColor<<"Creating PyGlideSolver"<<defaultColor<<std::endl;
        throw std::runtime_error("PyGlideSolver must be compiled with pybind11");
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd PyGlideSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        return Eigen::VectorXd();
    }
    #endif


    template class PyGlideSolver<DislocationNetwork<3,0>>;
    
}
#endif

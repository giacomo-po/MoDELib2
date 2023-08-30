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

#ifndef model_DislocationGlideSolverFactory_cpp_
#define model_DislocationGlideSolverFactory_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationGlideSolverFactory.h>
#include <GalerkinGlideSolver.h>
#include <PyGlideSolver.h>

namespace model
{
    
    template <typename DislocationNetworkType>
    DislocationGlideSolverBase<DislocationNetworkType>::DislocationGlideSolverBase(const DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    {
        
    }

template <typename DislocationNetworkType>
std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>> DislocationGlideSolverFactory<DislocationNetworkType>::getGlideSolver(const DislocationNetworkType& DN,const std::string& solverType)
{
        if(solverType=="Galerkin" || solverType=="galerkin")
        {
            return std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>>(new GalerkinGlideSolver<DislocationNetworkType>(DN));
        }
        else if(solverType=="Pybind11" || solverType=="pybind11")
        {
            return std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>>(new PyGlideSolver<DislocationNetworkType>(DN));
        }
    else
    {
        throw std::runtime_error("Unknown glide solver type '"+solverType+"'");
        return nullptr;
    }
    
}
    
    template struct DislocationGlideSolverBase<DislocationNetwork<3,0>>;
    template struct DislocationGlideSolverFactory<DislocationNetwork<3,0>>;

}
#endif

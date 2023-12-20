/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SecondPhase_cpp_
#define model_SecondPhase_cpp_

#include <memory>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <SecondPhase.h>
#include <TerminalColors.h>

namespace model
{
    
    template<int dim>
    SecondPhase<dim>::SecondPhase(const std::string& _name,
                                  const std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>>& _gsMap) :
    /* init */ name(_name)
    /* init */,gsMap(_gsMap)
    {
        
        std::cout<<greenBoldColor<<"Creating SecondPhase "<<name<<", phaseID= "<<this->sID<<defaultColor<<std::endl;
        
    }

    template<int dim>
    double SecondPhase<dim>::misfitEnergy(const VectorDim& b,const GlidePlaneBase* const gpb) const
    {
        
//        std::cout<<"inPlane="<<n<<std::endl;
//
//        std::cout<<"planes"<<std::endl;
//        for(const auto& pair : gsMap)
//        {
//            std::cout<<pair.first<<std::endl;
//        }
        
        const auto gammaIter(gsMap.find(gpb));
        if(gammaIter!=gsMap.end())
        {
            const VectorDim bL(gpb->G2L*b);
            if(std::fabs(bL(2))>FLT_EPSILON)
            {
                throw std::runtime_error("SLIP VECTOR NOT ON GlidePlaneBase");
            }
            return gammaIter->second->operator()(bL.template segment<2>(0));
        }
        else
        {
            return 0.0;
        }

    }

    
template struct SecondPhase<3>;
}
#endif

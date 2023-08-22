/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationFieldBase_h_
#define _model_DislocationFieldBase_h_

#include <Eigen/Dense>

namespace model
{
    
    
    template<int dim>
    struct DislocationFieldBase
    {
        
        //! Dislocation core size
        inline static  double a=1;
        
        //! Dislocation core size squared
        inline static  double a2=1;
        
        //! dim x dim identity matrix
        inline static const Eigen::Matrix<double,dim,dim> I=Eigen::Matrix<double,dim,dim>::Identity();
        static void initFromFile(const std::string& fileName);
        
        
    };
    
    
    //! Dislocation core size
//    template<int dim>
//    double DislocationFieldBase<dim>::a=1.0;
//
//    //! Dislocation core size squared
//    template<int dim>
//    double DislocationFieldBase<dim>::a2=1.0;
//
//    //! Identity matrix
//    template<int dim>
//    const Eigen::Matrix<double,dim,dim> DislocationFieldBase<dim>::I=Eigen::Matrix<double,dim,dim>::Identity();
    
}	// close namespace
#endif




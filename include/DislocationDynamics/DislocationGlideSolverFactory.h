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

#ifndef model_DislocationGlideSolverFactory_h_
#define model_DislocationGlideSolverFactory_h_

#include <memory>
#include <Eigen/Dense>

namespace model
{
    
    template <typename DislocationNetworkType>
    struct DislocationGlideSolverBase
    {
        const DislocationNetworkType& DN;
        DislocationGlideSolverBase(const DislocationNetworkType& );
        virtual Eigen::VectorXd getNodeVelocities() const = 0;
    };

    template <typename DislocationNetworkType>
    struct DislocationGlideSolverFactory
    {
        static std::shared_ptr<DislocationGlideSolverBase<DislocationNetworkType>> getGlideSolver(const DislocationNetworkType& DN, const std::string& solverType);
    };
    
}
#endif

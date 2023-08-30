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

#ifndef model_GalerkinGlideSolver_h_
#define model_GalerkinGlideSolver_h_

//#include <PlanarDislocationNode.h>


#ifndef NDEBUG
#define VerboseGalerkinGlideSolver(N,x) if(verboseGalerkinGlideSolver>=N){std::cout<<x;}
#else
#define VerboseGalerkinGlideSolver(N,x)
#endif

#include <DislocationDynamicsModule.h>
#include <DislocationGlideSolverFactory.h>
namespace model
{
    
    template <typename DislocationNetworkType>
//    template<int dim,int corder>
    class GalerkinGlideSolver : public DislocationGlideSolverBase<DislocationNetworkType>
    {
//        typedef DislocationNetwork<dim,corder> DislocationNetworkType;
        typedef typename DislocationNetworkType::VectorDim   VectorDim;
        typedef typename DislocationNetworkType::MatrixDim   MatrixDim;
        typedef typename DislocationNetworkType::NetworkNodeType   NetworkNodeType;
        typedef typename DislocationNetworkType::LoopNodeType   LoopNodeType;
        typedef typename DislocationNetworkType::LoopType   LoopType;

        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        //        typedef std::vector<Eigen::Triplet<double> > TripletContainerType;
        typedef std::deque<Eigen::Triplet<double> > TripletContainerType;
//        enum {NdofXnode=NodeType::NdofXnode};
        constexpr static int NdofXnode=NetworkNodeType::NdofXnode;
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;


//        void lumpedSolve();
//        void solve(const size_t& runID);
        Eigen::VectorXd lumpedSolve() const;
        size_t assembleNCtriplets(TripletContainerType& kqqT, Eigen::VectorXd& Fq) const;
//        void storeNodeSolution(const Eigen::VectorXd& X);
        size_t assembleConstraintsforPeriodicSimulationsNULL(TripletContainerType& zT) const;

        
        public:
        
        GalerkinGlideSolver(const DislocationNetworkType& );
        Eigen::VectorXd getNodeVelocities() const override;

        
    };
    
}
#endif

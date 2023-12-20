/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui  <cuiyinan@g.ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_cpp_
#define model_DislocationNodeContraction_cpp_

#include <memory>
#include <tuple>
#include <Eigen/Dense>
#include <TypeTraits.h>
#include <GlidePlane.h>
//#include <BoundingMeshSegments.h>
#include <ConfinedDislocationObject.h>
#include <Grain.h>
#include <FiniteLineSegment.h>
#include <PlanesIntersection.h>
#include <DislocationNetwork.h>
#include <DislocationNodeContraction.h>

#ifndef NDEBUG
#define VerboseNodeContraction(N,x) if(verboseNodeContraction>=N){std::cout<<x;}
#else
#define VerboseNodeContraction(N,x)
#endif

namespace model
{

    template <typename DislocationNetworkType>
    DislocationNodeContraction<DislocationNetworkType>::DislocationNodeContraction(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    /* init */,verboseNodeContraction(TextFileParser(DN.ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseNodeContraction",true))
    {
        
    }

    template <typename DislocationNetworkType>
    std::vector<std::set<std::shared_ptr<typename DislocationNodeContraction<DislocationNetworkType>::NetworkNodeType>>> DislocationNodeContraction<DislocationNetworkType>::contractBoundary(std::shared_ptr<NetworkNodeType> nA,
                                                                                                                                                                  std::shared_ptr<NetworkNodeType> nB)
    {
        std::vector<std::set<std::shared_ptr<NetworkNodeType>>> bndNodesToContract;
        
        //See if you want to implement contract network node based on sIDs .. this will greatly help with contracting younger nodes in contract boundary nodes
        
        if (DN.ddBase.simulationParameters.isPeriodicSimulation())
        {
            VerboseNodeContraction(1,"Populating boundary nodes to be contracted"<<std::endl;);
            
            for (const auto &nALN : nA->loopNodes())
            {
                const auto nApPrev(nALN->periodicPrev());
                const auto nApNext(nALN->periodicNext());
                for (const auto &nBLN : nB->loopNodes())
                {
                    const auto nBpPrev(nBLN->periodicPrev());
                    const auto nBpNext(nBLN->periodicNext());
                    if ((nApPrev && nApNext) && (nBpPrev && nBpNext))
                    {
                        if (nApPrev->networkNode==nBpPrev->networkNode)
                        {
                            VerboseNodeContraction(1,"Case a"<<std::endl;);
                            // std::cout<<"Prev networkNOde "<<nApPrev->networkNode->tag()<<std::endl;
                            // std::cout<<" bnd next prev size "<<nALN->boundaryPrev().size()<<" "<<nBLN->boundaryPrev().size()<<std::endl;
                            //bnd prev needs to be contracted
                            for (const auto& bndPrevA : nALN->boundaryPrev())
                            {
                                for (const auto& bndPrevB : nBLN->boundaryPrev())
                                {
                                    // std::cout<<" bndPrevA->networkNode->tag()=>bndPrevB->networkNode->tag() "<<bndPrevA->networkNode->tag()
                                    // <<bndPrevB->networkNode->tag()<<std::endl;
                                    
                                    // if (bndPrevA->networkNode->isMovableTo(bndPrevB->networkNode->get_P()))
                                    if (bndPrevA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces() &&
                                        (bndPrevA->networkNode->get_P()-bndPrevB->networkNode->get_P()).norm()<FLT_EPSILON)
                                    {
                                        const auto minNode(std::min(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
                                                           bndPrevA->networkNode : bndPrevB->networkNode);
                                        const auto maxNode(std::max(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
                                                           bndPrevA->networkNode : bndPrevB->networkNode);
                                        //search in the container if the nodes are already present.. if yes place in that location
                                        bool nodeAlreadyPresentInSet(false);
                                        if (bndNodesToContract.size()==0)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        for (auto& setIter : bndNodesToContract)
                                        {
                                            if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                            {
                                                setIter.insert(minNode);
                                                setIter.insert(maxNode);
                                                nodeAlreadyPresentInSet=true;
                                                break;
                                            }
                                        }
                                        if (!nodeAlreadyPresentInSet)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                    }
                                }
                            }
                        }
                        if (nApNext->networkNode==nBpNext->networkNode)
                        {
                            VerboseNodeContraction(1,"Case b"<<std::endl;);
                            
                            //bnd next needs to be contracted
                            for (const auto &bndNextA : nALN->boundaryNext())
                            {
                                for (const auto &bndNextB : nBLN->boundaryNext())
                                {
                                    if (bndNextA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces() &&
                                        (bndNextA->networkNode->get_P()-bndNextB->networkNode->get_P()).norm()<FLT_EPSILON)
                                    {
                                        const auto minNode (std::min(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
                                                            bndNextA->networkNode : bndNextB->networkNode);
                                        const auto maxNode (std::max(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
                                                            bndNextA->networkNode : bndNextB->networkNode);
                                        //search in the container if the nodes are already present.. if yes place in that location
                                        bool nodeAlreadyPresentInSet(false);
                                        if (bndNodesToContract.size()==0)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        for (auto& setIter : bndNodesToContract)
                                        {
                                            if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                            {
                                                setIter.insert(minNode);
                                                setIter.insert(maxNode);
                                                nodeAlreadyPresentInSet=true;
                                                break;
                                            }
                                        }
                                        if (!nodeAlreadyPresentInSet)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        // bndNodesToContract.insert(std::make_pair(minNode,maxNode));
                                    }
                                }
                            }
                        }
                        if (nApPrev->networkNode==nBpNext->networkNode)
                        {
                            VerboseNodeContraction(1,"Case c"<<std::endl;);
                            
                            //bnd prev of A needs to be contracted to bndNext of B
                            for (const auto &bndPrevA : nALN->boundaryPrev())
                            {
                                for (const auto &bndNextB : nBLN->boundaryNext())
                                {
                                    if (bndPrevA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces() &&
                                        (bndPrevA->networkNode->get_P()-bndNextB->networkNode->get_P()).norm()<FLT_EPSILON)
                                    {
                                        const auto minNode(std::min(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
                                                           bndPrevA->networkNode : bndNextB->networkNode);
                                        const auto maxNode(std::max(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
                                                           bndPrevA->networkNode : bndNextB->networkNode);
                                        //search in the container if the nodes are already present.. if yes place in that location
                                        bool nodeAlreadyPresentInSet(false);
                                        if (bndNodesToContract.size()==0)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        for (auto& setIter : bndNodesToContract)
                                        {
                                            if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                            {
                                                setIter.insert(minNode);
                                                setIter.insert(maxNode);
                                                nodeAlreadyPresentInSet=true;
                                            }
                                        }
                                        if (!nodeAlreadyPresentInSet)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        // bndNodesToContract.insert(std::make_pair(minNode,maxNode));
                                    }
                                }
                            }
                        }
                        if (nApNext->networkNode==nBpPrev->networkNode)
                        {
                            VerboseNodeContraction(1,"Case d"<<std::endl;);
                            
                            //bnd next of A needs to be contracted to bndPrev of B
                            for (const auto &bndNextA : nALN->boundaryNext())
                            {
                                for (const auto &bndPrevB : nBLN->boundaryPrev())
                                {
                                    if (bndNextA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces() &&
                                        (bndNextA->networkNode->get_P()-bndPrevB->networkNode->get_P()).norm()<FLT_EPSILON)
                                    {
                                        const auto minNode(std::min(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
                                                           bndNextA->networkNode : bndPrevB->networkNode);
                                        const auto maxNode(std::max(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
                                                           bndNextA->networkNode : bndPrevB->networkNode);
                                        //search in the container if the nodes are already present.. if yes place in that location
                                        bool nodeAlreadyPresentInSet(false);
                                        if (bndNodesToContract.size()==0)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        for (auto& setIter : bndNodesToContract)
                                        {
                                            if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                            {
                                                setIter.insert(minNode);
                                                setIter.insert(maxNode);
                                                nodeAlreadyPresentInSet=true;
                                            }
                                        }
                                        if (!nodeAlreadyPresentInSet)
                                        {
                                            std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                            nodeSet.insert(minNode);
                                            nodeSet.insert(maxNode);
                                            bndNodesToContract.push_back(nodeSet);
                                        }
                                        // bndNodesToContract.insert(std::make_pair(minNode,maxNode));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return bndNodesToContract;
        
    }

    template <typename DislocationNetworkType>
    bool DislocationNodeContraction<DislocationNetworkType>::contractSecondAndVirtual(std::shared_ptr<NetworkNodeType> nA,
                                                                                      std::shared_ptr<NetworkNodeType> nB)
    {
        
        if ((nB->get_P()-nA->get_P()).squaredNorm()>FLT_EPSILON)
        {
            // std::cout<<"Na and Nb are "<<nA->tag()<<"=>"<<nB->tag()<<std::endl;
            assert(nB->isMovableTo(nA->get_P()));
            bool movedB(nB->trySet_P(nA->get_P()));
            if(!movedB)
            {
                throw std::runtime_error("Node B could not be moved to position.");
            }
            //                assert(movedB);
        }
        //            switch (DN.simulationParameters.simulationType)
        //            {
        //                case DefectiveCrystalParameters::FINITE_FEM:
        //                {
        //                    if(nB->virtualBoundaryNode())
        //                    {
        //                        assert(nA->virtualBoundaryNode());
        //                        DN.contractNetworkNodes(nA->virtualBoundaryNode(),nB->virtualBoundaryNode());
        //                    }
        //                    break;
        //                }
        //
        ////                case DefectiveCrystalParameters::PERIODIC:
        ////                {
        ////                 // FINISH HERE
        ////                    break;
        ////                }
        //            }
        
        auto bndNodesToContract(contractBoundary(nA,nB));
        
        bool tempContract(DN.contractNetworkNodes(nA,nB));  //First contract the internal nodes then the boundary (Essential for cutLoop)
        //Now contract the boundary nodes
        for (const auto& vecIter : bndNodesToContract)
        {
            for (const auto &setIter : vecIter)
            {
                if (setIter!= *vecIter.begin())
                {
                    VerboseNodeContraction(1, "Contracting Boundary"<<(*vecIter.begin())->sID<<", "<<setIter->sID << std::endl;);
                    DN.contractNetworkNodes(*vecIter.begin(),setIter);
                }
            }
        }
        return tempContract;
    }

    template <typename DislocationNetworkType>
    bool DislocationNodeContraction<DislocationNetworkType>::contractYoungest(std::shared_ptr<NetworkNodeType> nA,
                                                                              std::shared_ptr<NetworkNodeType> nB)
    {
        return nA->sID<nB->sID? contractSecondAndVirtual(nA,nB) : contractSecondAndVirtual(nB,nA);
        
    }

    template <typename DislocationNetworkType>
    bool DislocationNodeContraction<DislocationNetworkType>::contractToPosition(std::shared_ptr<NetworkNodeType> nA,
                                                                                std::shared_ptr<NetworkNodeType> nB,
                                                                                const VectorDim& X,
                                                                                const double& maxRange)
    {
        
        bool movedA=false;
        bool movedB=false;
        
        if(   nA->isMovableTo(X)
           && nB->isMovableTo(X)
           && (nA->get_P()-X).norm()+(nB->get_P()-X).norm()<maxRange)
        {
            // std::cout<<" Before motion positions are "<<nA->get_P().transpose()<<"=>"<<nB->get_P().transpose()<<std::endl;
            movedA=nA->trySet_P(X);
            movedB=nB->trySet_P(X);
            // std::cout<<" After motion positions are "<<nA->get_P().transpose()<<"=>"<<nB->get_P().transpose()<<std::endl;
            // std::cout<<" X is "<<X.transpose()<<std::endl;
            // movedA=nA->set_P(X);
            // movedB=nB->set_P(X);
            VerboseNodeContraction(2,"contractToPosition"<<std::endl;);
            VerboseNodeContraction(2,"movedA="<<movedA<<std::endl;);
            VerboseNodeContraction(2,"movedB="<<movedB<<std::endl;);
            assert(movedA && movedB && "COULD NOT MOVE NODES");
        }
        
        return (movedA && movedB)? contractYoungest(nA,nB) : false;
    }



    template <typename DislocationNetworkType>
    bool DislocationNodeContraction<DislocationNetworkType>::contract(std::shared_ptr<NetworkNodeType> nA,
                                                                      std::shared_ptr<NetworkNodeType> nB)
    {
        VerboseNodeContraction(1, "DislocationNodeContraction::contract " << nA->sID << " " << nB->sID << std::endl;);
        if (   !nA->isBoundaryNode()
            && !nA->isGrainBoundaryNode()
            && !nB->isBoundaryNode()
            && !nB->isGrainBoundaryNode()
            ) //Added for avoiding node contraction with boundary node
        {
            
            const bool nAisMovable = nA->isMovableTo(nB->get_P());
            const bool nBisMovable = nB->isMovableTo(nA->get_P());
            
            if (nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 1a" << std::endl;);
                return contractYoungest(nA, nB);
            }
            else if (nAisMovable && !nBisMovable)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 1b" << std::endl;);
                return contractSecondAndVirtual(nB, nA);
            }
            else if (!nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 1c" << std::endl;);
                return contractSecondAndVirtual(nA, nB);
            }
            else
            {// nA and nB cannot be moved to each other. The calculation of a third point is necessary
                VerboseNodeContraction(1, "DislocationNodeContraction case 1d" << std::endl;);
                const auto planesA(nA->glidePlanes());
                const auto planesB(nB->glidePlanes());
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> N(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,planesA.size()+planesB.size()));
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,planesA.size()+planesB.size()));
                
                int k=0;
                for(const auto& plane : planesA)
                {
                    N.col(k)=plane->unitNormal;
                    P.col(k)=plane->P;
                    ++k;
                }
                for(const auto& plane : planesB)
                {
                    N.col(k)=plane->unitNormal;
                    P.col(k)=plane->P;
                    ++k;
                }
                
                PlanesIntersection<dim> pInt(N,P,FLT_EPSILON);
                const std::pair<bool,VectorDim> snapped(pInt.snap(0.5*(nA->get_P()+nB->get_P())));
                if(snapped.first)
                {
                    const double maxRange = 4.0 * (nA->get_P() - nB->get_P()).norm();
                    return contractToPosition(nA,nB,snapped.second,maxRange);
                }
                else
                {
                    return false;
                }
            }
        }
        else
        {
            if ((nA->get_P() - nB->get_P()).squaredNorm() < FLT_EPSILON)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 2a" << std::endl;);
                return contractYoungest(nA, nB);
            }
            else
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 2b" << std::endl;);
                return false;
            }
        }
    }

template class DislocationNodeContraction<DislocationNetwork<3,0>>;
}
#endif


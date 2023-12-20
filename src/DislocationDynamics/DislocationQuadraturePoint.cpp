/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

//Implements non-uniform open quadrature rule

#ifndef model_DislocationQuadraturePoint_cpp_
#define model_DislocationQuadraturePoint_cpp_

#include <DislocationQuadraturePoint.h>
#include <DislocationSegment.h>


namespace model
{

    template<int dim,int corder>
    DislocationQuadraturePoint<dim,corder>::DislocationQuadraturePoint(const LinkType& parentSegment,
                                                                       const int& q,const int& qOrder,
                                                                       const MatrixNcoeff& SFCH,
                                                                       const MatrixNcoeffDim& qH) :
    /* init */ sourceID(parentSegment.source->sID)
    /* init */,sinkID(parentSegment.sink->sID)
    /* init */,qID(q)
    /* init */,SF(QuadPowDynamicType::uPow(qOrder).row(qID)*SFCH)
    /* init */,r(SF*qH)
    /* init */,ru(QuadPowDynamicType::duPow(qOrder).row(qID)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH)
    /* init */,j(ru.norm())
    /* init */,rl(j > FLT_EPSILON ? (ru/j).eval() : VectorDim::Zero())
    /* init */,dL(j*QuadratureDynamicType::weight(qOrder,qID))
    /* init */,stress(MatrixDim::Zero())
    /* init */,pkForce(VectorDim::Zero())
    /* init */,stackingFaultForce(VectorDim::Zero())
    /* init */,lineTensionForce(VectorDim::Zero())
    /* init */,glideVelocity(VectorDim::Zero())
    /* init */,elasticEnergyPerLength(0.0)
    /* init */,coreEnergyPerLength(0.0)
    /* init */,inclusionID(-1)
    /* init */,slipVectors(std::make_pair(VectorDim::Zero(),VectorDim::Zero()))
    {
        for(const auto& inclusion : parentSegment.network().eshelbyInclusions() )
        {// Add EshelbyInclusions stress
            if(inclusion.second->contains(r))
            {
                inclusionID=inclusion.second->sID;
                break;
            }
        }
    }

    template<int dim,int corder>
    DislocationQuadraturePoint<dim,corder>::DislocationQuadraturePoint(const size_t sourceID_in,
                                                                       const size_t sinkID_in,
                                                                       const int& q,const int& qOrder,
                                                                       const MatrixNcoeff& SFCH,
                                                                       const MatrixNcoeffDim& qH) :
    /* init */ sourceID(sourceID_in)
    /* init */,sinkID(sinkID_in)
    /* init */,qID(q)
    /* init */,SF(QuadPowDynamicType::uPow(qOrder).row(qID)*SFCH)
    /* init */,r(SF*qH)
    /* init */,ru(QuadPowDynamicType::duPow(qOrder).row(qID)*SFCH.template block<Ncoeff-1,Ncoeff>(1,0)*qH)
    /* init */,j(ru.norm())
    /* init */,rl(j > FLT_EPSILON ? (ru/j).eval() : VectorDim::Zero())
    /* init */,dL(j*QuadratureDynamicType::weight(qOrder,qID))
    /* init */,stress(MatrixDim::Zero())
    /* init */,pkForce(VectorDim::Zero())
    /* init */,stackingFaultForce(VectorDim::Zero())
    /* init */,lineTensionForce(VectorDim::Zero())
    /* init */,glideVelocity(VectorDim::Zero())
    /* init */,elasticEnergyPerLength(0.0)
    /* init */,coreEnergyPerLength(0.0)
    /* init */,inclusionID(-1)
    /* init */,slipVectors(std::make_pair(VectorDim::Zero(),VectorDim::Zero()))
    {
        
    }


    template<int dim,int corder>
    DislocationQuadraturePoint<dim,corder>::DislocationQuadraturePoint() :
    /* init */ sourceID(0)
    /* init */,sinkID(0)
    /* init */,qID(0)
    /* init */,SF(Eigen::Matrix<double,1,Ncoeff>::Zero())
    /* init */,r(VectorDim::Zero())
    /* init */,ru(VectorDim::Zero())
    /* init */,j(0.0)
    /* init */,rl(VectorDim::Zero())
    /* init */,dL(0.0)
    /* init */,stress(MatrixDim::Zero())
    /* init */,pkForce(VectorDim::Zero())
    /* init */,stackingFaultForce(VectorDim::Zero())
    /* init */,lineTensionForce(VectorDim::Zero())
    /* init */,glideVelocity(VectorDim::Zero())
    /* init */,elasticEnergyPerLength(0.0)
    /* init */,coreEnergyPerLength(0.0)
    /* init */,inclusionID(-1)
    /* init */,slipVectors(std::make_pair(VectorDim::Zero(),VectorDim::Zero()))
    {
        
    }


    //
    //        template <class T>
    //        friend T& operator << (T& os, const DislocationQuadraturePoint<dim,corder>& dqp)
    //        {
    //            os  << dqp.sourceID<<" "<<dqp.sinkID<<" "<<dqp.qID<<" "
    //            /**/<< std::setprecision(15)<<std::scientific<< dqp.SF<<" "
    //            /**/<< dqp.r.transpose()<<" "
    //            /**/<< dqp.ru.transpose()<<" "
    //            /**/<< dqp.j<<" "
    //            /**/<< dqp.pkForce.transpose()<<" "
    //            /**/<< dqp.elasticEnergyPerLength;
    //            return os;
    //        }


    template<int dim,int corder>
    typename DislocationQuadraturePoint<dim,corder>::MatrixDimNdof DislocationQuadraturePoint<dim,corder>::SFgaussEx() const
    { /*! The MatrixDimNdof matrix of shape functions at the k-th quadrature point
       */
        MatrixDimNdof temp(MatrixDimNdof::Zero());
        for (size_t n=0;n<Ncoeff;++n)
        {
            temp.template block<dim,dim>(0,n*dim)=MatrixDim::Identity()*SF(n);
        }
        return temp;
    }


    //        template<typename LinkType>
    template<int dim,int corder>
    typename DislocationQuadraturePoint<dim,corder>::VectorDim  DislocationQuadraturePoint<dim,corder>::getGlideVelocity(const LinkType& parentSegment,
                                                                                                                         const VectorDim& , // position x
                                                                                                                         const VectorDim& fPK,
                                                                                                                         const MatrixDim& S,
                                                                                                                         const VectorDim& rl,
                                                                                                                         const double& dL,
                                                                                                                         const int& inclusionID
                                                                                                                         )
    {
        // std::cout<<"parentSegment "<<parentSegment.tag()<<" has slipSystem Compare to nullPtr"<<(parentSegment.slipSystem()==nullptr)<<std::endl;
        
        if(parentSegment.slipSystem())
        {
            VectorDim n(parentSegment.glidePlaneNormal()); // plane normal
            VectorDim b(parentSegment.burgers()); // Burgers vector
            VectorDim t(rl);            // tangent vector
            
            // Select right-handed normal whenever possible
            if(parentSegment.loopLinks().size()==1)
            {// pick right-handed normal for n
                const typename LinkType::LoopLinkType& loopLink(**parentSegment.loopLinks().begin());
                if(std::fabs(loopLink.loop->slippedArea())>FLT_EPSILON)
                {
                    if(parentSegment.source->sID!=loopLink.source->sID)
                    {// NetworkLink and LoopLink are oriented in opposite direction
                        b*=-1.0;
                        t*=-1.0;
                    }
                }
            }
            
            VectorDim glideForce = fPK-fPK.dot(n)*n;
            double glideForceNorm(glideForce.norm());
            //        std::cout<<"fPK="<<fPK.transpose()<<std::endl;
            //        std::cout<<"glideForceNorm="<<glideForceNorm<<std::endl;
            
            if(glideForceNorm<FLT_EPSILON && parentSegment.network().stochasticForceGenerator)
            {
                glideForce=parentSegment.chord().cross(n);
                glideForceNorm=glideForce.norm();
                if(glideForceNorm>FLT_EPSILON)
                {
                    glideForce/=glideForceNorm;
                }
            }
            
            VectorDim vv=VectorDim::Zero();
            if(glideForceNorm>FLT_EPSILON)
            {
                //            std::cout<<"a"<<std::endl;
                
                //                    double v =parentSegment.network().poly.mobility->velocity(S,b,t,n,
                double v =parentSegment.slipSystem()->mobility->velocity(S,b,t,n,
                                                                         parentSegment.network().ddBase.poly.T,
                                                                         dL,parentSegment.network().ddBase.simulationParameters.dt,parentSegment.network().stochasticForceGenerator);
                // std::cout<<"v is "<<v<<std::endl;
                if(v<0.0 && v>=-FLT_EPSILON)
                {
                    v=0.0; // kill roundoff errors for small negative velocities
                }
                
                if(inclusionID>=0)
                {
                    v*=parentSegment.network().eshelbyInclusions().at(inclusionID)->mobilityReduction;
                }
                
                assert((parentSegment.network().stochasticForceGenerator || v>= 0.0) && "Velocity must be a positive scalar");
                const bool useNonLinearVelocity=true;
                if(useNonLinearVelocity && v>FLT_EPSILON)
                {
                    v= 1.0-std::exp(-v);
                }
                
                //                    for(const auto& inclusion : parentSegment.network().eshelbyInclusions() )
                //                    {// Add EshelbyInclusions stress
                //                        if(inclusion.second.contains(r))
                //                        {
                //                            v*=inclusion.second.mobilityReduction;
                //                        }
                //                    }
                
                
                vv= v * glideForce/glideForceNorm;
                // std::cout<<"vv is "<<vv.transpose()<<std::endl;
                
            }
            return vv;
            
        }
        else
        {
            return VectorDim::Zero();
        }
        
    }

    //
    //        std::set<const EshelbyInclusion<dim>*> getInclusions(const std::map<size_t,EshelbyInclusion<dim>>& inclusionContainer,
    //                                                            const VectorDim& x)
    //        {
    //            std::set<const EshelbyInclusion<dim>*> temp;
    //            for(const auto& pair : inclusionContainer)
    //            {
    //                if(pair.second.cointains(x))
    //                {
    //                    temp.insert(&pair.second);
    //                }
    //            }
    //            return temp;
    //        }



    //        template<typename LinkType>
    template<int dim,int corder>
    void DislocationQuadraturePoint<dim,corder>::updateForcesAndVelocities(const LinkType& parentSegment)
    {
        
        const bool useSplineTangents(false);
        if(parentSegment.loopLinks().size()==1 && useSplineTangents)
        {
            const auto ll(*parentSegment.loopLinks().begin());
            const VectorDim prevNodePos(ll->prev->twin()? ll->prev->twin()->source->periodicPrev()->get_P()
                                        -ll->source->periodicPlanePatch()->shift
                                        +ll->prev->twin()->source->periodicPrev()->periodicPlanePatch()->shift
                                        : ll->source->periodicPrev()->get_P());
            const VectorDim nextNodePos(ll->next->twin()? ll->next->twin()->  sink->periodicNext()->get_P()
                                        -ll->  sink->periodicPlanePatch()->shift
                                        +ll->next->twin()->  sink->periodicNext()->periodicPlanePatch()->shift
                                        : ll->  sink->periodicNext()->get_P());
            
            
            const VectorDim node0(ll->source->networkNode==parentSegment.source? prevNodePos : nextNodePos);
            const VectorDim node1(ll->source->networkNode==parentSegment.source? ll->source->get_P() : ll->sink->get_P());
            const VectorDim node2(ll->source->networkNode==parentSegment.source? ll->sink->get_P() : ll->source->get_P());
            const VectorDim node3(ll->source->networkNode==parentSegment.source? nextNodePos : prevNodePos);
            
            CatmullRomSplineSegment<dim> cmSeg(node0,node1,node2,node3);
            const double paramU(QuadratureDynamicType::abscissa(parentSegment.quadraturePoints().size(),qID));
            const VectorDim llunitTangent(cmSeg.get_rl(paramU));
            
            pkForce=(stress*parentSegment.burgers()).cross(llunitTangent);
            const MatrixDim totalStress(stress+forceToStress(stackingFaultForce+lineTensionForce,llunitTangent,parentSegment));
            const VectorDim totalForce(pkForce+stackingFaultForce+lineTensionForce);
            glideVelocity=getGlideVelocity(parentSegment,r,totalForce,totalStress,llunitTangent,dL,inclusionID);
        }
        else
        {
            pkForce=(stress*parentSegment.burgers()).cross(rl);
            glideVelocity=getGlideVelocity(parentSegment,r,pkForce+stackingFaultForce+lineTensionForce,stress+forceToStress(stackingFaultForce+lineTensionForce,rl,parentSegment),rl,dL,inclusionID);
        }
    }

    //        template<typename LinkType>
    template<int dim,int corder>
    typename DislocationQuadraturePoint<dim,corder>::MatrixDim DislocationQuadraturePoint<dim,corder>::forceToStress(const VectorDim& force,const VectorDim& unitTangent,const LinkType& parentSegment) const
    {
        // std::cout<<" Curvature norm "<<cmSeg.get_rll(paramU).norm()<<std::endl;
        const double burgersNorm(parentSegment.burgers().norm());
        if(burgersNorm>FLT_EPSILON)
        {
            const double resolvedtensionStress (force.norm()/burgersNorm);
            const VectorDim unitBurgers(parentSegment.burgers()/burgersNorm);
            const MatrixDim tensionStress(resolvedtensionStress*(unitBurgers*parentSegment.glidePlaneNormal().transpose() + parentSegment.glidePlaneNormal()*unitBurgers.transpose()));
            const VectorDim eqPKForce ((tensionStress*parentSegment.burgers()).cross(unitTangent));
            const double eqForceNorm(eqPKForce.norm());
            if (eqForceNorm>FLT_EPSILON)
            {
                if (eqPKForce.dot(force) > 0.0)
                {
                    return tensionStress;
                }
                else
                {
                    return -tensionStress;
                }
            }
        }
        return MatrixDim::Zero();
        
        //            else
        //            {
        //                return MatrixDim::Zero();
        //            }
    }

    template struct DislocationQuadraturePoint<3,0>;





    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::VectorNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityLinearKernel(const int& k) const
    { /*!@param[in] k the current quadrature point
       *\returns The kernel N^T(k)*v(k)*j(k)
       */
        return quadraturePoint(k).SFgaussEx().transpose()*quadraturePoint(k).glideVelocity*quadraturePoint(k).j;
    }


    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::MatrixNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityBilinearKernel(const int& k) const
    { /*! @param[in] k the current quadrature point
       *  The stiffness matrix integrand evaluated at the k-th quadrature point.
       *	\f[
       *		\mathbf{K}^* = \mathbf{N}^T \mathbf{B} \mathbf{N} \frac{dl}{du}
       *	\f]
       */
        const MatrixDimNdof SFEx(quadraturePoint(k).SFgaussEx());
        return SFEx.transpose()*SFEx*quadraturePoint(k).j;
    }


    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::pkKernel(const int& k) const
    { /*!@param[in] k the current quadrature point
       *\returns dF_PK/du=dF_PK/dL*dL/du at quadrature point k, where
       * u in [0,1] is the spline parametrization
       */
        return quadraturePoint(k).pkForce*quadraturePoint(k).j;
    }


    //Added by Yash
    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::MatrixDim DislocationQuadraturePointContainer<dim,corder>::stressKernel(const int& k) const
    { /*!@param[in] k the current quadrature point
       *\returns d_sigma/du=d_sigma/dL*dL/du at quadrature point k, where
       * u in [0,1] is the spline parametrization
       */
        return quadraturePoint(k).stress*quadraturePoint(k).j;
    }

    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::glideVelocityKernel(const int& k) const
    {
        return quadraturePoint(k).glideVelocity*this->quadraturePoint(k).j;
    }


    //template<int dim,int corder>
    //void DislocationQuadraturePointContainer<dim,corder>::updateForcesAndVelocitiesKernelShifted(const LinkType& parentSegment,
    //                                                                                      const StraightDislocationSegment<dim>& ss,
    //                                                                                      const double& L0,
    //                                                                                      const VectorDim& c,
    //                                                                                      const VectorDim& shift,
    //                                                                                      const double& weight)
    //{
    //    SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,
    //                                    parentSegment.source->get_P()+shift,parentSegment.sink->get_P()+shift);
    //
    //    //                        const double dr(ssd.dMin/(L0+ss.length));
    //    const double dr(ssd.dMin/(L0));
    //
    //    if(dr<10.0)
    //    {// full interaction
    //        for (auto& qPoint : quadraturePoints())
    //        {
    //            qPoint.stress += weight*ss.stress(qPoint.r+shift);
    //        }
    //    }
    //    else if(dr<100.0)
    //    {// 2pt interpolation
    //        const MatrixDim stressSource(ss.stress(parentSegment.source->get_P()+shift));
    //        const MatrixDim stressSink(ss.stress(parentSegment.sink->get_P()+shift));
    //        for (auto& qPoint : quadraturePoints())
    //        {
    //            const double u(QuadratureDynamicType::abscissa(this->size(),qPoint.qID));
    //            qPoint.stress += (1.0-u)*weight*stressSource+u*weight*stressSink;
    //        }
    //    }
    //    else
    //    {// 1pt interpolation
    //        const MatrixDim stressC(ss.stress(c+shift));
    //        for (auto& qPoint : quadraturePoints())
    //        {
    //            qPoint.stress += weight*stressC;
    //        }
    //    }
    //
    //    // subtract periodicStressGradient
    ////        if(parentSegment.network().ddBase.simulationParameters.isPeriodicSimulation() && parentSegment.network().ddBase.simulationParameters.correctPeriodicGradients)
    ////        {
    ////            for (auto& qPoint : quadraturePoints())
    ////            {
    //////                const MatrixDim sOld(qPoint.stress);
    ////                const VectorDim dx(qPoint.r-0.5*(parentSegment.network().ddBase.mesh.xMax()+parentSegment.network().ddBase.mesh.xMin()));
    ////                for(int k=0;k<dim;++k)
    ////                {
    ////                    REMOVE THIS ABS
    ////                    qPoint.stress-= parentSegment.network().periodicStressGradient.template block<dim,dim>(k*dim,0)*std::abs(dx(k));
    ////                }
    //////                std::cout<<"\n"<<qPoint.stress-sOld<<std::endl;
    ////            }
    ////        }
    //
    //    if(parentSegment.network().computeElasticEnergyPerLength)
    //    {
    //            for (auto& qPoint : quadraturePoints())
    //            {
    //                qPoint.elasticEnergyPerLength += weight*ss.elasticInteractionEnergy(qPoint.r+shift,qPoint.rl,parentSegment.burgers());
    //            }
    //    }
    //}


//template<int dim,int corder>
//DislocationQuadraturePointContainer<dim,corder>::MatrixDim DislocationQuadraturePointContainer<dim,corder>::stressGradient(const LinkType& parentSegment,
//                                                                                      const StraightDislocationSegment<dim>& ss,
//                                                                                      const VectorDim& x,
//                                                                                      const VectorDim& shift) const
//{
//    MatrixDim temp(MatrixDim::Zero());
//    const VectorDim dx(x-parentSegment.network().ddBase.mesh.xCenter());
//    const Eigen::VectorXd alpha(parentSegment.network().ddBase.periodicLatticeReciprocalBasis.transpose()*dx); // coordinates with resspect to pariodic basis
////    std::cout<<"x="<<x.transpose()<<std::endl;
////    std::cout<<"alpha="<<alpha.transpose()<<std::endl;
//    //    for(size_t a=0;a<1;++a)
//    for(size_t a=0;a<parentSegment.network().ddBase.periodicBasis.size();++a)
//    {
//        Eigen::VectorXd alphaP(alpha);
//        alphaP(a)=+(parentSegment.network().ddBase.simulationParameters.periodicImageSize[a]+0.5);
////        alphaP(a)+=(parentSegment.network().ddBase.simulationParameters.periodicImageSize[a]+0.5);
//
////        alphaP(a)+=0.5;
//
//        const VectorDim xP(parentSegment.network().ddBase.periodicLatticeBasis*alphaP+parentSegment.network().ddBase.mesh.xCenter()+shift);
//        Eigen::VectorXd alphaM(alpha);
//        alphaM(a)=-(parentSegment.network().ddBase.simulationParameters.periodicImageSize[a]+0.5);
////        alphaM(a)-=(parentSegment.network().ddBase.simulationParameters.periodicImageSize[a]+0.5);
//
////        alphaM(a)-=0.5;
//        const VectorDim xM(parentSegment.network().ddBase.periodicLatticeBasis*alphaM+parentSegment.network().ddBase.mesh.xCenter()+shift);
//        temp+=-(ss.stress(xP)-ss.stress(xM))*alpha(a)/(2.0*parentSegment.network().ddBase.simulationParameters.periodicImageSize[a]+1)/4;
////        temp+=(ss.stress(xP)-ss.stress(xM))*std::fabs(alpha(a));
//
//    }
////    std::cout<<temp<<"\n"<<std::endl;
//    return temp;
//}

//template<int dim,int corder>
//double DislocationQuadraturePointContainer<dim,corder>::energyGradient(const LinkType& parentSegment,
//                                                                                      const StraightDislocationSegment<dim>& ss,
//                                                                                      const VectorDim& x,
//                                                                                      const VectorDim& rl,
//                                                                                      const VectorDim& shift) const
//{
//    double temp(0.0);
//    const VectorDim dx(x-parentSegment.network().ddBase.mesh.xCenter());
//    const Eigen::VectorXd alpha(parentSegment.network().ddBase.periodicLatticeReciprocalBasis.transpose()*dx); // coordinates with resspect to pariodic basis
//    for(size_t a=0;a<parentSegment.network().ddBase.periodicBasis.size();++a)
//    {
//        Eigen::VectorXd alphaP(alpha);
//        alphaP(a)=+0.5;
//        const VectorDim xP(parentSegment.network().ddBase.periodicLatticeBasis*alphaP+parentSegment.network().ddBase.mesh.xCenter()+shift);
//        Eigen::VectorXd alphaM(alpha);
//        alphaM(a)=-0.5;
//        const VectorDim xM(parentSegment.network().ddBase.periodicLatticeBasis*alphaM+parentSegment.network().ddBase.mesh.xCenter()+shift);
//        temp+=(ss.elasticInteractionEnergy(xP,rl,parentSegment.burgers())-ss.elasticInteractionEnergy(xM,rl,parentSegment.burgers()))*alpha(a);
//    }
//    return temp;
//}



    template<int dim,int corder>
    void DislocationQuadraturePointContainer<dim,corder>::updateForcesAndVelocitiesKernel(const LinkType& parentSegment,
                                                                                          const StraightDislocationSegment<dim>& ss,
                                                                                          const double& L0,
                                                                                          const VectorDim& c)
    {
        for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
        {
            SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,parentSegment.source->get_P()+shift,parentSegment.sink->get_P()+shift);
            //                        const double dr(ssd.dMin/(L0+ss.length));
            const double dr(ssd.dMin/(L0));

            if(dr<10.0)
            {// full interaction
                for (auto& qPoint : quadraturePoints())
                {
                    qPoint.stress += ss.stress(qPoint.r+shift);
//                    if(   parentSegment.network().ddBase.simulationParameters.isPeriodicSimulation()
//                       && parentSegment.network().ddBase.simulationParameters.correctPeriodicGradients)
//                    {
//                        qPoint.stress -= stressGradient(parentSegment,ss,qPoint.r,shift);
//                    }
                }
            }
            else if(dr<100.0)
            {// 2pt interpolation
                const MatrixDim stressSource(ss.stress(parentSegment.source->get_P()+shift));
                const MatrixDim stressSink(ss.stress(parentSegment.sink->get_P()+shift));
                for (auto& qPoint : quadraturePoints())
                {
                    const double u(QuadratureDynamicType::abscissa(this->size(),qPoint.qID));
                    qPoint.stress += (1.0-u)*stressSource+u*stressSink;
                }
//                if(   parentSegment.network().ddBase.simulationParameters.isPeriodicSimulation()
//                   && parentSegment.network().ddBase.simulationParameters.correctPeriodicGradients)
//                {
//                    const MatrixDim stressGradientSource(stressGradient(parentSegment,ss,parentSegment.source->get_P(),shift));
//                    const MatrixDim stressGradientSink(stressGradient(parentSegment,ss,parentSegment.sink->get_P(),shift));
//                    for (auto& qPoint : quadraturePoints())
//                    {
//                        const double u(QuadratureDynamicType::abscissa(this->size(),qPoint.qID));
//                        qPoint.stress -= ((1.0-u)*stressGradientSource+u*stressGradientSink);
//                    }
//                }
            }
            else
            {// 1pt interpolation
                const MatrixDim stressC(ss.stress(c+shift));
                for (auto& qPoint : quadraturePoints())
                {
                    qPoint.stress += stressC;
                }
//                if(   parentSegment.network().ddBase.simulationParameters.isPeriodicSimulation()
//                   && parentSegment.network().ddBase.simulationParameters.correctPeriodicGradients)
//                {
//                    const MatrixDim stressGradienC(stressGradient(parentSegment,ss,c,shift));
//                    for (auto& qPoint : quadraturePoints())
//                    {
//                        qPoint.stress -= stressGradienC;
//                    }
//                }
            }
                        
            if(parentSegment.network().computeElasticEnergyPerLength)
            {
                for (auto& qPoint : quadraturePoints())
                {
                    qPoint.elasticEnergyPerLength += ss.elasticInteractionEnergy(qPoint.r+shift,qPoint.rl,parentSegment.burgers());
//                    if(   parentSegment.network().ddBase.simulationParameters.isPeriodicSimulation()
//                       && parentSegment.network().ddBase.simulationParameters.correctPeriodicGradients)
//                    {
//                        qPoint.elasticEnergyPerLength -= energyGradient(parentSegment,ss,qPoint.r,qPoint.rl,shift);
//                    }
                }
            }
            
        }
    }

    //template<int dim,int corder>
    //void DislocationQuadraturePointContainer<dim,corder>::updateForcesAndVelocitiesKernel(const LinkType& parentSegment,
    //                                                                                      const StraightDislocationSegment<dim>& ss,
    //                                                                                      const double& L0,
    //                                                                                      const VectorDim& c)
    //{
    //
    //    if(parentSegment.network().ddBase.simulationParameters.periodicImageCentered)
    //    {// Determine which pair of images is the closest
    //
    //
    //        const Eigen::VectorXd alpha(parentSegment.network().ddBase.periodicLatticeReciprocalBasis*);
    //
    //        for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
    //        {
    //        }
    //
    ////        std::map<std::pair<float,float>,VectorDim> shiftMap; // order by distance, and for same distance pick smallest shift vector
    ////        SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,parentSegment.source->get_P(),parentSegment.sink->get_P());
    ////        shiftMap.emplace(std::pair<float,float>(ssd.dMin,0.0),VectorDim::Zero()); // force shift 0 to be default
    ////        for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
    ////        {
    //////            if(shift.squaredNorm()>FLT_EPSILON)
    //////            {
    ////                SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,parentSegment.source->get_P()+shift,parentSegment.sink->get_P()+shift);
    ////                shiftMap.emplace(std::pair<float,float>(ssd.dMin,shift.squaredNorm()),shift);
    //////            }
    ////        }
    ////        const VectorDim nearestShift(shiftMap.begin()->second);
    ////        const double nearestShiftDistance(shiftMap.begin()->first.first);
    ////        const double zeroShiftDistance(SegmentSegmentDistance<dim>(ss.P0,ss.P1,parentSegment.source->get_P(),parentSegment.sink->get_P()).dMin);
    ////        const double nearestShiftWeigth(std::pow(1.0/(nearestShiftDistance+1.0),0.5));
    ////        const double zeroShiftWeigth(std::pow(1.0/(zeroShiftDistance+1.0),0.5));
    ////        const double totalWeigth(nearestShiftWeigth+zeroShiftWeigth);
    ////
    ////
    //////        std::cout<<"nearestShift="<<nearestShift.transpose()<<std::endl;
    //////        std::cout<<"nearestShiftDistance="<<nearestShiftDistance<<std::endl;
    //////        std::cout<<"zeroShiftDistance="<<zeroShiftDistance<<std::endl;
    ////
    ////
    ////        for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
    ////        {
    ////            updateForcesAndVelocitiesKernelShifted(parentSegment,ss,L0,c,shift,zeroShiftWeigth/totalWeigth);
    ////            updateForcesAndVelocitiesKernelShifted(parentSegment,ss,L0,c,shift+nearestShift,nearestShiftWeigth/totalWeigth);
    //////            updateForcesAndVelocitiesKernelShifted(parentSegment,ss,L0,c,shift,0.0);
    //////            updateForcesAndVelocitiesKernelShifted(parentSegment,ss,L0,c,shift+nearestShift,1.0);
    ////
    ////        }
    ////
    //////        std::vector<std::pair<float,VectorDim>> shiftWeigths; // order by distance, and for same distance pick smallest shift vector
    //////        SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,parentSegment.source->get_P(),parentSegment.sink->get_P());
    ////////        shiftMap.emplace(std::pair<float,float>(ssd.dMin,0.0),VectorDim::Zero()); // force shift 0 to be default
    //////
    //////        for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
    //////        {
    //////            SegmentSegmentDistance<dim> ssd(ss.P0,ss.P1,parentSegment.source->get_P()+shift,parentSegment.sink->get_P()+shift);
    //////            shiftWeigths.emplace_back(ssd.dMin,shift);
    //////            //std::pair<float,float>(,shift.squaredNorm()),shift);
    //////        }
    //////
    //////        double totWeight;
    //////        for(const auto& baseShift : shiftWeigths)
    //////        {
    //////            totWeight+=baseShift.first;
    //////        }
    //////
    //////        for(const auto& baseShift : shiftWeigths)
    //////        {
    //////            for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
    //////            {
    //////                updateForcesAndVelocitiesKernelShifted(parentSegment,ss,L0,c,shift+baseShift.second,(totWeight-baseShift.first)/totWeight);
    //////            }
    //////        }
    //    }
    //    else
    //    {
    //        for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
    //        {
    //            updateForcesAndVelocitiesKernelShifted(parentSegment,ss,L0,c,shift,1.0);
    //        }
    //    }
    //
    ////    for(const auto& baseShift : parentSegment.network().ddBase.periodicShifts)
    //
    //
    ////    if(parentSegment.network().computeElasticEnergyPerLength)
    ////    {
    ////        for(const auto& shift : parentSegment.network().periodicShifts)
    ////        {
    ////            for (auto& qPoint : quadraturePoints())
    ////            {
    ////                qPoint.elasticEnergyPerLength += ss.elasticInteractionEnergy(qPoint.r+shift,qPoint.rl,parentSegment.burgers());
    ////            }
    ////        }
    ////    }
    //}

    template<int dim,int corder>
    void DislocationQuadraturePointContainer<dim,corder>::updateForcesAndVelocities(const LinkType& parentSegment,
                                                                                    const double&,
                                                                                    const double& )
    {
        // updateQuadraturePoints(parentSegment,quadPerLength,isClimbStep);
        
        
        const bool useLoopBasedKernels(false);
        
        if(this->size())
        {
            if(parentSegment.network().computeDDinteractions)
            {
                const double L0(parentSegment.chord().norm());
                const VectorDim c(0.5*(parentSegment.source->get_P()+parentSegment.sink->get_P()));
                if(useLoopBasedKernels)
                {
                    for(const auto& loop : parentSegment.network().loops())
                    {
                        for(const auto& patch : loop.second.lock()->patches().globalPatches())
                        {
                            for(size_t k0=0;k0<patch.second.size();++k0)
                            {
                                const size_t k1(k0+1<patch.second.size()? k0+1 : 0);
                                const VectorDim& P0(patch.second[k0]);
                                const VectorDim& P1(patch.second[k1]);
                                const VectorDim chord(P1-P0);
                                const double chordLength(chord.norm());
                                if(chordLength>FLT_EPSILON)
                                {
                                    const VectorDim tangent(chord/chordLength);
                                    //                                const StraightDislocationSegment<dim> ss(parentSegment.network().poly,P0-patch.first->shift,P1-patch.first->shift,loop.second.lock()->burgers(),chordLength,tangent);
                                    const StraightDislocationSegment<dim> ss(parentSegment.network().ddBase.poly,P0,P1,loop.second.lock()->burgers(),chordLength,tangent,parentSegment.network().ddBase.EwaldLength);
                                    updateForcesAndVelocitiesKernel(parentSegment,ss,L0,c);
                                }
                            }
                        }
                    }
                }
                else
                {
                    for(const auto& link : parentSegment.network().networkLinks())
                    {
                        if(   !link.second.lock()->hasZeroBurgers()
                           && !(link.second.lock()->isBoundarySegment() && parentSegment.network().ddBase.simulationParameters.simulationType==DDtraitsIO::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
                           //                           && !(link.second.lock()->isVirtualBoundarySegment() && parentSegment.network().simulationParameters.simulationType==DDtraitsIO::PERIODIC_IMAGES)
                           )
                        {
                            
                            const StraightDislocationSegment<dim>& ss(link.second.lock()->straight);
                            updateForcesAndVelocitiesKernel(parentSegment,ss,L0,c);
                        }
                    }
                }
            }
            
            
            // Stacking fault contribution in the matrix
            if(parentSegment.slipSystem() && parentSegment.glidePlanes().size()==1)
            {
                const auto& glidePlane(**parentSegment.glidePlanes().begin());
                const auto& slipSystem(*parentSegment.slipSystem());
                
                
                // Add stacking fault force
                const double eps=1.0e-2;
                VectorDim outDir(parentSegment.unitDirection().cross(slipSystem.unitNormal));
                const double outDirNorm(outDir.norm());
                if(outDirNorm>FLT_EPSILON)
                {
                    outDir/=outDirNorm;
                    std::vector<std::pair<VectorDim,VectorDim>> qPointSlip(quadraturePoints().size(),std::make_pair(VectorDim::Zero(),VectorDim::Zero())); // accumulated slip vectors (outside,inside) for each qPoint
                    
                    for(const auto& weakSourceLoop : parentSegment.network().loops())
                    {
                        const auto sourceLoop(weakSourceLoop.second.lock());
                        if(sourceLoop->slipSystem())
                        {
                            if(slipSystem.n==sourceLoop->slipSystem()->n) // "opposite" slip systems have same normal and opposite slip directions, so this is ok
                            {// same glide plane family
                                
                                for(auto& qPoint : quadraturePoints())
                                {
                                    if(qPoint.inclusionID>=0 || sourceLoop->slipSystem()->isPartial())
                                    {// point inside inclusion, or contribution of a partial slip system
                                        qPoint.slipVectors.first -=sourceLoop->patches().windingNumber(qPoint.r + eps*outDir)*sourceLoop->burgers(); // slip vector is negative the burgers vector. windingNumber will skip parallel planes
                                        qPoint.slipVectors.second-=sourceLoop->patches().windingNumber(qPoint.r - eps*outDir)*sourceLoop->burgers(); // slip vector is negative the burgers vector. windingNumber will skip parallel planes
                                    }
                                }
                            }
                        }
                    }
                    
                    for(auto& qPoint : quadraturePoints())
                    {
                        if((qPoint.slipVectors.first-qPoint.slipVectors.second).squaredNorm()>FLT_EPSILON)
                        {
                            //                                auto& qPoint(loopLink->networkLink()->quadraturePoints()[q]);
                            if(qPoint.inclusionID<0)
                            {// qPoint is not inside an inclusion, we use the matrix gamma-surface
                                
                                const double gamma1(slipSystem.n.misfitEnergy(qPoint.slipVectors.first));  // outer point
                                const double gamma2(slipSystem.n.misfitEnergy(qPoint.slipVectors.second)); // inner point
                                
                                double gammaNoise(0.0);
                                if(slipSystem.planeNoise)
                                {
                                    gammaNoise=std::get<2>(slipSystem.gridInterp(qPoint.r-glidePlane.P));
                                }
                                qPoint.stackingFaultForce+= -(gamma2-gamma1+gammaNoise)*outDir;
                            }
                            else
                            {// qPoint is inside an inclusion, we use the inclusion gamma-surface
                                
                                const auto& secondPhase(parentSegment.network().eshelbyInclusions().at(qPoint.inclusionID)->secondPhase);
                                if(secondPhase)
                                {
                                    const double gamma1(secondPhase->misfitEnergy(qPoint.slipVectors.first ,&slipSystem.n));  // outer point
                                    const double gamma2(secondPhase->misfitEnergy(qPoint.slipVectors.second,&slipSystem.n)); // inner point
                                    qPoint.stackingFaultForce+= -(gamma2-gamma1)*outDir;
                                }
                            }
                        }
                    }
                }
                
                
                
                // Add solid-soution noise
                if(slipSystem.planeNoise)
                {
                    for (auto& qPoint : quadraturePoints())
                    {
                        const auto noiseVal(slipSystem.gridInterp(qPoint.r-glidePlane.P));
                        qPoint.stress += std::get<0>(noiseVal);
                    }
                }
            }
            
            // Add other stress contributions, and compute PK force
            for (auto& qPoint : quadraturePoints())
            {
                
                if(parentSegment.network().externalLoadController)
                {// Add stress of externalLoadController
                    qPoint.stress += parentSegment.network().externalLoadController->stress(qPoint.r);
                }
                
                if(parentSegment.network().bvpSolver)
                {// Add BVP stress
                    qPoint.stress += parentSegment.network().bvpSolver->stress(qPoint.r,parentSegment.source->includingSimplex());
                }
                
                //                    for(const auto& sStraight : parentSegment.network().poly.grainBoundaryDislocations() )
                //                    {// Add GB stress
                //                        qPoint.stress += sStraight.stress(qPoint.r);
                //                    }
                
                for(const auto& inclusion : parentSegment.network().eshelbyInclusions() )
                {// Add EshelbyInclusions stress
                    for(const auto& shift : parentSegment.network().ddBase.periodicShifts)
                    {
                        qPoint.stress += inclusion.second->stress(qPoint.r+shift);
                    }
                }
                
                //Add line tension contribution
                if (parentSegment.network().alphaLineTension>0.0 && !parentSegment.hasZeroBurgers() && parentSegment.chordLength()>FLT_EPSILON)
                {
                    //Add the line tension contribution due to only non -zero segments
                    for (const auto &ll : parentSegment.loopLinks())
                    {
                        
                        //                    assert(ll->source->periodicPrev()!=nullptr && "PeriodicPrev must exist from network link for line tension contribution");
                        //                    assert(ll->sink->periodicNext()!=nullptr && "PeriodicNext must exist from network link for line tension contribution");
                        //
                        //                    const VectorDim prevNodePos(ll->prev->twin()? ll->prev->twin()->source->periodicPrev()->get_P()
                        //                                                -ll->source->periodicPlanePatch()->shift
                        //                                                +ll->prev->twin()->source->periodicPrev()->periodicPlanePatch()->shift
                        //                                                : ll->source->periodicPrev()->get_P());
                        //                    const VectorDim nextNodePos(ll->next->twin()? ll->next->twin()->  sink->periodicNext()->get_P()
                        //                                                -ll->  sink->periodicPlanePatch()->shift
                        //                                                +ll->next->twin()->  sink->periodicNext()->periodicPlanePatch()->shift
                        //                                                : ll->  sink->periodicNext()->get_P());
                        //                    CatmullRomSplineSegment<dim> spline(prevNodePos,ll->source->get_P(),ll->sink->get_P(),nextNodePos);
                        //                    const double paramUTemp ((qPoint.r-parentSegment.source->get_P()).norm()/parentSegment.chordLength());
                        const double paramUTemp(QuadratureDynamicType::abscissa(parentSegment.quadraturePoints().size(),qPoint.qID));
                        const double paramU (ll->source->networkNode==parentSegment.source ? paramUTemp : 1.0-paramUTemp);
                        const auto spline(ll->spline());
                        const VectorDim llunitTangent(spline.get_rl(paramU));
                        const double qPointEnergyDensity (parentSegment.network().alphaLineTension * parentSegment.network().ddBase.poly.C2 * (ll->loop->burgers().squaredNorm()-parentSegment.network().ddBase.poly.nu*(std::pow(llunitTangent.dot(ll->loop->burgers()),2))));
                        qPoint.coreEnergyPerLength+=qPointEnergyDensity;
                        
                        const VectorDim rll(spline.get_rll(paramU));
                        if(!std::isnan(rll.squaredNorm()))
                        {// curvature can be nan for straight lines
                            const VectorDim qPointForceVector (qPointEnergyDensity * rll);
                            qPoint.lineTensionForce+=qPointForceVector;
                        }
                    }
                }
                
                qPoint.updateForcesAndVelocities(parentSegment);
            }
            
        }
    }



    //        template<typename LinkType>
    template<int dim,int corder>
    void DislocationQuadraturePointContainer<dim,corder>::updateQuadraturePoints(const LinkType& parentSegment,
                                                                                 const double& quadPerLength,
                                                                                 const bool& isClimbStep)
    {
        
        this->clear();
        
        if(    !parentSegment.hasZeroBurgers()
           &&  !parentSegment.isBoundarySegment()
           &&  (!parentSegment.isSessile() || isClimbStep || parentSegment.network().computeElasticEnergyPerLength)
           //       &&  !parentSegment.isVirtualBoundarySegment()
           )
        {
            const int order=QuadPowDynamicType::lowerOrder(quadPerLength*parentSegment.chord().norm());
            const MatrixNcoeff  SFCH(parentSegment.sfCoeffs());
            const MatrixNcoeffDim qH(parentSegment.hermiteDofs());
            for(int q=0;q<order;++q)
            {
                this->emplace_back(parentSegment,q,order,SFCH,qH);
            }
        }
    }


    template<int dim,int corder>
    const typename DislocationQuadraturePointContainer<dim,corder>::DislocationQuadraturePointContainerType& DislocationQuadraturePointContainer<dim,corder>::quadraturePoints() const
    {
        return *this;
    }


    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::DislocationQuadraturePointContainerType& DislocationQuadraturePointContainer<dim,corder>::quadraturePoints()
    {
        return *this;
    }


    template<int dim,int corder>
    const typename DislocationQuadraturePointContainer<dim,corder>::DislocationQuadraturePointType& DislocationQuadraturePointContainer<dim,corder>::quadraturePoint(const int& k) const
    {
        return this->operator[](k);
    }




    // VectorNdof nodalVelocityVector() const
    // { /*\returns The segment-integrated nodal velocity vector int N^T*v dL
    //    */
    //     VectorNdof Fq(VectorNdof::Zero());
    //     QuadratureDynamicType::integrate(this->size(),this,Fq,&DislocationQuadraturePointContainerType::nodalVelocityLinearKernel);
    //     return Fq;
    // }

    //        template<typename LinkType>
    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::VectorNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityVector(const LinkType& parentSegment) const
    { /*\returns The segment-integrated nodal velocity vector int N^T*v dL
       */
        VectorNdof Fq(VectorNdof::Zero());
        if(parentSegment.chordLength()>FLT_EPSILON)
        {
            QuadratureDynamicType::integrate(this->size(),this,Fq,&DislocationQuadraturePointContainerType::nodalVelocityLinearKernel);
        }
        return Fq;
    }


    //        template<typename LinkType>
    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::MatrixNdof DislocationQuadraturePointContainer<dim,corder>::nodalVelocityMatrix(const LinkType& parentSegment) const
    {
        MatrixNdof Kqq(MatrixNdof::Zero());
        if(corder==0)
        {// Analytical integral can be performed
            const double L(parentSegment.chordLength());
            Kqq<<L/3.0,  0.0,  0.0, L/6.0,  0.0,  0.0,
            /**/   0.0,L/3.0,  0.0,   0.0,L/6.0,  0.0,
            /**/   0.0,  0.0,L/3.0,   0.0,  0.0,L/6.0,
            /**/ L/6.0,  0.0,  0.0, L/3.0,  0.0,  0.0,
            /**/   0.0,L/6.0,  0.0,   0.0,L/3.0,  0.0,
            /**/   0.0,  0.0,L/6.0,   0.0,  0.0,L/3.0;
        }
        else
        {// Numerical integral must be performed
            assert(0 && "WE NEED TO INCREASE qOrder FOR THE FOLLOWING INTEGRATION, SINCE EVEN FOR LINEAR SEGMENTS Kqq IS NOT INTEGRATED CORRECLTY FOR SMALL qOrder");
            assert(0 && "THIS MUST RETURN A NON_ZERO VALUE EVEN FOR SEGMENTS WHICH DONT HAAVE QUADRATURE POINTS TO CORRECTLY SOVE K*V=F");
            QuadratureDynamicType::integrate(this->size(),this,Kqq,&DislocationQuadraturePointContainerType::nodalVelocityBilinearKernel);
        }
        return Kqq;
    }

    /*************************************************************/
    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::glideVelocityIntegral() const
    {/*!\returns The integral of the glide velocity over the segment.
      */
        VectorDim V(VectorDim::Zero());
        QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,V,&DislocationQuadraturePointContainerType::glideVelocityKernel);
        return V;
    }


    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::VectorDim DislocationQuadraturePointContainer<dim,corder>::pkIntegral() const
    {/*!\returns The integral of the PK force over the segment.
      */
        VectorDim F(VectorDim::Zero());
        QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,F,&DislocationQuadraturePointContainerType::pkKernel);
        return F;
    }


    //Added by Yash
    template<int dim,int corder>
    typename DislocationQuadraturePointContainer<dim,corder>::MatrixDim DislocationQuadraturePointContainer<dim,corder>::stressIntegral() const
    {/*!\returns The integral of the stress over the segment.
      */
        MatrixDim stressInt(MatrixDim::Zero());
        QuadratureDynamicType::integrate(this->quadraturePoints().size(),this,stressInt,&DislocationQuadraturePointContainerType::stressKernel);
        return stressInt;
    }

    template class DislocationQuadraturePointContainer<3,0>;


}
#endif

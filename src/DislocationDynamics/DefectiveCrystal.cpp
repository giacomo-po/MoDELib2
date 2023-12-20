/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_cpp_
#define model_DefectiveCrystal_cpp_

#include <DefectiveCrystal.h>

namespace model
{
                
        template <int _dim, short unsigned int corder>
        std::unique_ptr<ExternalLoadControllerBase<_dim>> DefectiveCrystal<_dim,corder>::getExternalLoadController(const DislocationDynamicsBase<dim>& ddBase,const MatrixDim& plasticStrain_in)
        {
                     
            if(ddBase.simulationParameters.simulationType==DDtraitsIO::FINITE_FEM)
            {
                return std::unique_ptr<ExternalLoadControllerBase<_dim>>(nullptr);
            }
            else
            {
                if(ddBase.simulationParameters.externalLoadControllerName=="UniformExternalLoadController")
                {
                    return std::unique_ptr<ExternalLoadControllerBase<_dim>>(new UniformExternalLoadController<_dim>(ddBase,plasticStrain_in));
                }
                else
                {
                    std::cout<<"Unknown externalLoadController name "<<ddBase.simulationParameters.externalLoadControllerName<<"No controller applied."<<std::endl;
                    return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
                }
            }
        }
        
        template <int _dim, short unsigned int corder>
        void DefectiveCrystal<_dim,corder>::updateLoadControllers(const long int& runID, const bool& isClimbStep)
        {/*! Updates bvpSolver using the stress and displacement fields of the
          *  current DD configuration.
          */
            const int quadraturePerTriangle=37;
            if(bvpSolver)
            {
                if (!(runID%bvpSolver->stepsBetweenBVPupdates))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    std::cout<<"Updating bvpSolver ... "<<std::endl;
                    bvpSolver->template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*DN, isClimbStep);
                }
            }
            if (externalLoadController)
            {
                std::cout<<"Updating externalLoadController... "<<std::endl;
                externalLoadController->update(plasticStrain());
            }
        }
        
        template <int _dim, short unsigned int corder>
        DefectiveCrystal<_dim,corder>::DefectiveCrystal(DislocationDynamicsBase<_dim>& ddBase_in) :
        /* init */ ddBase(ddBase_in)
//        /* init */,DN(ddBase.simulationParameters.useDislocations? new DislocationNetworkType(ddBase.simulationParameters,ddBase.mesh,ddBase.poly,bvpSolver,externalLoadController,ddBase.periodicShifts,ddBase.simulationParameters.runID) : nullptr)
        /* init */,DN(ddBase.simulationParameters.useDislocations? new DislocationNetworkType(ddBase,bvpSolver,externalLoadController) : nullptr)
        /* init */,CS(ddBase.simulationParameters.useCracks? new CrackSystemType() : nullptr)
        /* init */,bvpSolver(ddBase.simulationParameters.simulationType==DDtraitsIO::FINITE_FEM? new BVPsolverType(ddBase.mesh,*DN) : nullptr)
        /* init */,externalLoadController(getExternalLoadController(ddBase,plasticStrain()))
        {

        }

        template <int _dim, short unsigned int corder>
        double DefectiveCrystal<_dim,corder>::getMaxVelocity() const
        {
            double vmax = 0.0;

            for (const auto &nodeIter : DN->networkNodes())
            {
                    const double vNorm(nodeIter.second.lock()->get_V().norm());
                    if (vNorm > vmax)
                    {
                        vmax = vNorm;
                    }
            }
            return vmax;
        }
        
        template <int _dim, short unsigned int corder>
        void DefectiveCrystal<_dim,corder>::singleGlideStep()
        {
            std::cout<<blueBoldColor<< "runID="<<ddBase.simulationParameters.runID<<" (of "<<ddBase.simulationParameters.Nsteps<<")"
            /*                    */<< ", time="<<ddBase.simulationParameters.totalTime;
            if(DN)
            {
                std::cout<< ": networkNodes="<<DN->networkNodes().size()
                /*                    */<< ", networkSegments="<<DN->networkLinks().size()
                /*                    */<< ", loopNodes="<<DN->loopNodes().size()
                /*                    */<< ", loopSegments="<<DN->loopLinks().size()
                /*                    */<< ", loops="<<DN->loops().size();
            }
            std::cout<< defaultColor<<std::endl;

            if(DN)
            {
                DislocationNode<dim,corder>::totalCappedNodes=0;
                DN->updateGeometry();
                updateLoadControllers(ddBase.simulationParameters.runID, false);
                const double maxVelocity(getMaxVelocity());
                DN->assembleGlide(ddBase.simulationParameters.runID, maxVelocity);
                DN->storeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
                DN->solveGlide();
                ddBase.simulationParameters.dt=DN->timeIntegrator.getGlideTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem
                DN->updateRates();
                DN->io().output(ddBase.simulationParameters.runID);
                DN->moveGlide(ddBase.simulationParameters.dt);
                DN->executeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
                if (DN->capMaxVelocity)
                {
                    std::cout<<redBoldColor<<"( "<<(DislocationNode<dim,corder>::totalCappedNodes)<<" total nodes capped "<<defaultColor<<std::endl;
                    std::cout<<redBoldColor<<", "<<(double(DislocationNode<dim,corder>::totalCappedNodes)/double(DN->networkNodes().size()))<<" fraction of nodes capped "
                    <<defaultColor<<" )"<<std::endl;
                }
            }
            ddBase.simulationParameters.totalTime+=ddBase.simulationParameters.dt;
            ++ddBase.simulationParameters.runID;
        }

        template <int _dim, short unsigned int corder>
        void DefectiveCrystal<_dim,corder>::runGlideSteps()
        {/*! Runs a number of simulation time steps defined by simulationParameters.Nsteps
          */
            const auto t0= std::chrono::system_clock::now();
            while (ddBase.simulationParameters.runID<ddBase.simulationParameters.Nsteps)
            {
                std::cout<<std::endl; // leave a blank line
                singleGlideStep();
            }
            if(DN)
            {
                DN->updateGeometry();
            }
            std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<ddBase.simulationParameters.Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }

        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticDistortion() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortion();
            }
            if(CS)
            {
                temp+=CS->plasticDistortion();
            }
            return temp;
        }

        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticStrain() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(plasticDistortion());
            return 0.5*(temp+temp.transpose());
        }

        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticDistortionRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortionRate();
            }
            if(CS)
            {
                temp+=CS->plasticDistortionRate();
            }
            return temp;
        }

        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticStrainRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(plasticDistortionRate());
            return 0.5*(temp+temp.transpose());
        }

    template class DefectiveCrystal <3,0>;
}
#endif

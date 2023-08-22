/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UniformExternalLoadController_cpp_
#define model_UniformExternalLoadController_cpp_


#include <UniformExternalLoadController.h>



namespace model
{

    template <int dim>
    UniformExternalLoadController<dim>::UniformExternalLoadController(const DislocationDynamicsBase<dim>& ddBase,const MatrixDim& plasticStrain_in) :
    /* init */ ExternalLoadControllerBase<dim>(ddBase,ddBase.simulationParameters.traitsIO.simulationFolder+"/inputFiles/uniformExternalLoadController.txt")
    /* init */,voigtTraits(getVoigtTraits(ddBase.simulationParameters.traitsIO.simulationFolder+"/inputFiles/uniformExternalLoadController.txt"))
    /* init */,enable(TextFileParser(this->inputFileName).readScalar<int>("enable",true))
    /* init */,relaxSteps(TextFileParser(this->inputFileName).readScalar<int>("relaxSteps",true))
    {
        std::cout<<greenColor<<"Initializing UniformExternalLoadController at runID="<<this->ddBase.simulationParameters.runID<<defaultColor<<std::endl;
                
        // initilize the default voigt machinestiffness
        this->MachineStiffnessRatio.block(0,0,1,dim)=this->MachineStiffnessRatio.block(0,0,1,dim)*(2+this->lambda); // using strategy 2 in test.m file  alpha*diag(C)  (first three)*(2*mu+this->lambda)  (first three)*mu
        Eigen::Matrix<double,voigtSize,voigtSize>  Cinv=Eigen::Matrix<double,voigtSize,voigtSize>::Identity();
        Cinv.block(0,0,dim,dim)<<this->nu_use/this->nu, -this->nu_use,       -this->nu_use,
        -this->nu_use,   this->nu_use/this->nu,     -this->nu_use,
        -this->nu_use,   -this->nu_use,       this->nu_use/this->nu;
        Eigen::Matrix<double,voigtSize,voigtSize>  machinestiffness=this->MachineStiffnessRatio.asDiagonal();;
        this->stressmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse();
        this->strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse()*machinestiffness;
        
        update(plasticStrain_in);
    }

    template <int dim>
    SymmetricVoigtTraits<dim> UniformExternalLoadController<dim>::getVoigtTraits(const std::string& )
    {
        Eigen::Matrix<size_t,SymmetricVoigtTraits<dim>::voigtSize,2>     voigtorder;
        size_t vnum=0;
        for (size_t i=0;i<dim;i++)
        {
            for (size_t j=i;j<dim;j++)
            {
                voigtorder.row(vnum)<<j-i,j;
                vnum++;
            }
        }
        std::cout<<" voigtorder="<<voigtorder<<std::endl;
        return SymmetricVoigtTraits<dim>(voigtorder);
    }


    template <int dim>
    typename UniformExternalLoadController<dim>::MatrixDim UniformExternalLoadController<dim>::stressconsidermachinestiffness(const MatrixDim& S_strain,const MatrixDim& S_stress) const
    {/*!For stress updating,
      * \f[
      * \dot{\sigma}=\frac{\alpha}{1+\alpha} C::(\dot{\epsilon}-\dot{\epsilon}^p)+ \frac{1}{1+\alpha} \dot{\sigma}
      * \f]
      * \alpha=this->MachineStiffnessRatio, which is the stiffness ratio between machine stiffness and sample stiffness.
      * More details see <Physical Review Letters 117, 155502, 2016>.
      */
        return voigtTraits.v2m(this->strainmultimachinestiffness*voigtTraits.m2v(S_strain,true)+this->stressmultimachinestiffness*voigtTraits.m2v(S_stress,false),false);
        
    }

    template <int dim>
    typename UniformExternalLoadController<dim>::MatrixDim UniformExternalLoadController<dim>::straininducedStress(const MatrixDim& dstrain, const double& lambda) const
    {/*!for isotropic linear case,
      * \f[
      * \sigma_{ij}=\this->lambda \epsilon_{kk} \delta_{ij}+ \mu (\epsilon_{ij}+\epsilon_{ji})
      * \f]
      */
        return this->lambda*dstrain.trace()*MatrixDim::Identity()+2.0*dstrain;
    }

    template <int dim>
    typename UniformExternalLoadController<dim>::MatrixDim UniformExternalLoadController<dim>::elasticstrain(const MatrixDim& _externalStress,const double& nu_use) const
    {/*!for isotropic linear case,
      * \f[
      * \epsilon_{ij}=\frac{1}{E}((1+\nu)\sigma_{ij}-\nu \sigma_{kk}\delta_{ij})
      * \f]
      */
        return -nu_use*_externalStress.trace()*MatrixDim::Identity()+_externalStress*0.5;
    }

    template <int dim>
    typename UniformExternalLoadController<dim>::MatrixDim UniformExternalLoadController<dim>::stress(const VectorDim&) const
    {
        return enable? this->ExternalStress : MatrixDim::Zero();
    }

    template <int dim>
    typename UniformExternalLoadController<dim>::MatrixDim UniformExternalLoadController<dim>::strain(const VectorDim&) const
    {
        return enable? this->ExternalStrain : MatrixDim::Zero();
    }

    template <int dim>
    void UniformExternalLoadController<dim>::update(const MatrixDim& plasticStrain_in)
    {
        if(this->ddBase.simulationParameters.runID>=relaxSteps)
        {
            this->plasticStrain=plasticStrain_in;
            this->ExternalStress=stressconsidermachinestiffness(this->ExternalStrain0+this->ExternalStrainRate*this->ddBase.simulationParameters.totalTime-this->plasticStrain,this->ExternalStress0+this->ExternalStressRate*this->ddBase.simulationParameters.totalTime);  //2017-12-7
            this->ExternalStrain=elasticstrain(this->ExternalStress,this->nu_use)+this->plasticStrain;
        }
    }

    template <int dim>
    void UniformExternalLoadController<dim>::output(std::ofstream& f_file,
                                                    std::ofstream& F_labels) const
    {
        f_file<<this->ExternalStrain.row(0)<<" "<<this->ExternalStrain.row(1)<<" "<<this->ExternalStrain.row(2)<<" "<<this->ExternalStress.row(0)<<" "<<this->ExternalStress.row(1)<<" "<<this->ExternalStress.row(2)<<" ";
        
        if(this->ddBase.simulationParameters.runID==0)
        {
            F_labels<<"e_11\n";
            F_labels<<"e_12\n";
            F_labels<<"e_13\n";
            F_labels<<"e_21\n";
            F_labels<<"e_22\n";
            F_labels<<"e_23\n";
            F_labels<<"e_31\n";
            F_labels<<"e_32\n";
            F_labels<<"e_33\n";
            F_labels<<"s_11 [mu]\n";
            F_labels<<"s_12 [mu]\n";
            F_labels<<"s_13 [mu]\n";
            F_labels<<"s_21 [mu]\n";
            F_labels<<"s_22 [mu]\n";
            F_labels<<"s_23 [mu]\n";
            F_labels<<"s_31 [mu]\n";
            F_labels<<"s_32 [mu]\n";
            F_labels<<"s_33 [mu]"<<std::endl;
        }
    }


    template class UniformExternalLoadController<3>;

}
#endif

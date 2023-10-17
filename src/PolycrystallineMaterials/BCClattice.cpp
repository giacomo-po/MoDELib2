/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_BCClattice_cpp_
#define model_BCClattice_cpp_

#include <BCClattice.h>
#include <DislocationMobilityBCC.h>

namespace model
{
    
        
        BCClattice<3>::BCClattice(const MatrixDim& Q,const PolycrystallineMaterialBase& material,const std::string& polyFile) :
        /* init */ SingleCrystalBase<dim>(getLatticeBasis(),Q)
        /* init */,PlaneNormalContainerType(getPlaneNormals(material,polyFile))
        /* init */,SlipSystemContainerType(getSlipSystems(material,polyFile,*this))
        /* init */,SecondPhaseContainerType(getSecondPhases(material,*this))
        {
            
        }
        
        Eigen::Matrix<double,3,3> BCClattice<3>::getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << -1.0,  1.0,  1.0,
            /*   */  1.0, -1.0,  1.0,
            /*   */  1.0,  1.0, -1.0;
            
            return temp/sqrt(3.0);
        }

    const typename BCClattice<3>::PlaneNormalContainerType& BCClattice<3>::planeNormals() const
    {
        return *this;
    }

    const typename BCClattice<3>::SlipSystemContainerType& BCClattice<3>::slipSystems() const
    {
        return *this;
    }

    const typename BCClattice<3>::SecondPhaseContainerType& BCClattice<3>::secondPhases() const
    {
        return *this;
    }


//    const typename BCClattice<3>::DislocationMobilityContainerType& dislocationMobilities() const
//    {
//        return *this;
//    }


        
        std::vector<std::shared_ptr<GlidePlaneBase>> BCClattice<3>::getPlaneNormals(const PolycrystallineMaterialBase& ,
                                                                                    const std::string& ) const
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the BCC lattice
          */
            
//            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
//            typedef LatticeVector<dim> LatticeVectorType;

            LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),*this);
            LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),*this);
            LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),*this);
            LatticeVectorType  y((VectorDimI()<<1,1,1).finished(),*this);

            std::vector<std::shared_ptr<GlidePlaneBase>> temp;
            temp.emplace_back(new GlidePlaneBase(a3,a1,nullptr));
            temp.emplace_back(new GlidePlaneBase( y,a2,nullptr));
            temp.emplace_back(new GlidePlaneBase(a2,a3,nullptr));
            temp.emplace_back(new GlidePlaneBase( y,a1,nullptr));
            temp.emplace_back(new GlidePlaneBase(a1,a2,nullptr));
            temp.emplace_back(new GlidePlaneBase( y,a3,nullptr));
            
            return temp;

        }
        
        std::vector<std::shared_ptr<SlipSystem>> BCClattice<3>::getSlipSystems(const PolycrystallineMaterialBase& material,
                                                                               const std::string& polyFile,
                                                                               const PlaneNormalContainerType& plN)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            
            const std::string dislocationMobilityType(TextFileParser(polyFile).readString("dislocationMobilityType",true));
            DislocationMobilitySelector mobilitySelector("BCC");
            const std::shared_ptr<DislocationMobilityBase> mobility110(mobilitySelector.getMobility(dislocationMobilityType,material));

            
//            const std::shared_ptr<DislocationMobilityBase> mobility110(new DislocationMobilityBCC(material));

            
            const int solidSolutionNoiseMode(TextFileParser(polyFile).readScalar<int>("solidSolutionNoiseMode",true));
            const int stackingFaultNoiseMode(TextFileParser(polyFile).readScalar<int>("stackingFaultNoiseMode",true));
            std::shared_ptr<GlidePlaneNoise> planeNoise((solidSolutionNoiseMode||stackingFaultNoiseMode)? new GlidePlaneNoise(polyFile,material) : nullptr);

            
//            typedef Eigen::Matrix<double,dim,1> VectorDimD;
            const double d110(this->reciprocalLatticeDirection(this->C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());

            std::vector<std::shared_ptr<SlipSystem>> temp;
            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase->planeSpacing()-d110)<FLT_EPSILON)
                {// a {110} plane
                    const auto& a1(planeBase->primitiveVectors.first);
                    const auto& a3(planeBase->primitiveVectors.second);
                    
                    std::vector<RationalLatticeDirection<3>> slipDirs;
                    
                    if(material.enabledSlipSystems.find("full")!=material.enabledSlipSystems.end())
                    {
                        slipDirs.emplace_back(Rational( 1,1),a1);
                        slipDirs.emplace_back(Rational(-1,1),a1);
                        slipDirs.emplace_back(Rational( 1,1),a3);
                        slipDirs.emplace_back(Rational(-1,1),a3);
                    }
                    
//                    temp.emplace_back(new SlipSystem(*planeBase, a1,mobility110,nullptr,planeNoise));
//                    temp.emplace_back(new SlipSystem(*planeBase,a1*(-1),mobility110,nullptr,planeNoise));
//                    temp.emplace_back(new SlipSystem(*planeBase, a3,mobility110,nullptr,planeNoise));
//                    temp.emplace_back(new SlipSystem(*planeBase,a3*(-1),mobility110,nullptr,planeNoise));
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace_back(new SlipSystem(*planeBase, slipDir,mobility110,planeNoise));
                    }

                }
            }

            return temp;
        }
        
        

    typename BCClattice<3>::SecondPhaseContainerType BCClattice<3>::getSecondPhases(const PolycrystallineMaterialBase& material,
                                                                                const PlaneNormalContainerType& planeNormals)
    {
        
//        std::vector<std::shared_ptr<SecondPhase<3>>> temp;
        SecondPhaseContainerType temp;

        for(const std::string& sp : material.enabledSecondPhases)
        {
            if(sp=="chi")
            {
                Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors110(TextFileParser(material.materialFile).readMatrixCols<double>("chiWaveVectors110",2,true));
                Eigen::Matrix<double,Eigen::Dynamic,3> f110(TextFileParser(material.materialFile).readMatrixCols<double>("chiGammaSurfacePoints110",3,true));
                const int rotSymm110(1);
                std::vector<Eigen::Matrix<double,2,1>> mirSymm110;
                mirSymm110.push_back((Eigen::Matrix<double,2,1>()<<1.0,0.0).finished()); // symm with respect to local y-axis
                mirSymm110.push_back((Eigen::Matrix<double,2,1>()<<0.0,1.0).finished()); // symm with respect to local x-axis
                const Eigen::Matrix<double,2,2> A110((Eigen::Matrix<double,2,2>()<< 1.0,-0.5,
                                                                                 0.0,0.5*std::sqrt(3.0)).finished());
                std::shared_ptr<GammaSurface> gammaSurface110(new GammaSurface(A110,waveVectors110,f110,rotSymm110,mirSymm110));
                const double d110(this->reciprocalLatticeDirection(this->C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& pn : planeNormals)
                {
                    if(std::abs(pn->planeSpacing()-d110)<FLT_EPSILON)
                    {// a 111 plane
                        gsMap.emplace(pn.get(),gammaSurface110);
                    }
                }
                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("chi",gsMap));
                temp.emplace(sp->sID,sp);
            }
            else if(sp=="sigma")
            {
                Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors110(TextFileParser(material.materialFile).readMatrixCols<double>("sigmaWaveVectors110",2,true));
                Eigen::Matrix<double,Eigen::Dynamic,3> f110(TextFileParser(material.materialFile).readMatrixCols<double>("sigmaGammaSurfacePoints110",3,true));
                const int rotSymm110(1);
                std::vector<Eigen::Matrix<double,2,1>> mirSymm110;
                mirSymm110.push_back((Eigen::Matrix<double,2,1>()<<1.0,0.0).finished()); // symm with respect to local y-axis
                mirSymm110.push_back((Eigen::Matrix<double,2,1>()<<0.0,1.0).finished()); // symm with respect to local x-axis
                const Eigen::Matrix<double,2,2> A110((Eigen::Matrix<double,2,2>()<< 1.0,-0.5,
                                                                                 0.0,0.5*std::sqrt(3.0)).finished());
                std::shared_ptr<GammaSurface> gammaSurface110(new GammaSurface(A110,waveVectors110,f110,rotSymm110,mirSymm110));
                const double d110(this->reciprocalLatticeDirection(this->C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& pn : planeNormals)
                {
                    if(std::abs(pn->planeSpacing()-d110)<FLT_EPSILON)
                    {// a 111 plane
                        gsMap.emplace(pn.get(),gammaSurface110);
                    }
                }
                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("sigma",gsMap));
                temp.emplace(sp->sID,sp);
            }

            else
            {
                throw std::runtime_error("Unnown SecondPhase "+sp+" in BCC crystals.");
            }
        }
        return temp;
    }
        
} // namespace model
#endif

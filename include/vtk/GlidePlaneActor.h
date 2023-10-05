/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneActor_h_
#define model_GlidePlaneActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QComboBox>
#include <QGroupBox>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>

#include <IDreader.h>
//#include <PlanarPolygon.h>
#include <DDauxIO.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <Polycrystal.h>
#include <GlidePlaneFactory.h>
#include <GlidePlaneNoise.h>
#include <ConfigurationFields.h>
#include <TriangularMesh.h>
#include <CompareVectorsByComponent.h>

namespace model
{

    struct SingleGlidePlaneActor
    {
        typedef CompareVectorsByComponent<double,2,float> CompareType;
        std::map<Eigen::Matrix<double,2,1>,int,CompareType> uniquePointsIDs;

        std::vector<Eigen::Matrix<double,2,1>> points;
        std::vector<Eigen::Matrix<int,2,1>> segments;
        
        const GlidePlane<3>& glidePlane;
        
        SingleGlidePlaneActor(const GlidePlane<3>& glidePlane_in);
        void appendClosedPolygon(const std::vector<Eigen::Matrix<double,2,1>>& newPoints);
        
    };

    struct GlidePlaneActor : public QWidget
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        vtkRenderer* const renderer;
//        const Polycrystal<3>& poly;
        ConfigurationFields<3>& configFields;

        QGridLayout* mainLayout;
        QGroupBox* glidePlanesGroup;
//        QCheckBox* showGlidePlanes;
        
        QGroupBox* glidePlanesNoiseGroup;

//        QCheckBox* showGlidePlanesNoise;
        QComboBox* glidePlanesNoiseBox;
        QComboBox* slipSystemNoiseBox;
        
        QLineEdit* ssNoiseMin;
        QLineEdit* ssNoiseMax;
        QLineEdit* sfNoiseMin;
        QLineEdit* sfNoiseMax;
        
        QGroupBox* glidePlaneMeshGroup;



        vtkSmartPointer<vtkPolyData> glidePlanePolydata;
        vtkSmartPointer<vtkPolyDataMapper> glidePlaneMapper;
        vtkSmartPointer<vtkActor> glidePlaneActor;
        
        
        vtkSmartPointer<vtkPolyData> meshPolydata;
        vtkSmartPointer<vtkPolyDataMapper> meshMapper;
        vtkSmartPointer<vtkActor> meshActor;

        
        std::map<size_t,std::vector<vtkSmartPointer<vtkDataSetMapper>>> noiseMappers;
        std::map<size_t,std::vector<vtkSmartPointer<vtkActor>>> noiseActors;


//        GlidePlaneFactory<3> glidePlaneFactory;
//        GlidePlaneNoise planeNoise;
        
        Eigen::Array<double,2,2> noiseLimits;
//        Eigen::Array<double,2,1> sfNoiseLimits;

//        std::map<LatticePlaneKey<3>,SingleGlidePlaneActor> singleGlidePlaneMap;

        
        public:
                        
//        GlidePlaneActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,const Polycrystal<3>& poly,const DDtraitsIO& traitsIO);
////        void updateConfiguration(const DDauxIO<3>& auxIO);
//        void updateConfiguration(const DDconfigIO<3>& configIO);

        GlidePlaneActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,ConfigurationFields<3>& configFields_in);
//        void updateConfiguration(const DDauxIO<3>& auxIO);
        void updateConfiguration();
        
    };
    
} // namespace model
#endif

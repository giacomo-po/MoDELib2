/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneActor_cpp_
#define model_GlidePlaneActor_cpp_

#include <numbers>
#include <vtkLine.h>
#include <vtkPlanes.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkStructuredGridAppend.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>

#include <GlidePlaneActor.h>
#include <Polygon2D.h>

namespace model
{

SingleGlidePlaneActor::SingleGlidePlaneActor(const GlidePlane<3>& glidePlane_in) :
/* init */ glidePlane(glidePlane_in)
{
    
    
}

void SingleGlidePlaneActor::appendClosedPolygon(const std::vector<Eigen::Matrix<double,2,1>>& newPoints)
{
    size_t nextPointID(points.size());
    
    for(size_t k=0;k<newPoints.size();++k)
    {
        
        const auto iterK(uniquePointsIDs.find(newPoints[k]));
        int kID(nextPointID);
        if(iterK!=uniquePointsIDs.end())
        {// newPoints[k] exists
            kID=iterK->second;
        }
        else
        {// newPoints[k] does not exists
            points.emplace_back(newPoints[k]);
            uniquePointsIDs.emplace(newPoints[k],nextPointID);
            nextPointID++;
        }
        
        const size_t k1(k<newPoints.size()-1? k+1 : 0);
        const auto iterK1(uniquePointsIDs.find(newPoints[k1]));
        int k1ID(nextPointID);
        if(iterK1!=uniquePointsIDs.end())
        {// newPoints[k1] exists
            k1ID=iterK1->second;
        }
        else
        {// newPoints[k1] does not exists
            points.emplace_back(newPoints[k1]);
            uniquePointsIDs.emplace(newPoints[k1],nextPointID);
            nextPointID++;
        }
        
        
        segments.emplace_back((Eigen::Matrix<int,2,1>()<<kID,k1ID).finished());
    }
    
}

GlidePlaneActor::GlidePlaneActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,ConfigurationFields<3>& configFields_in) :
/* init */ renderWindow(renWin)
/* init */,renderer(ren)
/* init */,configFields(configFields_in)
/* init */,mainLayout(new QGridLayout(this))
/* init */,glidePlanesGroup(new QGroupBox(tr("&Planes")))
/* init */,glidePlanesNoiseGroup(new QGroupBox(tr("&Noise")))
/* init */,glidePlanesNoiseBox(new QComboBox(this))
/* init */,slipSystemNoiseBox(new QComboBox(this))
/* init */,ssNoiseMin(new QLineEdit("0.0"))
/* init */,ssNoiseMax(new QLineEdit("0.0"))
/* init */,sfNoiseMin(new QLineEdit("0.0"))
/* init */,sfNoiseMax(new QLineEdit("0.0"))
/* init */,glidePlaneMeshGroup(new QGroupBox(tr("&Stacking Faults")))
/* init */,glidePlanePolydata(vtkSmartPointer<vtkPolyData>::New())
/* init */,glidePlaneMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,glidePlaneActor(vtkSmartPointer<vtkActor>::New())
/* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
/* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,meshActor(vtkSmartPointer<vtkActor>::New())
/* init */,noiseLimits(Eigen::Array<double,2,2>::Zero())
{


    
    glidePlanesGroup->setCheckable(true);
    glidePlanesGroup->setChecked(false);
    glidePlaneActor->SetVisibility(false);
    
    glidePlanesNoiseGroup->setCheckable(true);
    glidePlanesNoiseGroup->setChecked(false);

    glidePlaneMeshGroup->setCheckable(true);
    glidePlaneMeshGroup->setChecked(false);
    meshActor->SetVisibility(false);

    
//    showGlidePlanesNoise->setChecked(false);
//    showGlidePlanesNoise->setText("show GlidePlanesNoise");
//    glidePlanesNoiseBox->setEnabled(false);
//    slipSystemNoiseBox->setEnabled(false);
    
    glidePlanesNoiseBox->addItem("RSS");
    glidePlanesNoiseBox->addItem("ISF energy");
    const auto& grain(configFields.ddBase.poly.grains.begin()->second);
    
    for(size_t k=0; k<grain.singleCrystal->slipSystems().size();++k)
    {
        slipSystemNoiseBox->addItem(QString::fromStdString("slipSystem "+ std::to_string(k)));
    }
    
    QGridLayout *noiseLayout = new QGridLayout();
    glidePlanesNoiseGroup->setLayout(noiseLayout);
    noiseLayout->addWidget(glidePlanesNoiseBox,1,1,1,1);
    noiseLayout->addWidget(slipSystemNoiseBox,2,1,1,1);
    noiseLayout->addWidget(ssNoiseMin,3,0,1,1);
    noiseLayout->addWidget(ssNoiseMax,3,1,1,1);
    noiseLayout->addWidget(sfNoiseMin,4,0,1,1);
    noiseLayout->addWidget(sfNoiseMax,4,1,1,1);

    
    mainLayout->addWidget(glidePlanesGroup,0,0,1,1);
    mainLayout->addWidget(glidePlanesNoiseGroup,1,0,1,1);
    mainLayout->addWidget(glidePlaneMeshGroup,2,0,1,1);

    
    
    
    
    this->setLayout(mainLayout);
    connect(glidePlanesGroup,SIGNAL(toggled(bool)), this, SLOT(modify()));
    connect(glidePlanesNoiseGroup,SIGNAL(toggled(bool)), this, SLOT(modify()));
    connect(glidePlaneMeshGroup,SIGNAL(toggled(bool)), this, SLOT(modify()));

    connect(glidePlanesNoiseBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
    connect(slipSystemNoiseBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
    connect(ssNoiseMin,SIGNAL(returnPressed()), this, SLOT(modify()));
    connect(ssNoiseMax,SIGNAL(returnPressed()), this, SLOT(modify()));
    connect(sfNoiseMin,SIGNAL(returnPressed()), this, SLOT(modify()));
    connect(sfNoiseMax,SIGNAL(returnPressed()), this, SLOT(modify()));
    
    
    
    glidePlaneMapper->SetInputData(glidePlanePolydata);
    glidePlaneActor->SetMapper(glidePlaneMapper);
    glidePlaneActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)
    renderer->AddActor(glidePlaneActor);

    meshPolydata->Allocate();
    meshMapper->SetInputData(meshPolydata);
    meshActor->SetMapper ( meshMapper );
    meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.
    renderer->AddActor(meshActor);

    
    //        renderer->AddActor(solidSolutionNoiseActorXZ);
    
    
}

//void GlidePlaneActor::updateConfiguration(const DDconfigIO<3>& configIO)
void GlidePlaneActor::updateConfiguration()
{
    
    std::cout<<"Updating GlidePlanes..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    
    
    for(const auto& actorPair : noiseActors)
    {
        for(const auto& actor : actorPair.second)
        {
            renderer->RemoveActor(actor);
        }
    }
    noiseMappers.clear();
    noiseActors.clear();
    
    vtkSmartPointer<vtkPoints> glidePlanePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> glidePlaneCells(vtkSmartPointer<vtkCellArray>::New());
    //        vtkNew<vtkStructuredGridAppend> gridsAppend;
    
    
    size_t pointIDs(0);
    //        double noiseSSMin(std::numeric_limits<double>::max());
    //        double noiseSSMax(std::numeric_limits<double>::min());
    //        double noiseSSMin(std::numeric_limits<double>::max());
    //        double noiseSSMax(std::numeric_limits<double>::min());
    
    noiseLimits<<std::numeric_limits<double>::max(),std::numeric_limits<double>::min(),
    /*         */std::numeric_limits<double>::max(),std::numeric_limits<double>::min();
    

    for(const auto& patch : configFields.loopPatches())
    {
        for(const auto& pair : patch.second.localPatches())
        {
            const auto& glidePlane(pair.first->glidePlane);
            const auto& grain(glidePlane->grain);
    //    for(const auto& loopIO : configFields.configIO.loops())
//    {
        // Plot GlidePlaneBoundary
        //const auto& grain(configFields.ddBase.poly.grain(loopIO.grainID));
        //GlidePlaneKey<3> planeKey(loopIO.P, grain.singleCrystal->reciprocalLatticeDirection(loopIO.N));
        //const auto glidePlane(configFields.ddBase.glidePlaneFactory.getFromKey(planeKey));
        if(glidePlane->slipSystems().size())
        {
            
            // Plot intersections with mesh
            for(const auto& mshInt : glidePlane->meshIntersections)
            {
                glidePlanePoints->InsertNextPoint(mshInt->P0.data());
                
                glidePlanePoints->InsertNextPoint(mshInt->P1.data());
                
                vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                line->GetPointIds()->SetId(0, pointIDs); // the second 0 is the index of the Origin in linesPolyData's points
                line->GetPointIds()->SetId(1, pointIDs+1);
                glidePlaneCells->InsertNextCell(line);
                pointIDs+=2;
            }
            
            // Plot intersectins with SphericalInclusions
            for(const auto& siIO : configFields.configIO.sphericalInclusions())
            {
                const Eigen::Matrix<double,3,1> planeC(glidePlane->snapToPlane(siIO.C));
                const Eigen::Matrix<double,3,1> dC(planeC-siIO.C);
                const double R2(std::pow(siIO.a,2));
                const double dC2(dC.squaredNorm());
                if(dC2<R2)
                {// Plane intersects inclusion
                    const double r(std::sqrt(R2-dC2));
                    const int NP(50);
                    const double dTheta(2.0*std::numbers::pi/NP);
                    for(int k=0;k<NP;++k)
                    {
                        const Eigen::Matrix<double,3,1> circPoint0(glidePlane->globalPosition((Eigen::Matrix<double,2,1>()<<r*std::cos(    k*dTheta),r*std::sin(    k*dTheta)).finished())-glidePlane->P+planeC);
                        const Eigen::Matrix<double,3,1> circPoint1(glidePlane->globalPosition((Eigen::Matrix<double,2,1>()<<r*std::cos((k+1)*dTheta),r*std::sin((k+1)*dTheta)).finished())-glidePlane->P+planeC);
                        
                        glidePlanePoints->InsertNextPoint(circPoint0.data());
                        glidePlanePoints->InsertNextPoint(circPoint1.data());
                        vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
                        line->GetPointIds()->SetId(0, pointIDs); // the second 0 is the index of the Origin in linesPolyData's points
                        line->GetPointIds()->SetId(1, pointIDs+1);
                        glidePlaneCells->InsertNextCell(line);
                        pointIDs+=2;
                    }
                }
            }
            
            
            auto clipNormals = vtkSmartPointer<vtkDoubleArray>::New();
            clipNormals->SetNumberOfComponents(3);
            vtkNew<vtkPoints> clipPts;
            for(const auto& meshInt : glidePlane->meshIntersections)
            {
                clipNormals->InsertNextTuple(meshInt->outNormal.data());
                clipPts->InsertNextPoint(meshInt->P0.data());
            }
            vtkNew<vtkPlanes> clipPlanes;
            clipPlanes->SetNormals(clipNormals);
            clipPlanes->SetPoints(clipPts);
            
            // Plot solidSolution-noise
//            for(size_t slipSystemID=0; slipSystemID<configFields.ddBase.poly.grain(loopIO.grainID).singleCrystal->slipSystems().size();++slipSystemID)
            for(size_t slipSystemID=0; slipSystemID<grain.singleCrystal->slipSystems().size();++slipSystemID)
            {
//                const auto& slipSystem(configFields.ddBase.poly.grain(loopIO.grainID).singleCrystal->slipSystems()[slipSystemID]);
                const auto& slipSystem(grain.singleCrystal->slipSystems()[slipSystemID]);
                if(slipSystem->planeNoise)
                {
                    if(slipSystem->planeNoise->solidSolution || slipSystem->planeNoise->stackingFault)
                    {
                        if(glidePlane->n.cross(slipSystem->n).squaredNorm()==0)
                        {// glide plane includes slip system
                            
                            std::set<double> xPos;
                            std::set<double> yPos;
                            for(const auto& meshInt : glidePlane->meshIntersections)
                            {
                                const auto localPos(slipSystem->globalToLocal(meshInt->P0-glidePlane->P));
                                xPos.insert(localPos(0));
                                yPos.insert(localPos(1));
                            }
                            
                            const Eigen::Array<double,2,1> lCorner(*xPos. begin(),*yPos. begin());
                            const Eigen::Array<double,2,1> hCorner(*xPos.rbegin(),*yPos.rbegin());
                            
                            const Eigen::Array<int,2,1> lIdx(slipSystem->planeNoise->posToIdx(lCorner).first);
                            const Eigen::Array<int,2,1> hIdx(slipSystem->planeNoise->posToIdx(hCorner).second);
                            
                            const int Nx(hIdx(0)-lIdx(0)+1);
                            const int Ny(hIdx(1)-lIdx(1)+1);
                            
                            vtkNew<vtkPoints> glidePlaneNoisePoints;
                            glidePlaneNoisePoints->SetNumberOfPoints(Nx*Ny);
                            
                            vtkNew<vtkDoubleArray> glidePlaneSSNoise;
                            glidePlaneSSNoise->SetNumberOfValues(Nx*Ny);
                            
                            vtkNew<vtkDoubleArray> glidePlaneSFNoise;
                            glidePlaneSFNoise->SetNumberOfValues(Nx*Ny);
                            
                            for(int i=0;i<Nx;++i)
                            {
                                const int gridi(i+lIdx(0));
                                
                                for(int j=0;j<Ny;++j)
                                {
                                    const int gridj(j+lIdx(1));
                                    //                                        const int storageIdx = Ny*i + j;
                                    const int storageIdx = Nx*j + i;
                                    
                                    const Eigen::Array<int,2,1> idx(gridi,gridj);
                                    const Eigen::Array<double,2,1> localPos(slipSystem->planeNoise->idxToPos(idx));
                                    
                                    const Eigen::Array<double,3,1> globalPos(slipSystem->localToGlobal(localPos.matrix())+glidePlane->P);
                                    glidePlaneNoisePoints->SetPoint(storageIdx,globalPos.data());
                                    
                                    const auto noiseVal(slipSystem->gridVal(idx));
                                    glidePlaneSSNoise->SetValue(storageIdx,std::get<1>(noiseVal));
                                    glidePlaneSFNoise->SetValue(storageIdx,std::get<2>(noiseVal));
                                    noiseLimits<<std::min(noiseLimits(0,0),std::get<1>(noiseVal)),std::max(noiseLimits(0,1),std::get<1>(noiseVal)),
                                    /*         */std::min(noiseLimits(1,0),std::get<2>(noiseVal)),std::max(noiseLimits(1,1),std::get<2>(noiseVal));
                                    
                                }
                            }
                            
                            
                            ssNoiseMin->setText(QString::number(noiseLimits(0,0)));
                            ssNoiseMax->setText(QString::number(noiseLimits(0,1)));
                            sfNoiseMin->setText(QString::number(noiseLimits(1,0)));
                            sfNoiseMax->setText(QString::number(noiseLimits(1,1)));
                            
                            
                            vtkNew<vtkStructuredGrid> glidePlaneNoiseGrid;
                            
                            glidePlaneNoiseGrid->SetDimensions(Nx,Ny,1);
                            glidePlaneNoiseGrid->SetPoints(glidePlaneNoisePoints);
                            glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneSSNoise);
                            glidePlaneNoiseGrid->GetPointData()->AddArray(glidePlaneSFNoise);
                            
                            
                            vtkNew<vtkTableBasedClipDataSet> clipper;
                            clipper->SetInputData(glidePlaneNoiseGrid);
                            clipper->SetClipFunction(clipPlanes);
                            clipper->InsideOutOn();
                            
                            noiseMappers[slipSystemID].push_back(vtkSmartPointer<vtkDataSetMapper>::New());
                            noiseMappers[slipSystemID].back()->SetInputConnection(clipper->GetOutputPort());
                            noiseMappers[slipSystemID].back()->SetScalarModeToUsePointFieldData();
                            noiseMappers[slipSystemID].back()->SelectColorArray(glidePlanesNoiseBox->currentIndex());
                            noiseMappers[slipSystemID].back()->SetScalarRange(noiseLimits(glidePlanesNoiseBox->currentIndex(),0), noiseLimits(glidePlanesNoiseBox->currentIndex(),1));
                            noiseMappers[slipSystemID].back()->ScalarVisibilityOn();
                            noiseActors[slipSystemID].push_back(vtkSmartPointer<vtkActor>::New());
                            noiseActors[slipSystemID].back()->SetMapper(noiseMappers[slipSystemID].back());
                            noiseActors[slipSystemID].back()->SetVisibility(glidePlanesNoiseGroup->isChecked() && slipSystemNoiseBox->currentIndex()==int(slipSystemID));
                            renderer->AddActor( noiseActors[slipSystemID].back());
                        }
                    }
                }
            }
            
            
        }
    }
}
    
    // Update StackingFaults
//    singleGlidePlaneMap.clear();
    if(glidePlaneMeshGroup->isChecked())
    {
        std::map<LatticePlaneKey<3>,SingleGlidePlaneActor> singleGlidePlaneMap;
        for(const auto& patch : configFields.loopPatches())
        {
            for(const auto& pair : patch.second.localPatches())
            {
                const auto& glidePlane(pair.first->glidePlane);
                const auto& key(glidePlane->key);
                const auto glidePlaneIter(singleGlidePlaneMap.find(key));
                if(glidePlaneIter!=singleGlidePlaneMap.end())
                {// glide plane not found
                    glidePlaneIter->second.appendClosedPolygon(pair.second);
                }
                else
                {// glide plane not found
                    const auto success(singleGlidePlaneMap.emplace(std::piecewise_construct,
                                                                   std::forward_as_tuple(key),
                                                                   std::forward_as_tuple(*glidePlane)));
                    if(success.second)
                    {
                        success.first->second.appendClosedPolygon(pair.second);
                    }
                    else
                    {
                        throw std::runtime_error("Cannot insert glide plane in singleGlidePlaneMap");
                    }
                }
            }
        }
        
        vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
        vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
        vtkSmartPointer<vtkUnsignedCharArray> meshColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
        meshColors->SetNumberOfComponents(3);

    //    std::cout<<"singleGlidePlaneMap.size() ="<<singleGlidePlaneMap.size()<<std::endl;
        size_t numVertices(0);
        for(auto& pair : singleGlidePlaneMap)
        {
                TriangularMesh triMesh;
                triMesh.reMesh(pair.second.points,pair.second.segments,100000.0);
                
                for(const auto& point2d : triMesh.vertices())
                {
                    const auto point3d(pair.second.glidePlane.globalPosition(point2d));
                    meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
                }
                
                for(const auto& tri : triMesh.triangles())
                {
                    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                    triangle->GetPointIds()->SetId (0,tri(0)+numVertices);
                    triangle->GetPointIds()->SetId (1,tri(1)+numVertices);
                    triangle->GetPointIds()->SetId (2,tri(2)+numVertices);
                    meshTriangles->InsertNextCell ( triangle );
                    const auto triColor(Eigen::Matrix<int,1,3>::Random()*255);
                    meshColors->InsertNextTuple3(triColor(0),triColor(1),triColor(2)); // use this to assig color to each vertex
                }
                numVertices+=triMesh.vertices().size();
        }
    //    std::cout<<"Done meshing"<<std::endl;
        
        meshPolydata->SetPoints ( meshPts );
        meshPolydata->SetPolys ( meshTriangles );
        meshPolydata->GetCellData()->SetScalars(meshColors);
        meshPolydata->Modified();
        meshMapper->SetScalarModeToUseCellData();
    //    renWin->Render();
    }
    
    glidePlanePolydata->SetPoints(glidePlanePoints);
    glidePlanePolydata->SetLines(glidePlaneCells);
    glidePlanePolydata->Modified();
    
    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void GlidePlaneActor::modify()
{
    glidePlaneActor->SetVisibility(glidePlanesGroup->isChecked());
    meshActor->SetVisibility(glidePlaneMeshGroup->isChecked());

    
//    glidePlanesNoiseBox->setEnabled(glidePlanesNoiseGroup->isChecked());
//    slipSystemNoiseBox->setEnabled(glidePlanesNoiseGroup->isChecked());
    
    // Select solidSolution or stacking-fault noise
    for(const auto& mapperPair : noiseMappers)
    {
        for(const auto& mapper : mapperPair.second)
        {
            mapper->SelectColorArray(glidePlanesNoiseBox->currentIndex());
            
            switch (glidePlanesNoiseBox->currentIndex())
            {
                case 0:
                    mapper->SetScalarRange(ssNoiseMin->text().toDouble(), ssNoiseMax->text().toDouble());
                    break;
                    
                case 1:
                    mapper->SetScalarRange(sfNoiseMin->text().toDouble(), sfNoiseMax->text().toDouble());
                    break;
                    
                default:
                    break;
            }
        }
    }
    
    const size_t slipSystemID(slipSystemNoiseBox->currentIndex());
    for(const auto& actorPair : noiseActors)
    {
        for(const auto& actor : actorPair.second)
        {
            actor->SetVisibility(actorPair.first==slipSystemID && glidePlanesNoiseGroup->isChecked());
        }
    }
    renderWindow->Render();
}


} // namespace model
#endif

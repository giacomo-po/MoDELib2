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

#ifndef model_PyGlideSolver_cpp_
#define model_PyGlideSolver_cpp_

#include <filesystem>
#include <PyGlideSolver.h>
#include <TerminalColors.h>


namespace model
{
    

        
    #ifdef _MODEL_PYBIND11_ // COMPLIED WITH PYBIND11

    template <typename DislocationNetworkType>
    PyGlideSolver<DislocationNetworkType>::PyGlideSolver(const DislocationNetworkType& DN_in) :
    /* init */ DislocationGlideSolverBase<DislocationNetworkType>(DN_in)
    /* init */,pyModuleName(this->DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/"+TextFileParser(this->DN.ddBase.simulationParameters.traitsIO.ddFile).readString("pyModuleName"))
    {
        std::filesystem::path pyModulePath(pyModuleName);
        std::string pyModuleDir(pyModulePath.parent_path());
        std::string pyModuleFile(pyModulePath.filename());
        std::cout<<greenBoldColor<<"Creating PyGlideSolver: dir= "<<pyModuleDir<<", file="<<pyModuleFile<<defaultColor<<std::endl;
//        pybind11::scoped_interpreter guard{};
        pybind11::module sys = pybind11::module::import("sys");
        pybind11::list path = sys.attr("path");
        path.append(pyModuleDir.c_str());

        pyModule=pybind11::module::import(pyModuleFile.c_str());
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd PyGlideSolver<DislocationNetworkType>::getNodeVelocities() const
    {
    const auto t0= std::chrono::system_clock::now();
    std::cout<<"PyGlideSolver solving "<<std::flush;
    
    auto configIO(this->DN.io().configIO()); //evl files
    const auto auxIO(this->DN.io().auxIO()); //DDaux files
    const auto Temperature(this->DN.ddBase.poly.T);

    // PULL FROM MEMORY
    
    //MAPS NEEDED
    std::map<size_t, size_t> nodeIDmap; // key=nodeID, value=row index in file
    size_t k = 0;
    for (const auto &node : configIO.nodes())
    {
        if (node.meshLocation!=-2)
        {
            nodeIDmap.emplace(node.sID, k);
            k++;
        }
    }
    std::map<std::pair<size_t, size_t>, std::map<size_t, DislocationQuadraturePointIO<3>>> qPointsMap;
    for (const auto &qp : auxIO.quadraturePoints())
    {
        const std::pair<size_t, size_t> outerKey(std::make_pair(qp.sourceID, qp.sinkID));
        const size_t innerKey(qp.qID);
        qPointsMap[outerKey].emplace(innerKey, qp);
    }
    
    // QP SIZE
    int qpSize(0);
    for (const auto &outerPair : qPointsMap)
    { // write all quadrature points positions
        qpSize+=outerPair.second.size();
//        for (const auto &innerPair : outerPair.second)
//        {
//            qpSize+=1;
//        }
    }
    
    
    // POSITIONS
    Eigen::MatrixXd positionNodes(Eigen::MatrixXd::Zero(this->DN.networkNodes().size()+qpSize,this->DN.NdofXnode));
    int nodeCount(0);
    for (const auto &node : configIO.nodes())
    {
        if (node.meshLocation!=-2)
        {
            positionNodes.row(nodeCount) << node.P.transpose()*this->DN.ddBase.poly.b_SI*1.0e10;
            nodeCount+=1;
        }
    }
    for (const auto &outerPair : qPointsMap)
    { // write all quadrature points positions
        for (const auto &innerPair : outerPair.second)
        {
            const auto &qp(innerPair.second);
            positionNodes.row(nodeCount) << qp.r.transpose()*this->DN.ddBase.poly.b_SI*1.0e10;
            nodeCount+=1;
        }
    }
//    std::cout << std::endl << "Size of node positions: " << positionNodes.size() << std::endl;
    //std::cout << positionNodes << std::endl;

    
    // STRESSES
    Eigen::MatrixXd stressNodes(Eigen::MatrixXd::Zero(this->DN.networkNodes().size()+qpSize,this->DN.NdofXnode*this->DN.NdofXnode));
    int stressCount(0);
    for (const auto &node : nodeIDmap)
    {
        stressNodes.row(stressCount) << Eigen::Matrix<double,3,3>::Zero().reshaped<Eigen::RowMajor>().transpose();
        stressCount+=1;
    }
    for (const auto &outerPair : qPointsMap)
    {
        for (const auto &innerPair : outerPair.second)
        {
            const auto &qp(innerPair.second);
            Eigen::Matrix<double, 3, 3> tempStress(Eigen::Matrix<double,3,3>::Zero());
            tempStress <<  qp.stress*this->DN.ddBase.poly.mu_SI*1e-9;
            stressNodes.row(stressCount) << tempStress.reshaped<Eigen::RowMajor>().transpose();
            stressCount+=1;
        }
    }
//    std::cout << "Size of node stresses: " << stressNodes.size() <<std::endl;
    //std::cout<< stressNodes << std::endl;

    
    // BURGERS VECTOR
    configIO.finalize();
    const auto segmentMap(configIO.segments());
    Eigen::MatrixXd burgersNodes(Eigen::MatrixXd::Zero(this->DN.networkNodes().size(),this->DN.NdofXnode));
    int burgersCount(0);
    for (const auto &outerPair : qPointsMap)
    {
        burgersNodes.row(burgersCount) << segmentMap.at(outerPair.first).b.transpose()*this->DN.ddBase.poly.b_SI*1.0e10;
        //burgersNodes.push_back(segmentMap.at(outerPair.first).b.transpose()*this->DN.poly.b_SI*1.0e10);
        burgersCount+=1;
    }
//    std::cout << "Size of node burgers vectors: " << burgersNodes.size() <<std::endl;
    //std::cout << burgersNodes << std::endl;

    
    // CONNECTIONS
    std::list<std::list<int>> connections;
    std::list<int> single_list;
//    std::cout << "qPointsmap size:" << qPointsMap.size() << std::endl;
    size_t qpCounter(0);
    for (const auto &pair : qPointsMap)
    {
        const auto &sourceNodeID(pair.first.first);
        const auto &sinkNodeID(pair.first.second);
        const auto sourceNodeRow(nodeIDmap.at(sourceNodeID));
        const auto sinkNodeRow(nodeIDmap.at(sinkNodeID));
        //std::cout << "Source Node Row" << sourceNodeRow <<std::endl;
        single_list.push_back(sourceNodeRow);
    
        for (size_t k = 0; k < pair.second.size(); ++k)
        {
            single_list.push_back(qpCounter + k + nodeIDmap.size());
        }
        //std::cout << "Sink Node Row" << sinkNodeRow <<std::endl;
        single_list.push_back(sinkNodeRow);
        qpCounter += pair.second.size();
        connections.push_back(single_list); //Append single list into the connections
        single_list.erase(single_list.begin(),single_list.end());
    }
    //printNestedList(connections);

    
    // CALL TO ML MODEL WITH PYBIND
    Eigen::MatrixXd nodePythonVelocity(pyModule.attr("GlideSolver")().attr("getNodeVelocities")(positionNodes,stressNodes,burgersNodes,connections,Temperature).template cast<Eigen::MatrixXd>());

    std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
   
    return nodePythonVelocity.reshaped<Eigen::RowMajor>();
}

    #else // COMPLIED WITHOUT PYBIND11

    template <typename DislocationNetworkType>
    PyGlideSolver<DislocationNetworkType>::PyGlideSolver(const DislocationNetworkType& DN_in) :
    /* init */ DislocationGlideSolverBase<DislocationNetworkType>(DN_in)
    {
        std::cout<<greenBoldColor<<"Creating PyGlideSolver"<<defaultColor<<std::endl;
        throw std::runtime_error("PyGlideSolver must be compiled with pybind11");
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd PyGlideSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        return Eigen::VectorXd();
    }
    #endif


    template class PyGlideSolver<DislocationNetwork<3,0>>;
    
}
#endif

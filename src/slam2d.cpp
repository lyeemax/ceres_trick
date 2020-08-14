//
// Created by unicorn on 2020/8/14.
//

#include "readg2o.h"
#include <ceres/ceres.h>
#include <posegraph2d_factor.h>
#include "posegraph2d_factor.h"

using namespace std;

Vertex_SE2T Vertexs;
EDGE_SE2T Edges;
void assignToArr(Vertex_SE2T &Vertexs,double PQ[][3],int id){
    assert(Vertexs.find(id)!=Vertexs.end());
    PQ[id][0]=Vertexs[id].p.x();
    PQ[id][1]=Vertexs[id].p.y();
    PQ[id][2]=Vertexs[id].angle;
}
int main() {
    const string file = "/home/unicorn/CLionProjects/sfm/dataset/city10000.g2o";
    auto res = readg2oSE2(file, Vertexs, Edges);
    if (res) {
        cout << Vertexs.size() << endl;
        cout << Edges.size() << endl;
    }
    int N=Vertexs.size();
    double PQ[N][3];
    ceres::Problem problem;
    ceres::LossFunction *cauchy=new ceres::CauchyLoss(1.0);
    AngleLocalParamter *lp=new AngleLocalParamter();
    //ceres::LocalParameterization *lp=AngleLocalParameterization::Create();

    auto iter=Edges.begin();
    for(;iter!=Edges.end();iter++){
        auto id1=iter->second.ID1;
        assert(id1==iter->first);
        auto id2=iter->second.ID2;
        assignToArr(Vertexs,PQ,id1);
        assignToArr(Vertexs,PQ,id2);
        //auto sqrtInfo=iter->second.InfoMatrix.llt().matrixL();
        ceres::CostFunction *poseGraph=new PoseGraph2D(iter->second.t12.x(),iter->second.t12.y(),iter->second.R12);
        ceres::CostFunction *poseGraphAuto=PoseGraph2DAutoDiff::Create(iter->second.t12.x(),iter->second.t12.y(),iter->second.R12);
        double *param[2]={PQ[id1],PQ[id2]};
        Eigen::Vector3d res1,res2;
        res1.setZero();
        res2.setZero();
        Eigen::Matrix<double,9,1> autojac1,autojac2;
        Eigen::Matrix<double,9,1> jac1,jac2;
        double *Jac1[2]={autojac1.data(),autojac2.data()};
        double *Jac2[2]={jac1.data(),jac2.data()};
        poseGraphAuto->Evaluate(param,res1.data(),Jac1);
        //poseGraph->Evaluate(param,res2.data(),Jac2);


       // cout<<"auto jac is "<<endl<<Eigen::Map<Eigen::Matrix<double,9,1>,Eigen::RowMajor>(Jac1[0]).transpose()<<endl<<Eigen::Map<Eigen::Matrix<double,9,1>,Eigen::RowMajor>(Jac1[1]).transpose()<<endl;

//        cout<<"--------"<<endl;
//        cout<<"auto res"<<endl<<res1<<endl;
//        cout<<"my res "<<endl<<res2<<endl;
        problem.AddResidualBlock(poseGraph,cauchy,PQ[id1],PQ[id2]);
        problem.SetParameterization(PQ[id1],lp);
        problem.SetParameterization(PQ[id2],lp);
    }
    problem.SetParameterBlockConstant(PQ[0]);

    savePose2D("before.txt",PQ,N);
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    //options.num_threads = 2;
    //options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 200;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout<<summary.FullReport()<<endl;

    savePose2D("after.txt",PQ,N);
}
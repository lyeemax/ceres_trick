//
// Created by unicorn on 2020/8/12.
//
#include <iostream>
#include "readg2o.h"
#include <string>
#include <ceres/ceres.h>
#include "PoseLocalParameter.h"
#include "posegraph_factor.h"
#include <vector>
using namespace std;

void assignToArr(Vertex_SE3T &Vertex3s,double PQ[][7],int id){
    assert(Vertex3s.find(id)!=Vertex3s.end());
        PQ[id][0]=Vertex3s[id].p.x();
        PQ[id][1]=Vertex3s[id].p.y();
        PQ[id][2]=Vertex3s[id].p.z();
        PQ[id][3]=Vertex3s[id].q.x();
        PQ[id][4]=Vertex3s[id].q.y();
        PQ[id][5]=Vertex3s[id].q.z();
        PQ[id][6]=Vertex3s[id].q.w();

}
int main(){
    Vertex_SE3T Vertex3s;
    EDGE_SE3T Edge3s;
    const string se3file="/home/unicorn/CLionProjects/sfm/dataset/sphere2500.g2o";
    auto res3=readg2oSE3(se3file,Vertex3s,Edge3s);
    if(!res3){
        cout<<"read g2o file failed"<<endl;
    }
    int N=Vertex3s.size();
    double PQ[N][7];
    ceres::Problem problem;
    ceres::LossFunction *cauchy=new ceres::CauchyLoss(1.0);
    PoseLocalParameterizationMath *lp=new PoseLocalParameterizationMath();

    auto iter=Edge3s.begin();
    for(;iter!=Edge3s.end();iter++){
        auto id1=iter->second.ID1;
        assert(id1==iter->first);
        auto id2=iter->second.ID2;
        assignToArr(Vertex3s,PQ,id1);
        assignToArr(Vertex3s,PQ,id2);
        auto sqrtInfo=iter->second.InfoMatrix.llt().matrixL();
        ceres::CostFunction *poseGraph=PoseGraphAutoDiff::create(iter->second.t12,iter->second.R12,sqrtInfo);
        problem.AddResidualBlock(poseGraph,cauchy,PQ[id1],PQ[id2]);
        problem.SetParameterization(PQ[id1],lp);
        problem.SetParameterization(PQ[id2],lp);
    }
    problem.SetParameterBlockConstant(PQ[0]);

    savePose("before.txt",PQ,N);
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    //options.num_threads = 2;
    //options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 200;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout<<summary.FullReport()<<endl;

    savePose("after.txt",PQ,N);
}

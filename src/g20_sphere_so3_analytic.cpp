//
// Created by unicorn on 2020/8/16.
//

#ifndef CERES_DEBUG_TRICK_G20_SPHERE_SO3_ANALYTIC_SRC
#define CERES_DEBUG_TRICK_G20_SPHERE_SO3_ANALYTIC_SRC
#include <iostream>
#include "readg2o.h"
#include <string>
#include <ceres/ceres.h>
#include "so3_utils.h"
#include "posegraph_factor.h"
#include <vector>
using namespace std;

void assignToArr(Vertex_SE3T &Vertex3s,double PQ[][6],int id){
    assert(Vertex3s.find(id)!=Vertex3s.end());
    PQ[id][0]=Vertex3s[id].p.x();
    PQ[id][1]=Vertex3s[id].p.y();
    PQ[id][2]=Vertex3s[id].p.z();
    auto rvec=toRotationVec(Vertex3s[id].q);
    PQ[id][3]=rvec.x();
    PQ[id][4]=rvec.y();
    PQ[id][5]=rvec.z();
}
int main(){
    Vertex_SE3T Vertex3s;
    EDGE_SE3T Edge3s;
    const string se3file="../dataset/sphere2500.g2o";
    auto res3=readg2oSE3(se3file,Vertex3s,Edge3s);
    if(!res3){
        cout<<"read g2o file failed"<<endl;
    }
    int N=Vertex3s.size();
    double PQ[N][6];
    ceres::Problem problem;
    ceres::LossFunction *cauchy=new ceres::CauchyLoss(1.0);

    auto iter=Edge3s.begin();
    for(;iter!=Edge3s.end();iter++){
        auto id1=iter->second.ID1;
        assert(id1==iter->first);
        auto id2=iter->second.ID2;
        assignToArr(Vertex3s,PQ,id1);
        assignToArr(Vertex3s,PQ,id2);
        auto sqrtInfo=iter->second.InfoMatrix.llt().matrixL();
        PoseGraphSO3 *poseGraph=new PoseGraphSO3(iter->second.t12,iter->second.R12,sqrtInfo);
        problem.AddResidualBlock(poseGraph,cauchy,PQ[id1],PQ[id2]);
    }
    problem.SetParameterBlockConstant(PQ[0]);

    savePoseSO3("before.txt",PQ,N);
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    //options.num_threads = 2;
    //options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 200;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout<<summary.FullReport()<<endl;

    savePoseSO3("after.txt",PQ,N);


}

#endif //CERES_DEBUG_TRICK_G20_SPHERE_SO3_ANALYTIC_SRC

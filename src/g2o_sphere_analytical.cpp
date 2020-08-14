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

Vertex_SE3T Vertex3s;
EDGE_SE3T Edge3s;
double PQ[5000][7];
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
void test(){
    auto iter=Edge3s.begin();
    for(;iter!=Edge3s.end();iter++){
        auto id1=iter->second.ID1;
        auto id2=iter->second.ID2;
        assignToArr(Vertex3s,PQ,id1);
        assignToArr(Vertex3s,PQ,id2);
        auto sqrtInfo=iter->second.InfoMatrix.llt().matrixL();
        PoseGraph *poseGraph=new PoseGraph(iter->second.t12,iter->second.R12,sqrtInfo);
        double **param= new double*[2];
        param[0]=PQ[id1];
        param[1]=PQ[id2];
        poseGraph->check(param);
    }
}

int main(){
    const string se3file="../dataset/sphere2500.g2o";
    auto res3=readg2oSE3(se3file,Vertex3s,Edge3s);
    if(!res3){
        cout<<"read g2o file failed"<<endl;
    }
    int N=Vertex3s.size();
    //test();

    ceres::Problem problem;
    ceres::LossFunction *cauchy=new ceres::CauchyLoss(1.0);
    PoseLocalParameterization *lp=new PoseLocalParameterization();
    PoseLocalParameterizationMath *lp1=new PoseLocalParameterizationMath();

    auto iter=Edge3s.begin();
    for(;iter!=Edge3s.end();iter++){
        auto id1=iter->second.ID1;
        auto id2=iter->second.ID2;
        assignToArr(Vertex3s,PQ,id1);
        assignToArr(Vertex3s,PQ,id2);
        auto sqrtInfo=iter->second.InfoMatrix.llt().matrixL();
        //ceres::CostFunction *poseGrapgAuto=PoseGraphAutoDiff::create(iter->second.t12,iter->second.R12,sqrtInfo);

        PoseGraph *poseGraph=new PoseGraph(iter->second.t12,iter->second.R12,sqrtInfo);

//        Eigen::Matrix<double,7,6> jacqt1,jacqt2;
//        double *Jacqt[2]={jacqt1.data(),jacqt2.data()};
//        lp1->ComputeJacobian(PQ[id1],Jacqt[0]);
//        lp1->ComputeJacobian(PQ[id2],Jacqt[1]);
//
//        Eigen::Matrix<double,7,6> jac_analytic1,jac_analytic2;
//        double *JacAn[2]={jac_analytic1.data(),jac_analytic2.data()};
//        lp->ComputeJacobian(PQ[id1],JacAn[0]);
//        lp->ComputeJacobian(PQ[id2],JacAn[1]);
//
//        double *param[2]={PQ[id1],PQ[id2]};
//        Eigen::Matrix<double,6,1> res1,res2;
//        res1.setZero();
//        res2.setZero();
//        Eigen::Matrix<double,6,7> autojac1,autojac2;
//        Eigen::Matrix<double,6,7> jac1,jac2;
//        double *Jac1[2]={autojac1.data(),autojac2.data()};
//        double *Jac2[2]={jac1.data(),jac2.data()};
//        poseGrapgAuto->Evaluate(param,res1.data(),Jac1);
//        poseGraph->Evaluate(param,res2.data(),Jac2);
//
//        Eigen::Map<Eigen::Matrix<double,6,7,Eigen::RowMajor>> AJ1(Jac1[0]);
//        auto trueAJ1=AJ1*Eigen::Map<Eigen::Matrix<double,7,6,Eigen::RowMajor>>(Jacqt[0]);
//
//        Eigen::Map<Eigen::Matrix<double,6,7,Eigen::RowMajor>> AJ2(Jac1[1]);
//        auto trueAJ2=AJ2*Eigen::Map<Eigen::Matrix<double,7,6,Eigen::RowMajor>>(Jacqt[1]);
//        cout<<"auto jac is "<<endl<< trueAJ1<<endl
//        << trueAJ2<<endl;
//
//        cout<<"my jac is "<<endl<<Eigen::Map<Eigen::Matrix<double,6,7,Eigen::RowMajor>>(Jac2[0])*
//                Eigen::Map<Eigen::Matrix<double,7,6,Eigen::RowMajor>>(JacAn[0])<<endl
//            <<Eigen::Map<Eigen::Matrix<double,6,7,Eigen::RowMajor>>(Jac2[1])*
//                    Eigen::Map<Eigen::Matrix<double,7,6,Eigen::RowMajor>>(JacAn[1])<<endl;
//
//        cout<<"--------"<<endl;
//        cout<<"auto res"<<endl<<res1.transpose()<<endl;
//        cout<<"my res "<<endl<<res2.transpose()<<endl;

        problem.AddResidualBlock(poseGraph,cauchy,PQ[id1],PQ[id2]);
        problem.SetParameterization(PQ[id1],lp);
        problem.SetParameterization(PQ[id2],lp);
    }
    problem.SetParameterBlockConstant(PQ[0]);
    savePose("before.txt",PQ,N);
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 200;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout<<summary.FullReport()<<endl;


    savePose("after.txt",PQ,N);

}

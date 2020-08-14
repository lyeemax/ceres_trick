//
// Created by unicorn on 2020/8/12.
//
#include <iostream>
#include <curve_factor.h>
using namespace std;

int main() {
    Eigen::Matrix3d para;
    para << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    double px[]={15,22,32};
    double py[]={41,52,63};
    double pz[]={72,81,99};
    std::vector<Eigen::Vector3d> dataNosise;
    std::vector<Eigen::Vector3d> data;
    std::vector<Eigen::Vector3d> variables;
    gernerateData(100, para, dataNosise, data, variables);

    ceres::Problem problem;
    ceres::LossFunction * cauchy=new ceres::CauchyLoss(1.0);
    for (int i = 0; i <100 ; ++i) {
        CurveEdge *edge=new CurveEdge(variables[i],dataNosise[i]);
        problem.AddResidualBlock(edge,cauchy,px,py,pz);
    }
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 1000;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout<<summary.FullReport()<<endl;
    cout<<"px"<<px[0]<<" "<<px[1]<<" "<<px[2]<<endl;
    cout<<"py"<<py[0]<<" "<<py[1]<<" "<<py[2]<<endl;
    cout<<"pz"<<pz[0]<<" "<<pz[1]<<" "<<pz[2]<<endl;

}

//
// Created by unicorn on 2020/8/12.
//

#ifndef SFM_CURVE_FACTOR_H
#define SFM_CURVE_FACTOR_H

#include <ceres/ceres.h>
#include <random>
#include <Eigen/Core>
#include <vector>
#include <math.h>
//generate data with
//x         a3
//y =param* b2
//z         c1
void gernerateData(int n,Eigen::Matrix3d &para,
        std::vector<Eigen::Vector3d> &dataNosise,
        std::vector<Eigen::Vector3d> &data,
        std::vector<Eigen::Vector3d> &variables)
        {
            std::random_device rd;
            std::uniform_real_distribution<> dis(-20,20);
            std::normal_distribution<double> nd(0,3);
            for (int i = 0; i <n ; ++i) {
                Eigen::Vector3d abc(dis(rd),dis(rd),dis(rd));
                variables.push_back(abc);
                auto res=para*abc;
                data.push_back(res);
                auto nres=res+Eigen::Vector3d(nd(rd),nd(rd),nd(rd));
                dataNosise.push_back(nres);
            }
        }

class CurveEdge:public ceres::SizedCostFunction<3,3,3,3>{
public:
    CurveEdge(Eigen::Vector3d var,Eigen::Vector3d result){
        var_=var;
        result_=result;
    }
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const{
        Eigen::Map<Eigen::Vector3d> paramX(const_cast<double *>(parameters[0]));
        Eigen::Map<Eigen::Vector3d> paramY(const_cast<double *>(parameters[1]));
        Eigen::Map<Eigen::Vector3d> paramZ(const_cast<double *>(parameters[2]));
        Eigen::Matrix3d param;
        param<<paramX.transpose(),paramY.transpose(),paramZ.transpose();
        Eigen::Map<Eigen::Vector3d> res(residuals);
        res=result_-param*var_;

        if(jacobians){
            if(jacobians[0]){
                Eigen::Map<Eigen::Matrix<double,1,9>,Eigen::RowMajor> jac0(jacobians[0]);
                jac0.setZero();
                jac0.block<1,3>(0,0)=-var_.transpose();
            }
            if(jacobians[1]){
                Eigen::Map<Eigen::Matrix<double,1,9>,Eigen::RowMajor> jac1(jacobians[1]);
                jac1.setZero();
                jac1.block<1,3>(0,3)=-var_.transpose();
            }
            if(jacobians[2]){
                Eigen::Map<Eigen::Matrix<double,1,9>,Eigen::RowMajor> jac2(jacobians[2]);
                jac2.setZero();
                jac2.block<1,3>(0,6)=-var_.transpose();
            }
        }
        return true;
    }

    Eigen::Vector3d var_,result_;

};
#endif //SFM_CURVE_FACTOR_H

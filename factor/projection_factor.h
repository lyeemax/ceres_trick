//
// Created by unicorn on 2020/8/5.
//

#ifndef SFM_PROJECTION_H
#define SFM_PROJECTION_H

#include <ceres/ceres.h>
#include <ceres/problem.h>
#include "quaternion_utils.h"
using namespace ceres;
namespace ceres{
    class ProjectionEdge:public SizedCostFunction<3,7>{
    public:
        ProjectionEdge(const Eigen::Vector3d& pts_i,const Eigen::Vector3d& pts_j) {
            pts_i_=pts_i;
            pts_j_=pts_j;
        }
        virtual bool Evaluate(double const *const *parameters,double *residuals,double **jacobians)const{
            Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
            Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
            Eigen::Map<Eigen::Vector3d> residual(residuals);
            residual=pts_i_-Qi.toRotationMatrix()*pts_j_-Pi;

            if(jacobians){
                if(jacobians[0]){
                    Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobianRT(jacobians[0]);
                    Eigen::Matrix3d Jr=Qi.toRotationMatrix()*skew(pts_j_);
                    Eigen::Matrix<double,3,3> Jt=-Eigen::Matrix3d::Identity();
                    jacobianRT.leftCols(3)=Jr;
                    jacobianRT.block<3,3>(0,3)=Jt;
                    jacobianRT.rightCols(1).setZero();
                }
            }
            return true;
        }
        Eigen::Vector3d pts_i_,pts_j_;
    };
}


#endif //SFM_PROJECTION_H

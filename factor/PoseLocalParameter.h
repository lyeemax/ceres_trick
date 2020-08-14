//
// Created by unicorn on 2020/8/10.
//

#ifndef SFM_POSELOCALPARAMETER_H
#define SFM_POSELOCALPARAMETER_H
#include <eigen3/Eigen/Dense>
#include <ceres/ceres.h>
#include "quaternion_utils.h"

class PoseLocalParameterization : public ceres::LocalParameterization
{
public:
    bool Plus(const double *x, const double *delta, double *x_plus_delta) const {
        Eigen::Map<const Eigen::Vector3d> _p(x);
        Eigen::Map<const Eigen::Quaterniond> _q(x + 3);

        Eigen::Map<const Eigen::Vector3d> dp(delta);

        Eigen::Quaterniond dq = deltaQ(Eigen::Map<const Eigen::Vector3d>(delta + 3));

        Eigen::Map<Eigen::Vector3d> p(x_plus_delta);
        Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3);

        p = _p +_q.toRotationMatrix()* dp;
        q = (_q * dq).normalized();

        return true;
    }
    bool ComputeJacobian(const double *x, double *jacobian) const{
        Eigen::Map<Eigen::Matrix<double, 7, 6,Eigen::RowMajor>> j(jacobian);
        j.topRows<6>().setIdentity();
        j.bottomRows<1>().setZero();
        return true;
    }
    virtual int GlobalSize() const { return 7; };
    virtual int LocalSize() const { return 6; };
};

class PoseLocalParameterizationMath:public ceres::LocalParameterization{

public:
    bool Plus(const double *x, const double *delta, double *x_plus_delta)const {
        Eigen::Map<const Eigen::Vector3d> _p(x);
        Eigen::Map<const Eigen::Quaterniond> _q(x + 3);

        Eigen::Map<const Eigen::Vector3d> dp(delta);

        Eigen::Quaterniond dq = deltaQ(Eigen::Map<const Eigen::Vector3d>(delta + 3));

        Eigen::Map<Eigen::Vector3d> p(x_plus_delta);
        Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3);

        p = _p + _q.toRotationMatrix()*dp;
        q = (_q * dq).normalized();

        return true;
    }
    bool ComputeJacobian(const double *x, double *jacobian) const{
        Eigen::Map<Eigen::Matrix<double, 7, 6,Eigen::RowMajor>> j(jacobian);
        j.setZero();
        Eigen::Map< Eigen::Quaterniond> q(const_cast<double *>(x+3));
        Eigen::Matrix<double,4,3> Rxyz=Eigen::Matrix<double,4,3>::Zero();
        Rxyz.block<3,3>(1,0)=0.5*Eigen::Matrix3d::Identity();

        Eigen::Matrix<double,4,3> Jphi=Left(q)*Rxyz;
        Eigen::Matrix<double,1,3> Jphiw=Jphi.row(0);
        Jphi.block<3,3>(0,0)=Jphi.block<3,3>(1,0);
        Jphi.bottomRows(1)=Jphiw;
        j.block<3,3>(0,0)=q.toRotationMatrix();
        j.block<3,3>(0,3)=Eigen::Matrix3d::Zero();
        j.block<4,3>(3,0)=Eigen::Matrix<double,4,3>::Zero();
        j.block<4,3>(3,3)=Jphi;
        //std::cout<<"jac is "<<std::endl<<j<<std::endl;
        return true;
    }
    virtual int GlobalSize() const { return 7; };
    virtual int LocalSize() const { return 6; };
};
#endif //SFM_POSELOCALPARAMETER_H

//
// Created by unicorn on 2020/8/14.
//

#ifndef SFM_POSEGRAPH2D_FACTOR_H
#define SFM_POSEGRAPH2D_FACTOR_H

#include <iostream>
#include "readg2o.h"
#include <ceres/ceres.h>
#include "quaternion_utils.h"

double Normalize(double theta){
    return theta-(2.0*M_PI)*floor((theta+M_PI)/(2.0*M_PI));
}

template <typename T>
inline T NormalizeAngle(const T& angle_radians) {
    // Use ceres::floor because it is specialized for double and Jet types.
    T two_pi(2.0 * M_PI);
    return angle_radians -two_pi * ceres::floor((angle_radians + T(M_PI)) / two_pi);
}

class AngleLocalParamter:public ceres::LocalParameterization{
    virtual bool Plus(const double* x,
                      const double* delta,
                      double* x_plus_delta) const{
        x_plus_delta[0]=x[0]+delta[0];
        x_plus_delta[1]=x[1]+delta[1];
        x_plus_delta[2]=Normalize(x[2]+delta[2]);
        return true;
    }
    virtual bool ComputeJacobian(const double* x, double* jacobian) const{
        Eigen::Map<Eigen::Matrix3d,Eigen::RowMajor> j(jacobian);
        j.setIdentity();
        return true;
    }
    virtual int GlobalSize() const{
        return 3;
    };

    virtual int LocalSize() const{
        return 3;
    }

};

class AngleLocalParameterization {
public:

    template <typename T>
    bool operator()(const T* theta_radians, const T* delta_theta_radians,
                    T* theta_radians_plus_delta) const {
        theta_radians_plus_delta[0] =theta_radians[0] + delta_theta_radians[0];
        theta_radians_plus_delta[1] =theta_radians[1] + delta_theta_radians[1];
        theta_radians_plus_delta[2] =
                NormalizeAngle(theta_radians[2] + delta_theta_radians[2]);

        return true;
    }

    static ceres::LocalParameterization* Create() {
        return (new ceres::AutoDiffLocalParameterization<AngleLocalParameterization,
                3, 3>);
    }
};


class PoseGraph2D:public ceres::SizedCostFunction<3,3,3>{
public:
    PoseGraph2D(double pabx,double paby,double rab){
        pabx_=pabx;
        paby_=paby;
        rab_=rab;
    }
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const override {
        double pax=parameters[0][0],pay=parameters[0][1],ra=parameters[0][2];
        double pbx=parameters[1][0],pby=parameters[1][1],rb=parameters[1][2];
        double deltax=pbx-pax,deltay=pby-pay;
        residuals[0]=cos(ra)*deltax+sin(ra)*deltay-pabx_;
        residuals[1]=-sin(ra)*deltax+cos(ra)*deltay-paby_;
        residuals[2]=Normalize(rb-ra-rab_);
       // std::cout<<"my jac "<<std::endl;
        if(jacobians){
            if(jacobians[0]){
                Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> jac1(jacobians[0]);
                jac1(0,0)=-cos(ra);jac1(0,1)=-sin(ra);jac1(0,2)=-sin(ra)*deltax+cos(ra)*deltay;
                jac1(1,0)=sin(ra); jac1(1,1)=-cos(ra);jac1(1,2)=-cos(ra)*deltax-sin(ra)*deltay;
                jac1(2,0)=0;       jac1(2,1)=0;       jac1(2,2)=-1.0;

            }
            if(jacobians[1]){
                Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> jac2(jacobians[1]);
                jac2(0,0)=cos(ra); jac2(0,1)=sin(ra); jac2(0,2)=0;
                jac2(1,0)=-sin(ra);jac2(1,1)=cos(ra); jac2(1,2)=0;
                jac2(2,0)=0;       jac2(2,1)=0;       jac2(2,2)=1.0;

//                std::cout<<jac2<<std::endl;
            }
        }
        return true;

    }

    double pabx_,paby_,rab_;
};



class PoseGraph2DAutoDiff{
public:
    PoseGraph2DAutoDiff(double pabx,double paby,double rab){
        pabx_=pabx;
        paby_=paby;
        rab_=rab;
    }
    template <typename T>
    bool operator()( const T* const parameters1,const T* const parameters2,
                          T* residuals) const  {
        T pax=parameters1[0],pay=parameters1[1],ra=parameters1[2];
        T pbx=parameters2[0],pby=parameters2[1],rb=parameters2[2];
        T deltax=pbx-pax,deltay=pby-pay;
        residuals[0]=ceres::cos(ra)*deltax+ceres::sin(ra)*deltay-pabx_;
        residuals[1]=-ceres::sin(ra)*deltax+ceres::cos(ra)*deltay-paby_;
        residuals[2]=NormalizeAngle(rb-ra-rab_);

        return true;

    }

    static ceres::CostFunction* Create(double pabx,double paby,double rab) {
        return (new ceres::NumericDiffCostFunction<PoseGraph2DAutoDiff,ceres::CENTRAL, 3, 3, 3>(new PoseGraph2DAutoDiff(pabx,paby,rab)));
    }

    double pabx_,paby_,rab_;
};
#endif //SFM_POSEGRAPH2D_FACTOR_H

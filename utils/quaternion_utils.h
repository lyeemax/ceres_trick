//
// Created by unicorn on 2020/8/12.
//

#ifndef SFM_QUATERNION_UTILS_H
#define SFM_QUATERNION_UTILS_H

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <ceres/ceres.h>


Eigen::Matrix3d skew(Eigen::Vector3d v){
    Eigen::Matrix3d Vhat;
    Vhat<<0,-v[2],v[1],
        v[2],0,-v[0],
        -v[1],v[0],0;
    return Vhat;
}

Eigen::Matrix4d Left(Eigen::Quaterniond q){
    double qw=q.w();
    Eigen::Vector3d qv(q.x(),q.y(),q.z());
    Eigen::Matrix4d L=Eigen::Matrix4d::Zero();
    L+=qw*Eigen::Matrix4d::Identity();
    L.block<1,3>(0,1)+=-qv.transpose();
    L.block<3,1>(1,0)+=qv;
    L.block<3,3>(1,1)+=skew(qv);
    return L;
}

Eigen::Matrix4d Right(Eigen::Quaterniond q){
    double qw=q.w();
    Eigen::Vector3d qv(q.x(),q.y(),q.z());
    Eigen::Matrix4d R=Eigen::Matrix4d::Zero();
    R+=qw*Eigen::Matrix4d::Identity();
    R.block<1,3>(0,1)+=-qv.transpose();
    R.block<3,1>(1,0)+=qv;
    R.block<3,3>(1,1)+=-skew(qv);
    return R;

}
Eigen::Quaterniond deltaQ(Eigen::Vector3d theta){
    const double norm = theta.norm();
    double sinDelta= sin(norm / 2.0) / norm;
    Eigen::Quaterniond dq(cos(norm / 2.0), sinDelta * theta[0], sinDelta *theta[1],
                              sinDelta * theta[2]);

    dq.normalize();
    return dq;
}

//Eigen::Matrix<double,3,4> QuaternionsDerviate(Eigen::Quaterniond q,Eigen::Vector3d p){
//    Eigen::Matrix<double,3,4> jac;
//    double w=q.w();
//    Eigen::Vector3d qv(q.x(),q.y(),q.z());
//    jac.leftCols(1)=2*w*p+qv.cross(qv);
//    jac.rightCols(3)=2*(qv.transpose()*p*Eigen::Matrix3d::Identity()+qv*p.transpose()-p*qv.transpose()-w*skew(p));
//    return jac;
//}


#endif //SFM_QUATERNION_UTILS_H

//
// Created by unicorn on 2020/8/16.
//

#ifndef CERES_DEBUG_TRICK_SO3_UTILS_H
#define CERES_DEBUG_TRICK_SO3_UTILS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "quaternion_utils.h"

Eigen::Matrix3d toRotationMatrix(Eigen::Vector3d rvec){
    Eigen::Matrix3d R;
    double theta=rvec.norm();
    if(abs(theta)<1e-6) return Eigen::Matrix3d::Identity();
    Eigen::Vector3d v=rvec/theta;
    R=cos(theta)*Eigen::Matrix3d::Identity()+(1.0-cos(theta))*v*v.transpose()+sin(theta)*skew(v);
    return R;
}
Eigen::Vector3d toRotationVec(Eigen::Matrix3d R){
    Eigen::AngleAxisd rvec;
    rvec.fromRotationMatrix(R);
    Eigen::Vector3d v;
    double a=rvec.angle();
    v=a*rvec.axis();
    return v;
}

Eigen::Vector3d toRotationVec(Eigen::Quaterniond q){
    Eigen::Matrix3d R=q.toRotationMatrix();
    Eigen::AngleAxisd rvec;
    rvec.fromRotationMatrix(R);
    Eigen::Vector3d v;
    double a=rvec.angle();
    v=a*rvec.axis();
    return v;
}

Eigen::Quaterniond toQua(Eigen::Vector3d v){
    auto R=toRotationMatrix(v);
    return Eigen::Quaterniond(R);
}

Eigen::Matrix3d leftJacobianInverse(Eigen::Vector3d rvec){
    Eigen::Matrix3d J;
    double theta=rvec.norm();
    Eigen::Vector3d v=rvec/theta;
    double halftheta=theta/2.0;
    J=halftheta*(1.0/tan(halftheta))*Eigen::Matrix3d::Identity()+(1.0-halftheta*(1.0/tan(halftheta)))*v*v.transpose()-halftheta*skew(v);
    return J;
}
Eigen::Matrix3d leftJacobianInverse(Eigen::Matrix3d R){
    auto rvec=toRotationVec(R);
    return leftJacobianInverse(rvec);
}

Eigen::Matrix3d RightJacobianInverse(Eigen::Vector3d rvec){
    Eigen::Matrix3d J;
    double theta=rvec.norm();
    Eigen::Vector3d v=rvec/theta;
    double halftheta=theta/2.0;
    J=halftheta*(1.0/tan(halftheta))*Eigen::Matrix3d::Identity()+(1.0-halftheta*(1.0/tan(halftheta)))*v*v.transpose()+halftheta*skew(v);
    return J;
}

Eigen::Matrix3d RightJacobianInverse(Eigen::Matrix3d R){
    auto rvec=toRotationVec(R);
    return RightJacobianInverse(rvec);
}
#endif //CERES_DEBUG_TRICK_SO3_UTILS_H

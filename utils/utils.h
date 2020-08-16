//
// Created by unicorn on 2020/8/5.
//

#ifndef SFM_UTILS_H
#define SFM_UTILS_H

#include <iostream>
#include <opencv2/opencv.hpp>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <ceres/ceres.h>


Eigen::Matrix3d cvToRotation(cv::Mat R){
    Eigen::Matrix3d Rot;
    Rot<<R.at<double>(0,0),R.at<double>(0,1),R.at<double>(0,2),
            R.at<double>(1,0),R.at<double>(1,1),R.at<double>(1,2),
            R.at<double>(2,0),R.at<double>(2,1),R.at<double>(2,2);
    return Rot;
}


#endif //SFM_UTILS_H

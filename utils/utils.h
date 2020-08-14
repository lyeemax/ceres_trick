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

using namespace std;
void readTUM(const string &imgFile,const string &configFile,vector<pair<double,string>> &rgbdSets,cv::Mat &K,cv::Mat &distCoef){

    ifstream ifs(imgFile);
    if(!ifs){
        cout<<"Empty file"<<endl;
        return;
    }
    string s;
    while(getline(ifs,s)){
        stringstream ss(s);
        double t1,t2;
        string s1,s2;
        ss>>t1>>s1>>t2>>s2;
        rgbdSets.emplace_back(make_pair(t1,s1));
        rgbdSets.emplace_back(make_pair(t2,s2));
    }
    ifs.close();
    const string configYaml=configFile;
    cv::FileStorage file(configYaml,CV_STORAGE_READ);
    double fx=file["Camera.fx"];
    double fy=file["Camera.fy"];
    double cx=file["Camera.cx"];
    double cy=file["Camera.cy"];
    double k1=file["Camera.k1"];
    double k2=file["Camera.k2"];
    double k3=file["Camera.k3"];
    double p1=file["Camera.p1"];
    double p2=file["Camera.p2"];
    cv::Mat_<double> t(3,3) ;
    t<<fx,0,cx,0,fy,cy,0,0,1.0;
    K=t;
    cv::Mat_<double> t1(5,1);
    t1<<k1,k2,p1,p2,k3;
    distCoef=t1;
}

Eigen::Matrix3d skew(Eigen::Vector3d v){
    Eigen::Matrix3d Vhat;
    Vhat<<0,-v[2],v[1],
        v[2],0,-v[0],
        -v[1],v[0],0;
    return Vhat;
}

Eigen::Matrix<double,3,4> QuaternionsDerviate(Eigen::Quaterniond q,Eigen::Vector3d p){
    Eigen::Matrix<double,3,4> jac;
    double w=q.w();
    Eigen::Vector3d qv(q.x(),q.y(),q.z());
    jac.leftCols(1)=2*w*p;
    jac.rightCols(3)=2*(qv.transpose()*p*Eigen::Matrix3d::Identity()+qv*p.transpose()-p*qv.transpose()-w*skew(p));
    return jac;
}

Eigen::Matrix3d cvToRotation(cv::Mat R){
    Eigen::Matrix3d Rot;
    Rot<<R.at<double>(0,0),R.at<double>(0,1),R.at<double>(0,2),
            R.at<double>(1,0),R.at<double>(1,1),R.at<double>(1,2),
            R.at<double>(2,0),R.at<double>(2,1),R.at<double>(2,2);
    return Rot;
}


#endif //SFM_UTILS_H

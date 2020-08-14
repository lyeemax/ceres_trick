#include <iostream>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/filter.h>

#include <opencv2/opencv.hpp>
#include <fstream>
#include <vector>
#include <map>
#include "../utils.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>
#include "projection_factor.h"
#include "PoseLocalParameter.h"
#include <ceres/gradient_checker.h>


using namespace std;
int main() {

    string file="/home/unicorn/dataset/rgbd_dataset_freiburg1_rpy/association.txt";
    string dir="/home/unicorn/dataset/rgbd_dataset_freiburg1_rpy/";
    string yaml="/home/unicorn/dataset/rgbd_dataset_freiburg1_rpy/TUM1.yaml";
    vector<pair<double,string>> rgbdSets;
    cv::Mat K,disCoeff;
    readTUM(file,yaml,rgbdSets,K,disCoeff);
    vector<cv::Mat> rgbImgs,depthImgs;

    for (int i = 0; i <rgbdSets.size() ; ++i) {
        cv::Mat img=cv::imread(string(dir+rgbdSets[i].second),-1);
        if(img.channels()==3){
            cv::Mat undisImg;
            //cv::undistort(img,undisImg,K,disCoeff);
            rgbImgs.emplace_back(img);
        } else{
            //cout<<img.type()<<endl;
            depthImgs.emplace_back(img);
        }
    }
    //cout<<"total read "<<rgbImgs.size()<<endl;
    assert(depthImgs.size()==rgbImgs.size());

    double cx=K.at<double>(0,2),cy=K.at<double>(1,2);
    double invfx=1.0/K.at<double>(0,0),invfy=1.0/K.at<double>(1,1);
    vector<cv::Mat> Rset;
    vector<cv::Mat> Tset;
    cv::Ptr<cv::ORB> orb=cv::ORB::create(100);
    cv::BFMatcher matcher;
    cv::Mat mask;
    vector<vector<cv::Point3d>> Pci,Pcj;
    for (int j = 1; j <10 ; ++j) {
        int i=j-1;
        vector<cv::KeyPoint> keypointsi,keypointsj;
        keypointsi.clear();
        keypointsj.clear();
        cv::Mat desi,desj;
        orb->detectAndCompute(rgbImgs[i],mask,keypointsi,desi);
        orb->detectAndCompute(rgbImgs[j],mask,keypointsj,desj);
        vector<cv::DMatch> matches;
        matches.clear();
        matcher.match(desi,desj,matches);

        cout<<"find matches "<<matches.size()<<" in i "<<keypointsi.size()<<" in j "<<keypointsj.size()<<endl;
        vector<cv::Point2d> KPsi,KPsj;
        KPsi.clear();KPsj.clear();
        for (int m = 0; m <matches.size(); ++m) {
            auto pi=keypointsi[matches[m].queryIdx];
            auto pj=keypointsj[matches[m].trainIdx];
            KPsi.push_back(pi.pt);
            KPsj.push_back(pj.pt);
        }
        vector<bool> status(matches.size(),true);
        cv::Mat F=cv::findFundamentalMat(KPsi,KPsj);
        F.convertTo(F,CV_32F);
        for (int k = 0; k <matches.size(); ++k) {
            auto pi=KPsi[k];
            auto pj=KPsj[k];
            cv::Mat_<float> piNorm(3,1);
            piNorm<<pi.x,pi.y,1.0;
            cv::Mat_<float> pjNorm(3,1);
            pjNorm<<pj.x,pj.y,1.0;
            cv::Mat err=piNorm.t()*F*pjNorm;
            double errNorm=err.rowRange(0,1).dot(err.rowRange(0,1));
            if(errNorm>5) status[k]=false;
            //cout<<"error is "<<errNorm<<endl;
        }
        auto kiit=KPsi.begin(),kjit=KPsj.begin();
        auto mit=matches.begin();
        vector<cv::Point3d> pairi,pairj;
        pairi.clear();pairj.clear();
        for (int l = 0; l <status.size() ; ++l) {
            if(!status[l]){
                KPsi.erase(kiit+l);
                KPsj.erase(kjit+l);
                matches.erase(mit+l);
            }else{
                continue;
            }
        }

        for (int l = 0; l <matches.size() ; ++l) {
            double zi=double(depthImgs[i].at<unsigned short>(KPsi[l]))/5000.0;
            auto xi=(KPsi[l].x-cx)*zi*invfx,yi=(KPsi[l].y-cy)*zi*invfy;
            cv::Point3d pci(xi,yi,zi);
            double zj=double(depthImgs[i].at<unsigned short>(KPsj[l]))/5000.0;
            auto xj=(KPsj[l].x-cx)*zj*invfx,yj=(KPsj[l].y-cy)*zj*invfy;
            cv::Point3d pcj(xj,yj,zj);
            pairi.push_back(pci);
            pairj.push_back(pcj);
        }
        Pci.push_back(pairi);
        Pcj.push_back(pairj);
//        cv::Mat outImg;
//        cv::drawMatches(rgbImgs[i],keypointsi,rgbImgs[j],keypointsj,matches,outImg);
//        cv::imshow("compare",outImg);
//        cv::waitKey(0);

        cv::Mat Rr,tr;
        if(KPsi.size()>7){
            cv::Mat E=cv::findEssentialMat(KPsi,KPsj,K);
            cv::recoverPose(E,KPsi,KPsj,K,Rr,tr);
        }else{
            Rr=cv::Mat::eye(3,3,CV_32F);
            tr=cv::Mat_<double>(3,1).zeros(cv::Size(1,3));
        }
        Rset.push_back(Rr);
        Tset.push_back(tr);
        //cout<<"tr is "<<tr<<endl;
    }

    cout<<"size "<<Rset.size()<<"  "<<Tset.size()<<"  "<<Pci.size()<<"  "<<Pcj.size()<<"  "<<endl;
    int N=Rset.size();
    double qt[N][7];
    for(int i=0;i<N;i++){
        Eigen::Quaterniond q(cvToRotation(Rset[i]));
        qt[i][6]=q.w();
        qt[i][3]=q.x();
        qt[i][4]=q.y();
        qt[i][5]=q.z();
        qt[i][0]=Tset[i].at<double>(0,0);
        qt[i][1]=Tset[i].at<double>(1,0);
        qt[i][2]=Tset[i].at<double>(2,0);
    }

    return 0;
}

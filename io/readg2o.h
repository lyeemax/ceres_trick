//
// Created by unicorn on 2020/8/11.
//

#ifndef SFM_READG2O_H
#define SFM_READG2O_H

#include <iostream>
#include <fstream>
#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

// Reads a file in the g2o filename format that describes a pose graph
// problem. The g2o format consists of two entries, vertices and constraints.
//
// In 2D, a vertex is defined as follows:
//
// VERTEX_SE2 ID x_meters y_meters yaw_radians
//
// A constraint is defined as follows:
//
// EDGE_SE2 ID_A ID_B A_x_B A_y_B A_yaw_B I_11 I_12 I_13 I_22 I_23 I_33
//
// where I_ij is the (i, j)-th entry of the information matrix for the
// measurement.
//
//
// In 3D, a vertex is defined as follows:
//
// VERTEX_SE3:QUAT ID x y z q_x q_y q_z q_w
//
// where the quaternion is in Hamilton form.
// A constraint is defined as follows:
//
// EDGE_SE3:QUAT ID_a ID_b x_ab y_ab z_ab q_x_ab q_y_ab q_z_ab q_w_ab I_11 I_12 I_13 ... I_16 I_22 I_23 ... I_26 ... I_66 // NOLINT
//
// where I_ij is the (i, j)-th entry of the information matrix for the
// measurement. Only the upper-triangular part is stored. The measurement order
// is the delta position followed by the delta orientation.
struct Vertex_SE2{
    unsigned int ID;
    Eigen::Vector2d p;
    double angle;
};

struct Vertex_SE3{
    unsigned int ID;
    Eigen::Vector3d p;
    Eigen::Quaterniond q;
    Eigen::Matrix<double,7,1> pq;
};
struct EDGE_SE2{
    unsigned int ID1;
    unsigned int ID2;
    Eigen::Vector2d t12;
    double R12;
    Eigen::Matrix3d InfoMatrix;
    EDGE_SE2(){
        InfoMatrix.setZero();
    }
};


struct EDGE_SE3{
    unsigned int ID1;
    unsigned int ID2;
    Eigen::Vector3d t12;
    Eigen::Quaterniond R12;
    Eigen::Matrix<double,6,6> InfoMatrix;
    EDGE_SE3(){
        InfoMatrix.setZero();
    }
};

typedef std::map<int,Vertex_SE2> Vertex_SE2T;
typedef std::map<int,Vertex_SE3> Vertex_SE3T;
typedef std::multimap<int,EDGE_SE3> EDGE_SE3T;
typedef std::multimap<int,EDGE_SE2> EDGE_SE2T;

bool readg2oSE2(const std::string &filename, Vertex_SE2T &Vertexs,EDGE_SE2T &Edges){
    std::string VSE2T="";
    std::string ESE2T="";
    std::fstream f(filename);
    if(!f.is_open()){
        std::cout<<"file does not exits"<<std::endl;
    }
    std::string s;
    while(std::getline(f,s)){
        std::string Type;
        Type.clear();
        Type=s.substr(0,s.find_first_of(' '));
        std::stringstream ss(s);
        if(Type!="VERTEX_SE2" && Type!="EDGE_SE2"){
            std::cout<<"File is not type of SE2"<<std::endl;
            f.close();
            return false;
        }
        if(Type=="VERTEX_SE2"){
            Vertex_SE2 v;
            std::string ns;
            ss>>ns>>v.ID>>v.p.x()>>v.p.y()>>v.angle;
            Vertexs.insert(std::make_pair(v.ID,v));
        }else if (Type=="EDGE_SE2"){
            EDGE_SE2 e;
            std::string ns;
            ss>>ns>>e.ID1>>e.ID2>>e.t12.x()>>e.t12.y()>>e.R12
            >>e.InfoMatrix(0,0)>>e.InfoMatrix(0,1)>>e.InfoMatrix(0,2)
                                >>e.InfoMatrix(1,1)>>e.InfoMatrix(1,2)
                                                    >>e.InfoMatrix(2,2);
            for (int i = 0; i <3 ; ++i) {
                for (int j = 0; j <i ; ++j) {
                    e.InfoMatrix(i,j)=e.InfoMatrix(j,i);
                }
            }
            //std::cout<<e.InfoMatrix<<std::endl;
            Edges.insert(std::make_pair(e.ID1,e));
        }
    }
    f.close();
    return true;

}


bool readg2oSE3(const std::string &filename, Vertex_SE3T &Vertexs,EDGE_SE3T &Edges){
    std::fstream f(filename);
    if(!f.is_open()){
        std::cout<<"file does not exits"<<std::endl;
    }
    std::string s;
    while(std::getline(f,s)){
        std::string Type=s.substr(0,s.find_first_of(' '));
        std::stringstream ss(s);
        if(Type!="VERTEX_SE3:QUAT" && Type!="EDGE_SE3:QUAT"){
            std::cout<<"File is not type of SE3"<<std::endl;
            f.close();
            return false;
        }
        if(Type=="VERTEX_SE3:QUAT"){
            Vertex_SE3 v;
            std::string ns;
            ss>>ns>>v.ID>>v.p.x()>>v.p.y()>>v.p.z()>>v.q.x()>>v.q.y()>>v.q.z()>>v.q.w();
            v.pq.topRows(3)=v.p;
            v.pq.bottomRows(4)=v.q.coeffs();
            Vertexs.insert(std::make_pair(v.ID,v));
        } else if (Type=="EDGE_SE3:QUAT"){
            EDGE_SE3 e;
            std::string ns;
            ss>>ns>>e.ID1>>e.ID2>>e.t12.x()>>e.t12.y()>>e.t12.z()>>e.R12.x()>>e.R12.y()>>e.R12.z()>>e.R12.w()
                    >>e.InfoMatrix(0,0)>>e.InfoMatrix(0,1)>>e.InfoMatrix(0,2)>>e.InfoMatrix(0,3)>>e.InfoMatrix(0,4)>>e.InfoMatrix(0,5)
                    >>e.InfoMatrix(1,1)>>e.InfoMatrix(1,2)>>e.InfoMatrix(1,3)>>e.InfoMatrix(1,4)>>e.InfoMatrix(1,5)
                    >>e.InfoMatrix(2,2)>>e.InfoMatrix(2,3)>>e.InfoMatrix(2,4)>>e.InfoMatrix(2,5)
                    >>e.InfoMatrix(3,3)>>e.InfoMatrix(3,4)>>e.InfoMatrix(3,5)
                    >>e.InfoMatrix(4,4)>>e.InfoMatrix(4,5)
                    >>e.InfoMatrix(5,5);

            for (int i = 0; i <6 ; ++i) {
                for (int j = 0; j <i ; ++j) {
                    e.InfoMatrix(i,j)=e.InfoMatrix(j,i);
                }
            }
            //std::cout<<e.InfoMatrix<<std::endl;


            Edges.insert(std::make_pair(e.ID1,e));
        }
    }
    f.close();
    return true;

}

bool savePose(const std::string& filename,double PQ[][7],int n) {
    std::fstream file;
    file.open(filename.c_str(),std::istream::out);
    if(!file){
        std::cout<<"can not write file"<<std::endl;
        return false;
    }
    for (int i = 0; i <n ; ++i) {
        file<<i<<" "<<PQ[i][0]<<" "<<PQ[i][1]<<" "<<PQ[i][2]<<" "<<PQ[i][3]<<" "<<PQ[i][4]<<" "<<PQ[i][5]<<" "<<PQ[i][6]<<'\n';
    }
    return true;
}

bool savePose2D(const std::string& filename,double PQ[][3],int n) {
    std::fstream file;
    file.open(filename.c_str(),std::istream::out);
    if(!file){
        std::cout<<"can not write file"<<std::endl;
        return false;
    }
    for (int i = 0; i <n ; ++i) {
        file<<i<<" "<<PQ[i][0]<<" "<<PQ[i][1]<<" "<<PQ[i][2]<<'\n';
    }
    return true;
}

#endif //SFM_READG2O_H

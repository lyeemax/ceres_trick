//
// Created by unicorn on 2020/8/12.
//

#ifndef SFM_POSEGRAPH_FACTOR_H
#define SFM_POSEGRAPH_FACTOR_H

#include <iostream>
#include "readg2o.h"
#include <ceres/ceres.h>
#include "quaternion_utils.h"

class PoseGraph:public ceres::SizedCostFunction<6,7,7>{
public:
    PoseGraph(Eigen::Vector3d tab,Eigen::Quaterniond qab,Eigen::Matrix<double,6,6> info){
        qab_=qab;
        tab_=tab;
        Infomatrix=info;
    }
     bool Evaluate(double const* const* parameters,double* residuals,double** jacobians) const{
        Eigen::Vector3d Pa(parameters[0][0],parameters[0][1],parameters[0][2]);
        Eigen::Quaterniond Qa(parameters[0][6],parameters[0][3],parameters[0][4],parameters[0][5]);
        Qa.normalize();

        Eigen::Vector3d Pb(parameters[1][0],parameters[1][1],parameters[1][2]);
        Eigen::Quaterniond Qb(parameters[1][6],parameters[1][3],parameters[1][4],parameters[1][5]);
        Qb.normalize();

        Eigen::Map<Eigen::Matrix<double,6,1>,Eigen::RowMajor> res(residuals);

        Eigen::Vector3d rt=Qa.toRotationMatrix().transpose()*(Pb-Pa)-tab_;
        Eigen::Quaterniond qn=(Qa.inverse()*Qb*qab_.inverse()).normalized();
        Eigen::Vector3d rq=2.0*qn.vec();

        res.block<3,1>(0,0)=rt;
        res.block<3,1>(3,0)=rq;
        res=Infomatrix*res;
        //std::cout<<"res is "<<res.transpose()<<std::endl;

        if(jacobians){
            if(jacobians[0]){
                Eigen::Matrix3d rtpa=-Qa.toRotationMatrix().transpose();
                Eigen::Matrix3d zero=Eigen::Matrix3d::Zero();
                Eigen::Matrix<double,3,4> Lxzy=Eigen::Matrix<double,3,4>::Zero();
                Lxzy.block<3,3>(0,1)=Eigen::Matrix3d::Identity();
                Eigen::Matrix<double,4,3> Rxyz=Eigen::Matrix<double,4,3>::Zero();
                Rxyz.block<3,3>(1,0)=0.5*Eigen::Matrix3d::Identity();
                Eigen::Matrix3d rqqa=-2.0*Lxzy*Left(qab_*Qb.inverse()*Qa)*Rxyz;
                Eigen::Matrix3d rtqa=skew(Qa.inverse().toRotationMatrix()*(Pb-Pa));
                Eigen::Map<Eigen::Matrix<double,6,7>,Eigen::ColMajor> Jac0(jacobians[0]);
                Jac0.block<3,3>(0,0)=rtpa;
                Jac0.block<3,3>(0,3)=rtqa;
                Jac0.block<3,3>(3,0)=zero;
                Jac0.block<3,3>(3,3)=rqqa;
                Jac0.rightCols<1>().setZero();
                //std::cout<<"original A Jac "<<std::endl<<Jac0<<std::endl;
            }
            if(jacobians[1]){
                Eigen::Matrix3d zero=Eigen::Matrix3d::Zero();
                Eigen::Matrix<double,3,4> Lxzy=Eigen::Matrix<double,3,4>::Zero();
                Lxzy.block<3,3>(0,1)=Eigen::Matrix3d::Identity();
                Eigen::Matrix<double,4,3> Rxyz=Eigen::Matrix<double,4,3>::Zero();
                Rxyz.block<3,3>(1,0)=0.5*Eigen::Matrix3d::Identity();
                Eigen::Matrix3d rqqb=2.0*Lxzy*Left(Qa.inverse()*Qb)*Right(qab_.inverse())*Rxyz;
                Eigen::Map<Eigen::Matrix<double,6,7>,Eigen::ColMajor> Jac1(jacobians[1]);
                Jac1.block<3,3>(0,0)=Qa.toRotationMatrix().transpose();
                Jac1.block<3,3>(0,3)=zero;
                Jac1.block<3,3>(3,0)=zero;
                Jac1.block<3,3>(3,3)=rqqb;
                Jac1.rightCols<1>().setZero();
                //std::cout<<"original B Jac "<<std::endl<<Jac1<<std::endl;

            }
        }
        return true;
    }
    bool EvaluateSO3(double const* const* parameters,
                  double* residuals,
                  double** jacobians) const{
        //#TODO
        return true;
    }
    void check(double **parameters)
    {
        double *res = new double[6];
        double **jaco = new double *[2];
        jaco[0] = new double[6 * 7];
        jaco[1] = new double[6 * 7];
        Evaluate(parameters, res, jaco);
        puts("check begins");
        Eigen::Matrix<double, 7, 6,Eigen::RowMajor> Trans;
        Trans.topRows<6>().setIdentity();
        Trans.bottomRows<1>().setZero();

        std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7>,Eigen::RowMajor>(jaco[0])*Trans << std::endl
                  << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7>,Eigen::RowMajor>(jaco[1])*Trans << std::endl
                  << std::endl;


        Eigen::Vector3d Pa(parameters[0][0],parameters[0][1],parameters[0][2]);
        Eigen::Quaterniond Qa(parameters[0][6],parameters[0][3],parameters[0][4],parameters[0][5]);
        Qa.normalize();

        Eigen::Vector3d Pb(parameters[1][0],parameters[1][1],parameters[1][2]);
        Eigen::Quaterniond Qb(parameters[1][6],parameters[1][3],parameters[1][4],parameters[1][5]);
        Qb.normalize();

        Eigen::Matrix<double,6,1> origin_res;

        Eigen::Vector3d rt=Qa.toRotationMatrix().transpose()*(Pb-Pa)-tab_;
        Eigen::Quaterniond qn=(Qa.inverse()*Qb*qab_.inverse()).normalized();
        Eigen::Vector3d rq=2.0*qn.vec();

        origin_res.block<3,1>(0,0)=rt;
        origin_res.block<3,1>(3,0)=rq;

        const double eps = 1e-6;
        Eigen::Matrix<double, 6, 14> num_jacobian;
        for (int k = 0; k < 14; k++)
        {
            Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
            Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
            Qi.normalize();

            Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
            Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
            Qj.normalize();
            int a = k / 3, b = k % 3;
            Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

            if (a == 0)
                Pi += delta;
            else if (a == 1)
                Qi = Qi * deltaQ(delta);
            else if (a == 2)
                Pj += delta;
            else if (a == 3)
                Qj = Qj * deltaQ(delta);

            Eigen::Matrix<double,6,1> tmp_residual;
            Eigen::Vector3d rtt=Qi.toRotationMatrix().transpose()*(Pj-Pi)-tab_;
            Eigen::Quaterniond qnt=(Qi.inverse()*Qj*qab_.inverse()).normalized();
            Eigen::Vector3d rqt=2.0*qnt.vec();

            tmp_residual.block<3,1>(0,0)=rtt;
            tmp_residual.block<3,1>(3,0)=rqt;
            num_jacobian.col(k) = (tmp_residual - origin_res) / eps;
        }
        std::cout <<"Numercial A "<<std::endl<< (num_jacobian.leftCols(6)) << std::endl;
        std::cout <<"Numercial B "<<std::endl<<(num_jacobian.block<6,6>(0,6)) << std::endl;
    }

    Eigen::Quaterniond qab_;
    Eigen::Vector3d tab_;
    Eigen::Matrix<double,6,6> Infomatrix;

};

class PoseGraphAutoDiff{
public:
    PoseGraphAutoDiff(Eigen::Vector3d tab,Eigen::Quaterniond qab,Eigen::Matrix<double,6,6> info){
        qab_=qab.normalized();
        tab_=tab;
        Infomatrix=info;
    }
    template <typename T>
    bool operator()(const T* const pqa,const T* const pqb,T* residual) const{
        Eigen::Map<const Eigen::Matrix<T,3,1>> Pa(pqa);
        Eigen::Map<const Eigen::Quaternion<T>> Qa(pqa+3);

        Eigen::Map<const Eigen::Matrix<T,3,1>> Pb(pqb);
        Eigen::Map<const Eigen::Quaternion<T>> Qb(pqb+3);

        Eigen::Map<Eigen::Matrix<T,6,1>> res(residual);

        Eigen::Matrix<T,3,1> rt=Qa.toRotationMatrix().transpose()*(Pb-Pa)-tab_.template cast<T>();
        Eigen::Quaternion<T> qn=(Qa.inverse()*Qb*qab_.template cast<T>().inverse()).normalized();
        Eigen::Matrix<T,3,1> rq=T(2.0)*qn.vec();

        res.topRows(3)=rt;
        res.bottomRows(3)=rq;
        res=Infomatrix.template cast<T>()*res;
        //std::cout<<"res is "<<res<<std::endl;
        return true;
    }
    static ceres::CostFunction* create( Eigen::Vector3d tab, Eigen::Quaterniond qab, Eigen::Matrix<double,6,6> info) {
        return (new ceres::AutoDiffCostFunction<PoseGraphAutoDiff, 6, 7, 7>(
                new PoseGraphAutoDiff(tab,qab,info)));
    }
    Eigen::Quaterniond qab_;
    Eigen::Vector3d tab_;
    Eigen::Matrix<double,6,6> Infomatrix;



};

#endif //SFM_POSEGRAPH_FACTOR_H

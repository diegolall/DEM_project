//
// Created by Diego Lallemand on 16/04/2024.
//

#include "Bond.h"
#include "Disk.h"
#include "math.h"
#include "iostream"

Bond::Bond(Disk& i_d_1, Disk& i_d_2) {
    disk_1=&i_d_1;
    //std::cout<< &i_d_1 << std::endl;
    disk_2=&i_d_2;
    m_theta_1=(disk_1->theta());
    m_theta_2=(disk_2->theta());

    cp_1=disk_1->r()+Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_2),0);
    cp_2=disk_2->r()+Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);

    rij=cp_2-cp_1;
    n=(disk_2->r()-disk_1->r()).normalized();
    t=(rij - (rij.dot(n))*n).normalized();

    m_un=0.;
    m_ut=0.;
    v_rel=Eigen::Vector3d(0,0,0);

    kn=100.;
    kt=100.;
    ktheta=0.05;
    dtheta=0.;
    m_M=Eigen::Vector3d::Zero();
}

Bond::~Bond(){

}

void Bond::computeBond() {
    update_bond();
    m_un=rij.dot(n);
    m_ut=rij.dot(t);
    m_M=Eigen::Vector3d(0.,0.,-dtheta*ktheta);
    Fn=-kn*m_un*n;
    Ft=-kt*m_ut*t;

    //frottement
    frot_n =-5*v_rel.dot(n)*n;
    frot_t =-5*v_rel.dot(t)*t;
    disk_1->add_force(-frot_t-frot_n);
    disk_2->add_force(frot_t+frot_n);
    //forces
    Eigen::Vector3d F= Fn + Ft;
    disk_1->add_force(-F);
    disk_2->add_force(F);

    //moment de forces
    disk_1->add_momentum(-F.cross(l1));//moments au point de contacte
    disk_2->add_momentum(F.cross(l2));
    disk_1->add_momentum(-m_M);
    disk_2->add_momentum(m_M);


    //first disk only momentum
    /*
    if(disk_1->index() == 0){
        disk_1->add_momentum(m_M);
    }
    */
}

void Bond::update_bond() {
    m_theta_1=disk_1->theta();
    m_theta_2=disk_2->theta();

    l1=Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_2),0);
    l2=Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);

    cp_1=disk_1->r()+ l1;
    cp_2=disk_2->r()+ l2;

    v_rel=disk_2->v()-disk_1->v();

    rij=cp_2-cp_1;
    n=(disk_2->r()-disk_1->r()).normalized();
    /*
    double dot = rij.dot(n);
    double det = rij.x()*n.y() - rij.y()*n.x();
    dtheta=atan2(det,dot);
     */
    dtheta=m_theta_2-m_theta_1;
    t=(rij - (rij.dot(n))*n).normalized();
}
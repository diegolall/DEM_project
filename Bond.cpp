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
    n=disk_2->r()-disk_1->r();
    n.normalize();
    t=n.cross(Eigen::Vector3d(0.,0.,1));
    t.normalize();
    m_un=0.;
    m_ut=0.;

    kn=1000.;
    kt=1000.;
    ktheta=1.;
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
    Ft=kt*m_ut*t;
    if(disk_2->v().dot(t)-disk_1->v().dot(t)>0){
        Ft=-Ft;
    }
    Eigen::Vector3d Force= Fn + Ft;
    disk_2->add_force(Force);
    disk_2->add_momentum(Force.cross(cp_2));
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

    cp_1=disk_1->r()+ Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_2),0);
    cp_2=disk_2->r()+ Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);

    rij=cp_2-cp_1;
    n=(disk_2->r()-disk_1->r()).normalized();
    /*
    double dot = rij.dot(n);
    double det = rij.x()*n.y() - rij.y()*n.x();
    dtheta=atan2(det,dot);
     */
    dtheta=m_theta_2-m_theta_1;
    t=n.cross(Eigen::Vector3d(0.,0.,1.)).normalized();
}
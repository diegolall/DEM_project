//
// Created by Diego Lallemand on 16/04/2024.
//

#include "Bond.h"
#include "Disk.h"
#include "cmath"

Bond::Bond(Disk& i_d_1, Disk& i_d_2) {
    disk_1=&i_d_1;
    disk_2=&i_d_2;
    m_theta_1=(disk_1->theta());
    m_theta_2=(disk_2->theta());

    cp_1=disk_1->r()+Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_2),0);
    cp_2=disk_2->r()+Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);

    rij=cp_2-cp_1;
    n=disk_2->r()-disk_2->r();
    n.normalize();
    t=n.cross(Eigen::Vector3d(0.,0.,1));
    t.normalize();

    kn=10000.;
    kt=100000.;
}

Bond::~Bond(){

}

void Bond::computeBond() {
    update_bond();
    m_un=rij.dot(n);
    m_ut=rij.dot(t);
    disk_2->add_force(kn*m_un*n);
    disk_2->add_force(kt*m_ut*t);
}

void Bond::update_bond() {
    m_theta_1=disk_1->theta();
    m_theta_2=disk_2->theta();

    cp_1=disk_1->r()+Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_2),0);
    cp_2=disk_2->r()+Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);

    rij=cp_2-cp_1;
    n=disk_2->r()-disk_2->r();
    t=n.cross(Eigen::Vector3d(0.,0.,1));
}
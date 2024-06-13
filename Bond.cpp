//
// Created by Diego Lallemand on 16/04/2024.
//

#include "Bond.h"
#include "Disk.h"
#include "math.h"
#include "iostream"

//20mm de longueur et densité = 0.75mg/cm3 (épaisseur +- 350 micrometre)
//young modulus = 0.3 MPa

Bond::Bond(Disk& i_d_1, Disk& i_d_2) {
    //Bond variables
    disk_1=&i_d_1;
    disk_2=&i_d_2;
    m_theta_1=disk_1->theta();
    m_theta_2=disk_2->theta();

    cp_1=disk_1->r()+Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_1),0);
    cp_2=disk_2->r()+Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);
    rij=cp_2-cp_1;

    n=(disk_2->r()-disk_1->r()).normalized();
    t=(rij - (rij.dot(n))*n).normalized();

    m_un=0.;
    m_ut=0.;
    v_rel=Eigen::Vector3d(0,0,0);

    //Bond parameter
    E=disk_1->get_E();
    A=M_PI*(disk_1->radius()*disk_1->radius());
    l_b=2.*disk_1->radius();
    e=0.38;
    G=E/(2+2*e);

    kn=E*A/l_b;
    kt=G*A/l_b;
    //kn=100;
    //kt=10;

    ktheta=E*disk_1->inertia()/l_b;
    if(disk_1->index()==0) printf("kn = %lf, kt=%lf et ktheta= %lf\n",kn,kt,ktheta);
    dtheta=0.;
    nu=100.*disk_1->mass();//frottement visceux
    m_M=Eigen::Vector3d::Zero();
}

Bond::~Bond(){

}

void Bond::computeBond() {
    this->update_bond();
    m_un=rij.dot(n);
    m_ut=rij.dot(t);
    //std::cout << "un=" << m_un << " and ut = " << m_ut << std::endl;
    m_M=Eigen::Vector3d(0.,0.,-dtheta*ktheta);
    Fn=-kn*m_un*n;
    Ft=-kt*m_ut*t;

    //frottement
    frot_n =-nu*v_rel.dot(n)*n;
    frot_t =-nu*v_rel.dot(t)*t;
    disk_1->add_force(-frot_t-frot_n);
    disk_2->add_force(frot_t+frot_n);


    //forces
    Eigen::Vector3d F= Fn + Ft;
    std::cout.precision(10);
    disk_1->add_force(-F);
    disk_2->add_force(F);
   // if(disk_1->index()==1 && disk_2->index()==2) printf("force = %.12lf\n", F.norm());

    //moment de forces
    disk_1->add_momentum(F.cross(l1));//moments au point de contact
    disk_2->add_momentum(-F.cross(l2));
    disk_1->add_momentum(-m_M);
    disk_2->add_momentum(m_M);
    //disk_1->add_momentum(-nu*disk_1->w());
    //disk_2->add_momentum(nu*disk_2->w());
}

void Bond::update_bond() {
    m_theta_1=disk_1->theta();
    m_theta_2=disk_2->theta();

    l1=Eigen::Vector3d(disk_1->radius()*std::cos(m_theta_1),disk_1->radius()*std::sin(m_theta_1),0);
    l2=Eigen::Vector3d(-disk_2->radius()*std::cos(m_theta_2),-disk_2->radius()*std::sin(m_theta_2),0);

    //std::cout << cp_1 << std::endl;
    //std::cout << cp_2 << std::endl;
    cp_1=disk_1->r() +l1;
    cp_2=disk_2->r() +l2;


    v_rel=disk_2->v()-disk_1->v();//vitesse relative

    rij=cp_2-cp_1;//probleme ici
    //std::cout << rij << std::endl;
    n=(disk_2->r()-disk_1->r()).normalized();


    dtheta=m_theta_2-m_theta_1;
    t=(rij - (rij.dot(n))*n).normalized();
}
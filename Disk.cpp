#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Disk.h"
#include "Cell.h"

using namespace std;

double const mu0=4.*M_PI*pow(10,-7);

Disk::Disk(int i_index, double i_radius, double i_mass, double i_x, double i_y, double i_vx, double i_vy, Eigen::Vector3d B)
{
    m_index = i_index;
    m_linkedDisk = nullptr;
    m_linkedCell = nullptr;
    m_radius = i_radius;
    m_mass = i_mass;
    m_inertia = 0.5*m_mass*m_radius*m_radius;
    
    m_r = Eigen::Vector3d(i_x,i_y,0.);
    m_v = Eigen::Vector3d(i_vx,i_vy,0.);
    m_a = Eigen::Vector3d::Zero();
    m_F = Eigen::Vector3d::Zero();
    
    m_theta = 0.;
    m_w = Eigen::Vector3d::Zero();
    m_alpha = Eigen::Vector3d::Zero();
    m_M = Eigen::Vector3d::Zero();

    m_mu=B*(4./3.*M_PI*pow(m_radius,3))/mu0;//magnetic dipole
    m_chi=m_mu.dot(Eigen::Vector3d(1.,0.,0.));
}

Disk::~Disk(){

}

void Disk::reset_force()
{
    m_F = Eigen::Vector3d::Zero();
    m_M = Eigen::Vector3d::Zero();
}

void Disk::add_force(const Eigen::Vector3d& i_F)
{
    m_F += i_F;
}

void Disk::add_momentum(const Eigen::Vector3d& i_M)
{
    m_M += i_M;
}

void Disk::add_gravity_force(const Eigen::Vector3d& i_g)
{
    m_F += m_mass*i_g;
}
void Disk::add_ddi_force(Eigen::Vector3d n, Disk& i_disk,Disk& j_disk){
    Disk* a = &i_disk;
    Disk* b = &j_disk;

    //compute the dipole dipole force
    Eigen::Vector3d f=(3.*mu0/(4.*M_PI*pow(n.norm(),4)))*((n.cross(a->m_mu)).cross(b->m_mu)+(n.cross(b->m_mu)).cross(a->m_mu)
            -2.*n*(a->m_mu.dot(b->m_mu))+5.*n*((n.cross(a->m_mu)).dot(n.cross(b->m_mu))));
    m_F+=f;
}

void Disk::update_position(double dt)
{
    m_r += m_v*dt;
    m_theta += m_w.z()*dt;
}

void Disk::update_velocity(double dt)
{
    m_a = m_F/m_mass;
    m_v += m_a*dt;
    
    m_alpha = m_M/m_inertia;
    m_w += m_alpha*dt;
}

void Disk::set_linked_disk(Disk* i_linkedDisk)
{
    m_linkedDisk = i_linkedDisk;
}

void Disk::set_linked_cell(Cell* i_linkedCell)
{
    m_linkedCell = i_linkedCell;
}

void Disk::print(int i_num)
{
    ostringstream fileName;
    fileName<<"output/grain"<<i_num<<".txt";
    ofstream myfile;
    myfile.open(fileName.str(),ios::app);
    myfile.precision(10);
    myfile<<m_index<<"\t"<<m_r.x()<<"\t"<<m_r.y()<<"\t"<<m_v.x()<<"\t"<<m_v.y()<<"\t"<<m_theta<<"\t"<<m_radius<<"\t"<<m_F.norm()<<endl;
    myfile.close();
}

Disk* Disk::linked_disk()
{
    return m_linkedDisk;
}

Cell* Disk::linked_cell()
{
    return m_linkedCell;
}

int Disk::index()
{
    return m_index;
}

double Disk::radius()
{
    return m_radius;
}

double Disk::mass()
{
    return m_mass;
}

double Disk::theta() {
    return m_theta;
}

Eigen::Vector3d Disk::r()
{
    return m_r;
}

Eigen::Vector3d Disk::v()
{
    return m_v;
}

Eigen::Vector3d Disk::w()
{
    return m_w;
}

bool Disk::is_touching(double i_x,double i_y, double i_radius)
{
    return (m_r - Eigen::Vector3d(i_x,i_y,0.)).norm() < (m_radius+i_radius);
}

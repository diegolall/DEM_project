#include "Functions.h"
#include "Disk.h"
#include "Plan.h"

void compute_contact(Disk& i_disk, Disk* j_disk_ptr, double i_deltan, Eigen::Vector3d& i_n, double i_kn, double i_e, double i_mu)
{
    Disk* a = &i_disk;
    Disk* b = j_disk_ptr;
    
    //contact base
    i_n.normalize();
    Eigen::Vector3d v = a->v() - b->v() - a->radius()*a->w().cross(i_n) - b->radius()*b->w().cross(i_n);
    double vn_scalar = v.dot(i_n);
    Eigen::Vector3d vt = v - vn_scalar*i_n;
    double vt_scalar = vt.norm();
    Eigen::Vector3d t = (vt_scalar > 0.)? vt.normalized():Eigen::Vector3d(0.,1.,0.);
    
    //contact forces and torque
    double effectiveMass = (a->mass()*b->mass())/(a->mass()+b->mass());
    double eta = -2.*log(i_e)*sqrt(effectiveMass*i_kn/(log(i_e)*log(i_e)+M_PI*M_PI));
    
    //forces
    double Fn = -i_kn*i_deltan-eta*vn_scalar;
    if(Fn > 0.)
    {
        a->add_force(Fn*i_n);
        b->add_force(-Fn*i_n);
    }
    else
    {
        Fn=0.;
    }
    
    double Ft = -1000000.*vt_scalar;
    if(-Ft > i_mu*Fn)
    {
        Ft = -i_mu*Fn;
    }
    a->add_force(Ft*t);
    b->add_force(-Ft*t);
    
    //torque
    Eigen::Vector3d M = -Ft*a->radius()*i_n.cross(t);
    a->add_momentum(M);
    M = -Ft*b->radius()*i_n.cross(t);
    b->add_momentum(M);
}

void compute_contact(Disk& i_disk, double i_deltan, Eigen::Vector3d& i_n, double i_kn, double i_e, double i_mu)
{
    //contact base
    i_n.normalize();
    Eigen::Vector3d v = i_disk.v() - i_disk.radius()*i_disk.w().cross(i_n);
    double vn_scalar = v.dot(i_n);
    Eigen::Vector3d vt = v - vn_scalar*i_n;
    double vt_scalar = vt.norm();
    Eigen::Vector3d t = (vt_scalar > 0.)? vt.normalized():Eigen::Vector3d(0.,1.,0.);
    
    //contact forces and torque
    double effectiveMass = i_disk.mass();
    double eta = -2.*log(i_e)*sqrt(effectiveMass*i_kn/(log(i_e)*log(i_e)+M_PI*M_PI));
    
    //forces
    double Fn = -i_kn*i_deltan-eta*vn_scalar;
    if(Fn > 0.)
    {
        i_disk.add_force(Fn*i_n);
    }
    else
    {
        Fn=0.;
    }
    
    double Ft = -1000000.*vt_scalar;
    if(-Ft > i_mu*Fn)
    {
        Ft = -i_mu*Fn;
    }
    i_disk.add_force(Ft*t);
    
    //torque
    Eigen::Vector3d M = -Ft*i_disk.radius()*i_n.cross(t);
    i_disk.add_momentum(M);
}

#include "Eigen/Dense"

class Disk;

class Bond {
private:
    Disk* disk_1;
    Disk* disk_2;
    double m_theta_1;
    double m_theta_2;

    //bond parameters
    double E; //young modulus
    double A; //cross-sectional area
    double l_b;//bond lenght
    double e;//Poisson ration
    double G;

    //contact point
    Eigen::Vector3d cp_1;
    Eigen::Vector3d cp_2;
    Eigen::Vector3d l1;
    Eigen::Vector3d l2;

    //contact basis
    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d rij;
    Eigen::Vector3d Ft;
    Eigen::Vector3d Fn;
    Eigen::Vector3d v_rel; //vitesse relative
    Eigen::Vector3d frot_t;
    Eigen::Vector3d frot_n;
    double dtheta;

    //constant
    double kn;//spring constant
    double kt;
    double ktheta;
    double nu;//frottement visceux

    double m_un;
    double m_ut;
    Eigen::Vector3d m_M;

public:
    Bond(Disk&,Disk&);
    ~Bond();
    void computeBond();
    void update_bond();
    void set_nu(double i_nu);
};

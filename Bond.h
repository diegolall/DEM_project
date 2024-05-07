#include "Eigen/Dense"

class Disk;

class Bond {
private:
    Disk* disk_1;
    Disk* disk_2;
    double m_theta_1;
    double m_theta_2;

    //contact point
    Eigen::Vector3d cp_1;
    Eigen::Vector3d cp_2;
    Eigen::Vector3d lever;

    //contact basis
    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d rij;
    Eigen::Vector3d Ft;
    Eigen::Vector3d Fn;
    double dtheta;

    //constant
    double kn;//spring constant
    double kt;
    double ktheta;


    double m_un;
    double m_ut;
    Eigen::Vector3d m_M;

public:
    Bond(Disk&,Disk&);
    ~Bond();
    void computeBond();
    void update_bond();
};

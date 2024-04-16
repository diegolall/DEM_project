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

    //contact basis
    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d rij;

    //constant
    double kn;//spring constant
    double kt;

    double m_un;
    double m_ut;

public:
    Bond(Disk&,Disk&);
    ~Bond();
    void computeBond();
    void update_bond();
};

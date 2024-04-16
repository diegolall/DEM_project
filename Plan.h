#include "Eigen/Dense"

class Plan
{
    
private:
    
    Eigen::Vector3d m_r;
    Eigen::Vector3d m_n;

public:
    Plan(double,double,double,double);
    ~Plan();
    Eigen::Vector3d r();
    Eigen::Vector3d n();
};

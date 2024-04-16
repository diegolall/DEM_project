#include <math.h>
#include "Plan.h"

Plan::Plan(double i_x, double i_y, double i_nx, double i_ny)
{
    m_r = Eigen::Vector3d(i_x,i_y,0.);
    m_n = Eigen::Vector3d(i_nx,i_ny,0.);
}

Plan::~Plan()
{
    
}

Eigen::Vector3d Plan::r()
{
    return m_r;
}

Eigen::Vector3d Plan::n()
{
    return m_n;
}

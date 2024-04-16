#include "Eigen/Dense"
class Cell;
class Disk;
class Plan;

void compute_contact(Disk&,Disk*,double,Eigen::Vector3d&,double,double,double);
void compute_contact(Disk&,double,Eigen::Vector3d&,double,double,double);

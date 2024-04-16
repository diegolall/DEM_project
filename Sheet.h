//
// Created by Diego Lallemand on 08/04/2024.
//
#include <vector>
#include "Eigen/Dense"
class Disk;

class Sheet {
private:
    //sheet
    int m_number_of_particles;
    double m_size;

    //consitutive particles
    double m_radius, m_mass;


public:
    std::vector<Disk> m_particles;
    Sheet();
    Sheet(int i_particles,double i_size,double i_mass,double lx,double ly,Eigen::Vector3d B);
    ~Sheet();

    std::vector<Disk>& get_vector();
};


//
// Created by Diego Lallemand on 08/04/2024.
//
#include <vector>
#include "Eigen/Dense"
#include "Bond.h"
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
    Sheet(int i_particles,double i_long,double i_larg,double i_epai,double i_eta,Eigen::Vector3d B,double E);
    ~Sheet();

    std::vector<Disk>& get_vector();

    void print(std::string path);
    int n_part();
};


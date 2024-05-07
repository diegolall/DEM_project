//
// Created by Diego Lallemand on 08/04/2024.
//

#include "Sheet.h"
#include <iostream>
#include <stdlib.h>
#include "Disk.h"


Sheet::Sheet(int i_particles,double i_size,double i_mass,double lx,double ly,Eigen::Vector3d B) {
    m_number_of_particles=i_particles;
    double lx_shifted = lx*0.9;
    m_size=i_size;
    m_mass = i_mass;
    m_radius=m_size/((double)m_number_of_particles*2.);
    //printf("number of disk : %d\n",m_number_of_particles);

    m_particles.reserve(m_number_of_particles);
    for(int i=0;i<m_number_of_particles;i++){
        m_particles.emplace_back(i,m_radius,m_mass/(double)m_number_of_particles,-lx_shifted/2.+(double)i*2.*m_radius,0.,0.,0.,B);
    }
    printf("Sheet created with %d particles \n",m_particles.size());
}

Sheet::Sheet() {

}

Sheet::~Sheet() {
    //printf("sheet destroyed\n");
}

std::vector<Disk>& Sheet::get_vector() {
    return m_particles;
}

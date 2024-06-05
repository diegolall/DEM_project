//
// Created by Diego Lallemand on 08/04/2024.
//

#include "Sheet.h"
#include <iostream>
#include <stdlib.h>
#include "Disk.h"


Sheet::Sheet(int i_particles,double i_size,double e,double i_mass,Eigen::Vector3d B,double E) {
    m_number_of_particles=i_particles;
    m_size=i_size;
    m_mass = i_mass;
    m_radius=e/2.;

    printf("radius of disk : %lf\n",m_radius);

    m_particles.reserve(m_number_of_particles);
    for(int i=0;i<m_number_of_particles;i++){
        m_particles.emplace_back(i,m_radius,m_mass/(double)m_number_of_particles,-m_size/2.+(double)i*2.*m_radius,0.,0.,0.,B,E);
    }
    printf("Sheet created with %d particles \n",(int)m_particles.size());
}

Sheet::Sheet() {

}

Sheet::~Sheet() {
    //printf("sheet destroyed\n");
}

std::vector<Disk>& Sheet::get_vector() {
    return m_particles;
}

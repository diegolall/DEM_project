//
// Created by Diego Lallemand on 08/04/2024.
//

#include "Sheet.h"
#include <iostream>
#include <stdlib.h>
#include "Disk.h"
#include "fstream"


Sheet::Sheet(int i_particles,double i_long,double i_larg,double i_epai,double i_eta,Eigen::Vector3d B,double E) {
    m_number_of_particles=i_particles;
    m_size=i_long;
    m_mass = i_eta*i_long*i_larg*i_epai/(double)m_number_of_particles;
    m_radius=i_long/(double)m_number_of_particles;

    printf("radius of disk : %lf\n",m_radius);
    printf("mass of disk : %lf\n",m_mass);


    m_particles.reserve(m_number_of_particles);
    for(int i=0;i<m_number_of_particles;i++){
        if(i%3==0) {
            m_particles.emplace_back(i, m_radius, m_mass, -m_size / 2. + (double) i * 2. * m_radius, 0., 0., 0., B, E);
        }
        else{
            m_particles.emplace_back(i, m_radius, m_mass, -m_size / 2. + (double) i * 2. * m_radius, 0., 0., 0., Eigen::Vector3d(0.,0.,0.), E);
        }
    }
    //printf("Sheet created with %d particles \n",(int)m_particles.size());
}

Sheet::~Sheet() {
    //printf("sheet destroyed\n");
}

std::vector<Disk>& Sheet::get_vector() {
    return m_particles;
}

void Sheet::print(std::string path){
    std::ofstream file;
    file.open(path);
    file.precision(10);

    for(Disk& dsk : m_particles){
        file << dsk.index() << "\t" << dsk.x() << "\t" << dsk.y() << "\t" << dsk.x() << "\t" << dsk.y() << "\t" << dsk.theta() << "\t" << m_radius << "\t" << dsk.F().norm() << std::endl;
    }
    file.close();
}

int Sheet::n_part() {
    return m_number_of_particles;
}


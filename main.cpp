#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>
#include "Eigen/Dense"
#include <chrono>

#include "Cell.h"
#include "Disk.h"
#include "Plan.h"
#include "Functions.h"
#include "Sheet.h"

//use ./dem N_grain B_orientation(in degrees) intensity_of_B

int main(int argc, char *argv[])
{
    std::srand((unsigned int)time(nullptr));

    //sheet init
    double eta=1000;

    //double E=1.e5; //module de young
    double longueur =1.;
    double largeur =longueur*0.1;
    double epaisseur = largeur;
    //double mass =0.001;
    /*
    for (int i=0;i<1;i++){
        sheets.emplace_back(sheet_size[i],longueur,largeur,epaisseur,eta,B,E);
        sheets.at(i).set_index(i);
    }
    */
    //std::vector<int> size = {4,6,8,10,12,16};
    //std::vector<int> size={12};

    double E=1e7;

    int sheet_size=12;
    std::ostringstream path1;
    std::ostringstream path2;
    std::ostringstream path3;
    std::ostringstream path4;

    path1 << "output/test1.txt";
    path2 << "output/test2.txt";
    path3 << "output/test3.txt";
    path4 << "output/test4.txt";
    std::ofstream file1(path1.str());
    std::ofstream file2(path2.str());
    std::ofstream file3(path3.str());
    std::ofstream file4(path4.str());


    double theta_0=0.;
    double t_phi=0.;
    //magnetic interaction
    double phi=0.;
    double phi_rad = phi * M_PI / 180.;
    double B_norm=0.1285;
    Eigen::Vector3d B = Eigen::Vector3d(B_norm * std::cos(phi_rad), B_norm * std::sin(phi_rad), 0.);
    std::cout << "B: \n" << B << std::endl;
    Sheet sheet(sheet_size, longueur, largeur, epaisseur, eta, B, E);
    std::vector<Disk> vec = sheet.get_vector();

    std::vector<Bond> bonds;
    bonds.reserve(sheet_size);
    for (int i = 0; i < sheet_size - 1; i++) {
        bonds.emplace_back(vec.at(i), vec.at(i + 1));
    }
    //video
    int fps = 100;
    double tStartCapture = 0.;

    auto currentTime = std::chrono::system_clock::now();
    std::time_t currentTime_t = std::chrono::system_clock::to_time_t(currentTime);
    std::cout << "Current time: " << std::ctime(&currentTime_t) << std::endl;

    double totalTime = 7.;
    double recTime;

    //record initial state

    for (Disk &dsk: vec) {
        dsk.print(0);
    }

    //variables
    Eigen::Vector3d gravity(0., 0., 0.);
    Eigen::Vector3d r, n;

    double radius = vec.at(0).radius();
    double I = vec.at(0).inertia();
    double ft_max = (double)  E * I / pow((longueur - radius), 2);
    //double ft_max = 10;
    double ft = ft_max;

    double dt = 0.8165 * radius * std::sqrt(eta / E) / 10.;
    int total_step = (int) (totalTime / dt);
    printf("total step: %d\n", total_step);
    /*
   std::vector<Bond> bond_vec;

   for(int i=0;i<((int)sheets.at(0).get_vector().size()-1);i++){
       bond_vec.emplace_back(sheets.at(0).get_vector().at(i),sheets.at(0).get_vector().at(i+1));
       //std::cout << "\t"<<&disk_vec.at(i)<<std::endl;
   }
    */
    int step_count = 0;
    double ft_step = ft_max * dt;

    bool nu_done=true;
    double max_t=0.;
    for (double time = 0.; time < totalTime; time += dt) {
        //grains
        if(time>totalTime/8){
            if(nu_done)
            {
                for (Bond &bdn: bonds) {
                    bdn.set_nu(0);
                }
                nu_done=false;
                max_t=time;
            }
            ft=0;

        }
        gravity = Eigen::Vector3d(0., -ft, 0.);

        for (Disk &dsk: vec) {
            //update positions
            dsk.update_position(0.5 * dt);
            //reset force
            dsk.reset_force();
            //set gravity force
            if (dsk.index() == sheet.n_part() - 1) dsk.add_gravity_force(gravity); //seulement la dernière particule
            //dsk.add_gravity_force(gravity);
        }
        //compute ddi
        /*
        for (Disk &dsk: vec) {
            for (Disk &other_disk: vec) {
                if (dsk.index() != other_disk.index()) {
                    n = dsk.r() - other_disk.r();
                    dsk.add_ddi_force(n, dsk, other_disk);
                }
            }
            dsk.add_field_force(B);
        }
        */

        //compute bond
        for (Bond &bdn: bonds) {
            bdn.computeBond();
        }
        //update velocity and position
        for (Disk &dsk: vec) {
            if (dsk.index() != 0) {
                dsk.update_velocity(dt);
                dsk.update_position(0.5 * dt);
            }
        }
        //record

        recTime = time - tStartCapture;
        if (recTime >= 0.) {
            if ((int) ((recTime + dt) * fps) > (int) (recTime * fps)) {
                for (Disk &dsk: vec) {
                    dsk.print((int) ((recTime + dt) * fps));
                }
            }
        }

        if(!nu_done) {
            double T=0.745668;
            file1 << (time - max_t)/T << "\t" << vec.at(2).y() << std::endl;
            file2 << (time - max_t)/T << "\t" << vec.at(4).y()  << std::endl;
            file3 << (time - max_t)/T << "\t" << vec.at(7).y()  << std::endl;
            file4 << (time - max_t)/T << "\t" << vec.at(11).y()  << std::endl;
        }
        //printf("step %8d/%d\n", step_count, total_step);
        step_count++;
        //printf("step: %d\n",step_count);
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();

    return 0;
}


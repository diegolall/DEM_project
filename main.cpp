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
    double eta=1000.;

    double E=1.e4; //module de young
    double longueur =0.1;
    double largeur =longueur*0.1;
    double epaisseur = largeur;
    //double mass =0.001;

    //magnetic interaction
    double phi =(argc>1)?std::stod(argv[1]):90.;
    double phi_rad=phi*M_PI/180.;
    double B_norm=(argc>2)?std::stod(argv[2]):0.1; //in tesla (0.01 telsa= 100 gauss)
    Eigen::Vector3d B=Eigen::Vector3d(B_norm*std::sin(phi_rad),B_norm*std::cos(phi_rad),0.);
    /*
    for (int i=0;i<1;i++){
        sheets.emplace_back(sheet_size[i],longueur,largeur,epaisseur,eta,B,E);
        sheets.at(i).set_index(i);
    }
    */
    Sheet sheet(6,longueur,largeur,epaisseur,eta,B,E);
    std::vector<Disk> vec = sheet.get_vector();

    std::vector<Bond> bonds;
    bonds.reserve(5);
    for(int i=0;i<5;i++){
        bonds.emplace_back(vec.at(i),vec.at(i+1));
    }
    //video
    int fps = 100;
    double tStartCapture = 0.;

    auto currentTime = std::chrono::system_clock::now();
    std::time_t currentTime_t = std::chrono::system_clock::to_time_t(currentTime);
    std::cout << "Current time: " << std::ctime(&currentTime_t) << std::endl;

    double totalTime = 10.;
    double recTime;

    //record initial state

    for(Disk& dsk : vec)
    {
        dsk.print(0);
    }
    
    //variables
    Eigen::Vector3d gravity(0.,0.,0.);
    Eigen::Vector3d r,n;
    std::ostringstream path;
    path << "output/ft_" << 6 <<".txt";
    std::ofstream file(path.str());
    for (double a = 0; a < 1; a++) {
        //double I = vec.at(0).inertia();
        double radius = vec.at(0).radius();
        //double ft_max = (double) a * E * I / pow((longueur - radius), 2);
        double ft_max=1.;
        double ft_step = 4 * ft_max /totalTime;
        double ft = 0.;

        double dt = 0.8165*radius*std::sqrt(eta/E)/100.;
        int total_step=(int)(totalTime/dt);
        printf("total step: %d\n",total_step);
            /*
           std::vector<Bond> bond_vec;

           for(int i=0;i<((int)sheets.at(0).get_vector().size()-1);i++){
               bond_vec.emplace_back(sheets.at(0).get_vector().at(i),sheets.at(0).get_vector().at(i+1));
               //std::cout << "\t"<<&disk_vec.at(i)<<std::endl;
           }
            */
        int step_count = 0;

        for (double time = 0.; time < totalTime; time += dt) {
            gravity = Eigen::Vector3d(0., -9.81, 0.);
            //grains
            for (Disk &dsk: vec) {
                //update positions
                dsk.update_position(0.5 * dt);
                    //reset force
                dsk.reset_force();
                    //set gravity force
                    //if (dsk.index() == sht.n_part() - 1) dsk.add_gravity_force(gravity); //seulement la dernière particule
                dsk.add_gravity_force(gravity);
            }

            if(ft < ft_max) ft += ft_step;

                //compute ddi
                /*
                for(Disk& dsk: disk_vec){
                    for(Disk& other_disk: disk_vec){
                        if(dsk.index() != other_disk.index()){
                            n=dsk.r()-other_disk.r();
                            dsk.add_ddi_force(n,dsk,other_disk);
                        }
                    }
                }
                */

                //compute bond
            for (Bond &bdn:bonds) {
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
            recTime=time-tStartCapture;
            if(recTime>=0.)
            {
                if((int)((recTime+dt)*fps)>(int)(recTime*fps))
                {
                    for(Disk& dsk : vec)
                    {
                        dsk.print((int)((recTime+dt)*fps));
                    }
                }
            }

            /*
            if (step_count == total_step) {
                file << ft_max*pow(longueur-radius,2)/E*I << "\t" << vec.at(vec.n_part()-1).y() << std::endl;
            }
            */

                //printf("step %8d/%d\n", step_count, total_step);
            if (step_count % (total_step / 10) == 0) {
                std::cout << "progressing : " << (int) ((double) step_count / total_step * 100.) << "%"
                << std::endl;
            }
            step_count++;
                //printf("step: %d\n",step_count);
        }
    }
    file.close();
    return 0;
}


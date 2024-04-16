#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "Eigen/Dense"
#include <chrono>

#include "Cell.h"
#include "Disk.h"
#include "Plan.h"
#include "Functions.h"
#include "Sheet.h"
#include "Bond.h"

//use ./dem N_grain B_orientation(in degrees) intensity_of_B

int main(int argc, char *argv[])
{
    std::srand((unsigned int)time(NULL));

    double dt = 0.000001;
    
    //video
    int fps = 100;
    double tStartCapture = 0.;

    auto currentTime = std::chrono::system_clock::now();
    std::time_t currentTime_t = std::chrono::system_clock::to_time_t(currentTime);
    std::cout << "Current time: " << std::ctime(&currentTime_t) << std::endl;


    double totalTime = 0.5;
    double recTime;
    int total_step=(int)(totalTime/dt);
    printf("total step: %d\n",total_step);
    
    //container
    double lx = 0.1;
    double ly = 0.1;
    //magnetic interaction
    double phi =(argc>1)?std::stod(argv[1]):0.;
    double phi_rad=phi*M_PI/180.;
    double B_norm=(argc>2)?std::stod(argv[2]):0.01; //in tesla (0.01 telsa= 100 gauss)
    Eigen::Vector3d B=Eigen::Vector3d(B_norm*std::sin(phi_rad),B_norm*std::cos(phi_rad),0.);
    
    //sheet init
    double size =0.1;
    double mass = 0.01;
    Sheet sheet=Sheet(10,size,mass,lx,ly,B);
    std::vector<Disk> &disk_vec = sheet.get_vector();

    //bond init
    std::vector<Bond> bond_vec;


    //record initial state

    for(Disk& dsk : disk_vec)
    {
        dsk.print(0);
    }
    
    //variables
    Eigen::Vector3d gravity(0.,-9.81,0.);
    Eigen::Vector3d r,n;

    //ETA variables
    int step_count=0;
    
    for(double time = 0. ; time < totalTime ; time += dt)
    {
        //grains
        for(Disk& dsk : disk_vec)
        {
            //update positions
            dsk.update_position(0.5*dt);
            
            //reset force
            dsk.reset_force();
            
            //set gravity force
            //dsk.add_gravity_force(gravity);
        }

        //compute ddi
        for(Disk& dsk: disk_vec){
            for(Disk& other_disk: disk_vec){
                if(dsk.index() != other_disk.index()){
                    n=dsk.r()-other_disk.r();
                    //dsk.add_ddi_force(n,dsk,other_disk);
                }
            }
        }
        
        //update velocity and position
        for(Disk& dsk : disk_vec)
        {
            if(dsk.index()!=0) {
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
                for(Disk& dsk : disk_vec)
                {
                    dsk.print((int)((recTime+dt)*fps));
                }
            }
        }
        //printf("step %8d/%d\n", step_count, total_step);
        if(step_count%(total_step/10)==0) {
            std::cout << "progressing : " << (double)step_count/total_step*100. <<"%" << std::endl;
        }
        step_count++;
        //printf("step: %d\n",step_count);
    }
    return 0;
}


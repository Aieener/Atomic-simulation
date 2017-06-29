// main.cpp
// HW4 Molecular Dynamics
// Main function
// Author: Yuding Ai
// Penn ID: 31295008
// Data: 2017.03.21

#include "2dblock.h"

int main(){
    double start = clock();
    block bloc1(20,20,400);

    //---------------------------------------------------
    // Minimize the Potential energy by MD
    //---------------------------------------------------
    // bloc1.MD_mini(100000); //for delta_t = 1e-15
    bloc1.MD_mini(50000); //for delta_t = 1e-14
    bloc1.out_config("MD_minimization.txt");
    bloc1.rdf("rdf_MDmini.txt");

    // //---------------------------------------------------
    // // Simulate the temperature dependence by MD
    // //---------------------------------------------------
    // //-------- T_0 = 5K----------------------------------
    // bloc1.set_config("mini/MD_minimization.txt");
    // bloc1.MD(1000,5,0.0921569,"5K");
    // bloc1.out_config("5K/MD_minimization.txt");
    // bloc1.rdf("5K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 10K---------------------------------
    // bloc1.set_config("5K/MD_minimization.txt");
    // bloc1.MD(1000,10,5,"10K");
    // bloc1.out_config("10K/MD_minimization.txt");
    // bloc1.rdf("10K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 15K---------------------------------
    // bloc1.set_config("10K/MD_minimization.txt");
    // bloc1.MD(1000,15,10,"15K");
    //
    // bloc1.out_config("15K/MD_minimization.txt");
    // bloc1.rdf("15K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 20K---------------------------------
    // bloc1.set_config("15K/MD_minimization.txt");
    // bloc1.MD(1000,20,15,"20K");
    //
    // bloc1.out_config("20K/MD_minimization.txt");
    // bloc1.rdf("20K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 25K---------------------------------
    // bloc1.set_config("20K/MD_minimization.txt");
    // bloc1.MD(1000,25,20,"25K");
    //
    // bloc1.out_config("25K/MD_minimization.txt");
    // bloc1.rdf("25K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 30K---------------------------------
    // bloc1.set_config("25K/MD_minimization.txt");
    // bloc1.MD(1000,30,25,"30K");
    //
    // bloc1.out_config("30K/MD_minimization.txt");
    // bloc1.rdf("30K/rdf_MDmini.txt");
    //
    // //from now on, the phase is starting to change
    //
    // //-------- T_0 = 35K---------------------------------
    // bloc1.set_config("30K/MD_minimization.txt");
    // bloc1.MD(1000,35,30,"35K");
    //
    // bloc1.out_config("35K/MD_minimization.txt");
    // bloc1.rdf("35K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 40K---------------------------------
    // bloc1.set_config("35K/MD_minimization.txt");
    // bloc1.MD(1000,40,35,"40K");
    //
    // bloc1.out_config("40K/MD_minimization.txt");
    // bloc1.rdf("40/rdf_MDmini.txt");
    //
    // // //-------- T_0 = 45K---------------------------------
    // bloc1.set_config("40K/MD_minimization.txt");
    // bloc1.MD(1000,45,40,"45K");
    //
    // bloc1.out_config("45K/MD_minimization.txt");
    // bloc1.rdf("45K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 50K---------------------------------
    // bloc1.set_config("45K/MD_minimization.txt");
    // bloc1.MD(1000,50,45,"50K");
    //
    // bloc1.out_config("50K/MD_minimization.txt");
    // bloc1.rdf("50K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 55K---------------------------------
    // bloc1.set_config("50K/MD_minimization.txt");
    // bloc1.MD(1000,55,50,"55K");
    //
    // bloc1.out_config("55K/MD_minimization.txt");
    // bloc1.rdf("55K/rdf_MDmini.txt");
    // //-------- T_0 = 60K---------------------------------
    //
    // bloc1.set_config("55K/MD_minimization.txt");
    // bloc1.MD(1000,60,55,"60K");
    //
    // bloc1.out_config("60K/MD_minimization.txt");
    // bloc1.rdf("60K/rdf_MDmini.txt");
    //
    // // //-------- T_0 = 65K---------------------------------
    // bloc1.set_config("60K/MD_minimization.txt");
    // bloc1.MD(1000,65,60,"65K");
    //
    // bloc1.out_config("65K/MD_minimization.txt");
    // bloc1.rdf("65K/rdf_MDmini.txt");
    //
    // // //-------- T_0 = 70K---------------------------------
    // bloc1.set_config("65K/MD_minimization.txt");
    // bloc1.MD(1000,70,65,"70K");
    //
    // bloc1.out_config("70K/MD_minimization.txt");
    // bloc1.rdf("70K/rdf_MDmini.txt");
    //
    // // //-------- T_0 = 75K---------------------------------
    //  bloc1.set_config("70K/MD_minimization.txt");
    // bloc1.MD(1000,75,70,"75K");
    //
    // bloc1.out_config("75K/MD_minimization.txt");
    // bloc1.rdf("75K/rdf_MDmini.txt");
    //
    // //-------- T_0 = 80K---------------------------------
    // bloc1.set_config("75K/MD_minimization.txt");
    // bloc1.MD(1000,80,75,"80K");
    //
    // bloc1.out_config("80K/MD_minimization.txt");
    // bloc1.rdf("80K/rdf_MDmini.txt");
    //
    // //---------------auto correlation--------------------------
    // // takes 4345 secs to run auto correlation
    // // -------------------------------------------------------
    // bloc1.set_config("5K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"5K");
    //
    // bloc1.set_config("10K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"10K");
    //
    // bloc1.set_config("15K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"15K");
    //
    // bloc1.set_config("20K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"20K");
    //
    // bloc1.set_config("25K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"25K");
    //
    // bloc1.set_config("30K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"30K");
    //
    // bloc1.set_config("35K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"35K");
    //
    // bloc1.set_config("40K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"40K");
    //
    // bloc1.set_config("45K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"45K");
    //
    // bloc1.set_config("50K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"50K");
    //
    // bloc1.set_config("55K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"55K");
    //
    // bloc1.set_config("60K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"60K");
    //
    // bloc1.set_config("65K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"65K");
    //
    // bloc1.set_config("70K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"70K");
    //
    // bloc1.set_config("75K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"75K");
    //
    // bloc1.set_config("80K/MD_minimization.txt");
    // bloc1.autocorrelation(20000,"80K");
    //
    // // ------------------diffusion ----------------
    // stringstream s;
    // bloc1.set_config("50K/MD_minimization.txt");
    // s<<1.0/50<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // bloc1.set_config("55K/MD_minimization.txt");
    // s<<1.0/55<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // bloc1.set_config("60K/MD_minimization.txt");
    // s<<1.0/60<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // bloc1.set_config("65K/MD_minimization.txt");
    // s<<1.0/65<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // bloc1.set_config("70K/MD_minimization.txt");
    // s<<1.0/70<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // bloc1.set_config("75K/MD_minimization.txt");
    // s<<1.0/75<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // bloc1.set_config("80K/MD_minimization.txt");
    // s<<1.0/80<<" "<<bloc1.diffusion(20000)<<endl;
    //
    // string filename = "selfdiffusion_v2.txt";
    // ofstream myfile(filename);
    // string data = s.str();
    // myfile<< data;
    // myfile.close();
/*  */

    double end = clock();
    cout <<"This simulation takes: "<< (double(end-start)/CLOCKS_PER_SEC)<<" sec."<<endl;
}

// main.cpp
// HW5 Monte Carlo Simulation
// Main function
// Author: Yuding Ai
// Penn ID: 31295008
// Data: 2017.04.11

#include "2dblock.h"

int main(){
    double start = clock();
    block bloc1(20,20,400);

    //---------T = 5K ----------------------
    bloc1.MonteCarlo(5E5, 5, "5K");
    bloc1.out_config("5K/minimization.txt");
    bloc1.rdf("5K/rdf_MDmini.txt");

    //----------T = 10K --------------------
    bloc1.set_config("5K/minimization.txt");
    bloc1.MonteCarlo(5E5, 10, "10K");
    bloc1.out_config("10K/minimization.txt");
    bloc1.rdf("10K/rdf_MDmini.txt");

    //----------T = 15K --------------------
    bloc1.set_config("10K/minimization.txt");
    bloc1.MonteCarlo(5E5, 15, "15K");
    bloc1.out_config("15K/minimization.txt");
    bloc1.rdf("15K/rdf_MDmini.txt");

    //----------T = 20K --------------------
    bloc1.set_config("15K/minimization.txt");
    bloc1.MonteCarlo(5E5, 20, "20K");
    bloc1.out_config("20K/minimization.txt");
    bloc1.rdf("20K/rdf_MDmini.txt");

    //----------T = 25K --------------------
    bloc1.set_config("20K/minimization.txt");
    bloc1.MonteCarlo(5E5, 25, "25K");
    bloc1.out_config("25K/minimization.txt");
    bloc1.rdf("25K/rdf_MDmini.txt");

    //----------T = 30K --------------------
    bloc1.set_config("25K/minimization.txt");
    bloc1.MonteCarlo(5E5, 30, "30K");
    bloc1.out_config("30K/minimization.txt");
    bloc1.rdf("30K/rdf_MDmini.txt");

    //----------T = 35K --------------------
    bloc1.set_config("30K/minimization.txt");
    bloc1.MonteCarlo(5E5, 35, "35K");
    bloc1.out_config("35K/minimization.txt");
    bloc1.rdf("35K/rdf_MDmini.txt");

    //----------T = 40K --------------------
    bloc1.set_config("35K/minimization.txt");
    bloc1.MonteCarlo(5E5, 40, "40K");
    bloc1.out_config("40K/minimization.txt");
    bloc1.rdf("40K/rdf_MDmini.txt");

    //----------T = 45K --------------------
    bloc1.set_config("40K/minimization.txt");
    bloc1.MonteCarlo(5E5, 45, "45K");
    bloc1.out_config("45K/minimization.txt");
    bloc1.rdf("45K/rdf_MDmini.txt");

    //----------T = 50K --------------------
    bloc1.set_config("45K/minimization.txt");
    bloc1.MonteCarlo(5E5, 50, "50K");
    bloc1.out_config("50K/minimization.txt");
    bloc1.rdf("50K/rdf_MDmini.txt");

    //----------T = 55K --------------------
    bloc1.set_config("50K/minimization.txt");
    bloc1.MonteCarlo(5E5, 55, "55K");
    bloc1.out_config("55K/minimization.txt");
    bloc1.rdf("55K/rdf_MDmini.txt");

    //----------T = 60K --------------------
    bloc1.set_config("55K/minimization.txt");
    bloc1.MonteCarlo(5E5, 60, "60K");
    bloc1.out_config("60K/minimization.txt");
    bloc1.rdf("60K/rdf_MDmini.txt");

    //----------T = 65K --------------------
    bloc1.set_config("60K/minimization.txt");
    bloc1.MonteCarlo(5E5, 65, "65K");
    bloc1.out_config("65K/minimization.txt");
    bloc1.rdf("65K/rdf_MDmini.txt");

    //----------T = 70K --------------------
    bloc1.set_config("65K/minimization.txt");
    bloc1.MonteCarlo(5E5, 70, "70K");
    bloc1.out_config("70K/minimization.txt");
    bloc1.rdf("70K/rdf_MDmini.txt");

    //----------T = 75K --------------------
    bloc1.set_config("70K/minimization.txt");
    bloc1.MonteCarlo(5E5, 75, "75K");
    bloc1.out_config("75K/minimization.txt");
    bloc1.rdf("75K/rdf_MDmini.txt");

    //----------T = 80K --------------------
    bloc1.set_config("75K/minimization.txt");
    bloc1.MonteCarlo(5E5, 80, "80K");
    bloc1.out_config("80K/minimization.txt");
    bloc1.rdf("80K/rdf_MDmini.txt");


    // --- Now collect the equalibrium info from each state----
    // --- Corresponding to each temperature ------------------

    // ------------T = 5K --------------------
    bloc1.set_config("5K/minimization.txt");
    bloc1.CollectEqdata(1E4, 5, "5K");

    // ------------T = 10K --------------------
    bloc1.set_config("10K/minimization.txt");
    bloc1.CollectEqdata(1E4, 10, "10K");
    // ------------T = 15K --------------------
    bloc1.set_config("15K/minimization.txt");
    bloc1.CollectEqdata(1E4, 15, "15K");
    // ------------T = 20K --------------------
    bloc1.set_config("20K/minimization.txt");
    bloc1.CollectEqdata(1E4, 20, "20K");
    // ------------T = 25K --------------------
    bloc1.set_config("25K/minimization.txt");
    bloc1.CollectEqdata(1E4,25, "25K");
    // ------------T = 30K --------------------
    bloc1.set_config("30K/minimization.txt");
    bloc1.CollectEqdata(1E4,30, "30K");
    // ------------T = 35K --------------------
    bloc1.set_config("35K/minimization.txt");
    bloc1.CollectEqdata(1E4, 35, "35K");
    // ------------T = 40K --------------------
    bloc1.set_config("40K/minimization.txt");
    bloc1.CollectEqdata(1E4, 40, "40K");
    // ------------T = 45K --------------------
    bloc1.set_config("45K/minimization.txt");
    bloc1.CollectEqdata(1E4, 45, "45K");
    // ------------T = 50K --------------------
    bloc1.set_config("50K/minimization.txt");
    bloc1.CollectEqdata(1E4, 50, "50K");
    // ------------T = 55K --------------------
    bloc1.set_config("55K/minimization.txt");
    bloc1.CollectEqdata(1E4, 55, "55K");
    // ------------T = 60K --------------------
    bloc1.set_config("60K/minimization.txt");
    bloc1.CollectEqdata(1E4, 60, "60K");
    // ------------T = 65K --------------------
    bloc1.set_config("65K/minimization.txt");
    bloc1.CollectEqdata(1E4, 65, "65K");
    // ------------T = 70K --------------------
    bloc1.set_config("70K/minimization.txt");
    bloc1.CollectEqdata(1E4, 70, "70K");
    // ------------T = 75K --------------------
    bloc1.set_config("75K/minimization.txt");
    bloc1.CollectEqdata(1E4,75, "75K");
    // ------------T = 80K --------------------
    bloc1.set_config("80K/minimization.txt");
    bloc1.CollectEqdata(1E4, 80, "80K");

    double end = clock();
    cout <<"This simulation takes: "<< (double(end-start)/CLOCKS_PER_SEC)<<" sec."<<endl;
    
    return 0;
}



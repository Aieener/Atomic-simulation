// main.cpp
// HW3 Molecular statics relaxation
// Main function
// Author: Yuding Ai
// Penn ID: 31295008
// Data: 2017.03.01

#include "2dblock.h"


int main() {

    block bloc1(20,20,400);
    // output the original configuration and rdf
    bloc1.rdf("rdf.txt");
    bloc1.out_config("config.txt");

    // Now examine on lambda = 0.5,1,1.5 and 2 seperately

    // ---------------lambda = 0.5 -------------
    bloc1.SD(0.5,8000,"l=05_vsN.txt");
    bloc1.out_config("l=05_config.txt");
    bloc1.rdf("l=05_rdf.txt");

    // ---------------lambda = 1 ---------------
    block bloc2(20,20,400);
    bloc2.SD(1,8000,"l=1_vsN.txt");
    bloc2.out_config("l=1_config.txt");
    bloc2.rdf("l=1_rdf.txt");

    // ---------------lambda = 1.5 -------------
    block bloc3(20,20,400);
    bloc3.SD(1.5,8000,"l=15_vsN.txt");
    bloc3.out_config("l=15_config.txt");
    bloc3.rdf("l=15_rdf.txt");

    // ---------------lambda = 2 ---------------
    block bloc4(20,20,400);
    bloc4.SD(2,8000,"l=2_vsN.txt");
    bloc4.out_config("l=2_config.txt");
    bloc4.rdf("l=2_rdf.txt");

    // ---------------lambda = 8 ---------------
    block bloc5(20,20,400);
    bloc5.SD(8,8000,"l=8_vsN.txt");
    bloc5.out_config("l=8_config.txt");
    bloc5.rdf("l=8_rdf.txt");

    // --------study the evolution -------------
    // ---------------lambda = 1.5 -------------
    block bloc6(20,20,400);

    bloc6.SD(1.5,250,"250vsN.txt");
    bloc6.out_config("250config.txt");
    bloc6.rdf("250rdf.txt");

    bloc6.SD(1.5,250,"500vsN.txt");
    bloc6.out_config("500config.txt");
    bloc6.rdf("500rdf.txt");

    bloc6.SD(1.5,500,"1000vsN.txt");
    bloc6.out_config("1000config.txt");
    bloc6.rdf("1000rdf.txt");

    bloc6.SD(1.5,1000,"2000vsN.txt");
    bloc6.out_config("2000config.txt");
    bloc6.rdf("2000rdf.txt");

    return 0;
}



// 2dblock.h	
// HW3 Molecular statics relaxation
// 2-D Block
// Author: Yuding Ai
// Date: 2017.03.01

#ifndef BLOCK_H
#define BLOCK_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <vector>
using namespace std;

/* README
 * 
 * Since I misunderstood assignment 2 and kind of get lost last 
 * time, (I thought the block is mapping into a lattice space so 
 * every atom sits on a lattice site, like a pixel, now I realize 
 * it's completely wrong.) I have revised the whole structure of
 * my 2dblock. (hopefully it's getting better this time) 
 *  
 * Here is a brief description of it:
 * 
 * A. To compile the file: g++ -std=gnu++11 main.cpp 2dblock.cpp -o run
 * B. To run the compiled file: ./run
 * C. To get the plots: python hw3.py
 * 
 * 0.Part of the implementation (the optimized algorithm for 
 * get_neighbor_list for example) is inspired by Our TA Spencer's
 * nice reference code:assignment2Demo.f90; 
 * 
 * 1. Same as before, I treat this whole 2dblock as a single object
 * with certain amount of atoms in it. In this assignment, with 400 atoms.
 * The class 2dblock have 5 attributes:
 * w: width of the block; (#of atoms each raw at initial config); 
 * h: height ot the block;(#of atoms each col at initial config); 
 * n: number of atoms.
 * neighbor_list: a 'list' of neighborlist contains neighbors for each atom. 
 * atomlist: a list of atoms, contains all the infomations we need so far
 * for each atom.
 *
 * 2. Each atom is a '7D vector' (an array with 7 components) associated
 * with a index and stored into an atomlist; such that Atom(xposition, 
 * yposition, sigmaxx, sigmayy, sigmaxy, F_x, F_y).
 * so: std::array<double,7> atom;        
 * [Remark: I might want to make class for my atom for latter 
 * assignment if needed to, but for now, I will keep it as simple as 
 * a '7D vector']
 * 
 * 3.Methods (subroutine):
 * 
 * 3.1 method get_neighbor_list(r), where r= neighCut except for 
 * calculate rdf, in rdf we set r = rdfMax. Update attribute 
 * neighbor_list into current configuration.
 * 
 * 3.2 method get_atomic_stress(), update the stress for each atom
 * and corresponding to current configuration and load it to the 
 * class arribute: atomlist.
 * 
 * 3.3 method get_distance(idx1,idx2) return the distance of 2 atom
 *     method get_xdistance(idx1,idx2) return the xdistance of 2 atom
 *     method get_ydistance(idx1,idx2) return the ydistance of 2 atom
 * 
 * 3.4 method rdf(filename),calculate the rdf and store it into a 
 * txt file named filename.
 * 
 * 3.5 method out_config(filename), output the current configuration 
 * including xpos,ypos,sigma_xx,sigma_yy,sigma_xy for each atom 
 * for into a txt file named filename for latter plotting 
 * 
 * 3.6 method total_energy(), fist update the neighbor_list and then 
 * calculate the total energy under the corresponding configuration. 
 * 
 * 3.7 method calc_force(), calculate the force for all atoms and assign 
 *     it to each atom.
 *     method calc_single_force(idx), calculate the force for a single 
 *     atom with index idx
 *     method qcalc_mforce(), quickly calculate the maximum force of the 
 *     corresponding configuration.
 *     
 * 3.8 method calc_total_stress_pressure(), calculate the net stresses 
 * of the corresponding configuration.
 * 
 * 3.9 method psi(r) compute Lenard Jones potential between 2 atoms 
 *     that are separated by r.
 *     method dpsi(r) compute derivative of Lenard Jones potential 
 *     between 2 atoms that are separated by r.
 *     
 * 4.Core molecular simulation method SD(steepest descent): 
 * SD(lambda, iteration, filename) perform a simulation with SD method and
 * output the simulation data into a txtfile named filename
 * 
 * All in all, thanks very much for reading this script and a big 
 * special thanks to Spencer's reference code which is really helpful 
 * and answers a lot of questions I had before. 
 * 
 * Again, thanks!
 */ 


///-------------------------------------------------------------
//Constants  (following the class note and Spencer's reference code)
///-------------------------------------------------------------
const double epsilon = 0.010323;
const double sigma = 3.405;
const double rMin = 3.822;
const double rTail = 7.0;
const double rCut = 7.5;
const double A = -0.0068102128;
const double B = -0.0055640876;
const double rdfMin = 1.0;
const double rdfMax = 15.0;
const double deltaR = 0.01;
const double neighCut = 7.5;
const double PI = 3.1415926;


class block
{
    private:
        int w;           // Width; 
        int h;           // Height; 
        int n;           // number of atoms;

        // each atom is represented by an array  
        array<array<double,7>,400> atomlist; 
        // Once the idx is known, all the info about idxth atom is known
        // because we could reach such atom by: block[idx]

        //a list of neighbor lists of each atom
        array<vector<int>,400> neighbor_list;

    public:
        ///------------constructor--------------------------------------
        //initialize the block with n atoms and arranged to initial
        //configuration. (the one in assignment 2)
        //And set the 3 stresses into 0 by default
        ///-------------------------------------------------------------
        block(int x,int y,int N); 


        ///-------------------------------------------------------------
        /// getter (subroutines)
        ///-------------------------------------------------------------
        
        int getN();
        
        // With this operator[], once the idx is known, all the info about 
        // idxth atom is known because we could reach such atom by: block[idx]
        // for example: its position is (block[idx][0],block[idx][1]) and its
        // stresses are (block[idx][2],block[idx][3],block[idx][4])    
        array<double,7> &operator[](int index);


        // return the distance between two atoms ----r_{ik}
        double get_distance(int idx1, int idx2) const; 

        // return the distance in x axis between two atoms ----r_{ik}^{alpha}
        double get_xdistance(int idx1, int idx2) const; 

        // return the distance in y axis between two atoms ----r_{ik}^{beta}
        double get_ydistance(int idx1, int idx2) const; 

        //------------get_neighbor_list---------------------------------
        // take a single argument r, which is neighCut in this assignment 
        // and will update neighbor_list of current configuration. 
        // Using Spencer's trick to reduce the computational complicity
        // -------------------------------------------------------------
        void get_neighbor_list(double r); 

        //------------get_atomic_stress---------------------------------
        //take no argument and compute the atomic stresses of all the 
        //atoms and assign it to each atom.
        // -------------------------------------------------------------
        void get_atomic_stress(); 

        //------------rdf() --------------------------------------------
        //Get the radial distribution and output into a txt file
        // -------------------------------------------------------------
        void rdf(string filename);


        //------------out_config()--------------------------------------
        // output the configuration of atoms for latter plotting 
        // and visualization
        // -------------------------------------------------------------
        void out_config(string filename);

        //------------total_energy()------------------------------------
        // compute the total energy of a configuration by sum over
        // the potential energy of all atoms
        // E_tot = 1/2 sum_{i,j} psi(r_{i,j})
        // -------------------------------------------------------------
        double total_energy();

        //------------calc_force() & calc_single_force()----------------
        //compute the force of a configuration for each atom by:
        //F_i = -grad_{r_i}*psi(r_{i,j}); where j are neighbors of i
        //calc_force()store the force for each atom into a list and 
        //then find and return the maximum force; 
        //
        //calc_single_force()
        //compute the net force on a single atom and assign such
        //force onto atomlist[idx][5] and atomlist[idx][6]
        //
        //qcalc_mforce() calculate the maximum force a current atomlist
        //by simply sort the force component 
        // -------------------------------------------------------------
        double calc_force();
        void calc_single_force(int idx);
        double qcalc_mforce();
        

        //------------calc_total_stress_pressure()----------------------
        //compute the total stress sigma_xx, simga_yy and simga_xy and 
        //the hydrostatic pressure P. Store the 4 values into an array
        //such that (sum_simga_xx,sum_simga_yy,sum_simga_xy,P)
        //--------------------------------------------------------------
        array<double,4> calc_total_stress_pressure();

        //------------SD(double lambda)---------------------------------
        //perform a steepest decent minimization
        //in my case, for each SD, perform 8000 iterations
        // -------------------------------------------------------------
        void SD(double lambda,int it,string filename);

};

/// useful functions

// compute Lenard Jones potential
double psi(double r);

// compute the derivative of Lenard Jones potential
double dpsi(double r);
#endif /* 2DBLOCK_H */

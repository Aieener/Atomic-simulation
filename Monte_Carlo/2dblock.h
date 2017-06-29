// 2dblock.h	
// HW5 Monte Carlo Simulation
// 2-D Block
// Author: Yuding Ai
// Date: 2017.04.11

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
 * This implementation is build on top of my assignment 4
 * As always, thanks!
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
const double deltaR = 0.1;
const double neighCut = 7.5;
const double PI = 3.1415926;

const double k_b = 8.617*1E-5; //Boltzman's constant eV*K^-1
const double mass = 39.948*1.66E-27; // the mass of argon atom kg

// const double delta_t = 1E-15; // artifitially chosen the time step for MD
const double delta_t = 1E-14; // artifitially chosen the time step for MD


class block
{
    private:
        int w;           // Width; 
        int h;           // Height; 
        int n;           // number of atoms;

        // each atom is represented by an array  
        array<array<double,9>,400> atomlist; 
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
        array<double,9> &operator[](int index);


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

        //------------total_penergy()------------------------------------
        // compute the total potential energy of a configuration by sum over
        // the potential energy of all atoms
        // E_tot = 1/2 sum_{i,j} psi(r_{i,j})
        // -------------------------------------------------------------
        double total_penergy();

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

        //===============Methods for HW 4===============================
        
        /*----------- set_config(filename) -----------------------------
         * manully set the configuration for the block,so including 
         * position, velocities, Atomic Force and so on (9 value per atom)
         */
        void set_config(string filename);

        /*-----------update_position(int idx)---------------------------
         * update the position of each atom using verlet algorithm
         * then assign the velocity into each atom
         */
        void update_position();

        /*-----------update_velosity(int idx)---------------------------
         * update the velocity of each atom using verlet algorithm
         * then assign the velocity into each atom
         */
        void update_velosity(array<double,400> F_oldx, array<double,400> F_oldy);

        /* ---------- Molecular Dynamics MD_mini() : HW part1 -----------
         * perform a molecular dynamics simulation using verlet algorithm
         * verlet algorithm  to minimize the potential energy
         */
        void MD_mini(int it);


        /* ---------- Molecular Dynamics MD() : HW part2----------------
         * perform a molecular dynamics simulation using verlet algorithm
         * verlet algorithm to simulate the temperature dependence
         */
        void MD(int it,double desire_T,double old_T,string dir);

        /* ---------- autocorrelation() ---------------------------------
         * calculate the autocorrelation function for each Temperature
         * the configuration and initial velocity is predefined before
         * we do the molecular dynamics as to calcualte the autocorreation
         * function
         */
        void autocorrelation(int it,string dir);

        /* ---------- self-diffusion coefficient() -----------------------
         * calculate the self-diffusion coefficient D 
         */

        double diffusion(int it);

        //===============Methods for HW 5===============================
        
        /* --------------------- Move() ---------------------------------
         * Attempt to move the chosen atom that
         * x_new = x_old + xi_x*alpha;
         * y_new = y_old + xi_y*alpha;
         */
        void Move(double alpha, int &sidx,double &oriXpos, double &oriYpos);

        /* ----------------- redoMove() ---------------------------------
         * in the case of rejecting movement,
         * this method take the block back to its old configuration
         */
        void redoMove(int sidx, double oriXpos, double oriYpos);
        
        /* ---------------- MonteCarlo() --------------------------------
         * the master method for Monte carlo simulations
         */
        void MonteCarlo(long int steps, double T, string dir);

        /* ------------- CollectEqdata() --------------------------------
         * Once the equlibirum state has reached, we run an additional
         * Monte carlo simulation with fewer steps as to collect
         * the average value of potential energy and hydrostatic pressure
         */
        void CollectEqdata(long int steps,double T,string dir);

};

/// useful functions

// compute Lenard Jones potential
double psi(double r);

// compute the derivative of Lenard Jones potential
double dpsi(double r);
#endif /* 2DBLOCK_H */

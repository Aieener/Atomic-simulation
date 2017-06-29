// 2dblock.cpp
// HW3 Molecular statics relaxation
// 2-D Block
// Author: Yuding Ai
// Date: 2017.03.01

#include "2dblock.h"

block::block(int x, int y,int N){
    //set with and hight to x,y 
    w = x;
    h = y;
    n = N;

	//initialize the block with n atoms and arranged to initial
    //configuration. (the one in assignment 2)
    //And set the 3 stresses into 0 by default
    for (int i = 0; i < n; i++) { 
        double xpos,ypos;
        xpos = (i%w)* rMin;  // default spacing between atoms to be rMin
        ypos = (i/h)* rMin;
        array<double,7> atom;
        atom[0] = xpos;
        atom[1] = ypos;
        atom[2] = atom[3]=atom[4] = atom[5] = atom[6]=0.0;
        atomlist[i]=atom;
    } 

    // set the initial neighbor list 
    get_neighbor_list(neighCut); 

    // calculate the initial atomic stresses
    get_atomic_stress(); 

    // calculate the initial Force on each atom 
    calc_force();
}

///------------------------------------------------------------
/// getter (subroutines)
///------------------------------------------------------------

int block::getN(){return atomlist.size();}

array<double,7> &block::operator[](int index){return atomlist[index];}


double block::get_distance(int idx1, int idx2) const {
    array <double,2> coor1 = {{atomlist[idx1][0],atomlist[idx1][1]}};
    array <double,2> coor2 = {{atomlist[idx2][0],atomlist[idx2][1]}};
    double dis;
    dis =sqrt(pow(coor1[0]-coor2[0],2) + pow(coor1[1] - coor2[1],2));
    return dis;
}

double block::get_xdistance(int idx1, int idx2) const {
    array <double,2> coor1 = {{atomlist[idx1][0],atomlist[idx1][1]}};
    array <double,2> coor2 = {{atomlist[idx2][0],atomlist[idx2][1]}};
    double dis;
    dis =coor1[0]-coor2[0];
    return dis;
}

double block::get_ydistance(int idx1, int idx2) const {
    array <double,2> coor1 = {{atomlist[idx1][0],atomlist[idx1][1]}};
    array <double,2> coor2 = {{atomlist[idx2][0],atomlist[idx2][1]}};
    double dis;
    dis =coor1[1]-coor2[1];
    return dis;
}

void block::get_neighbor_list(double r){ 
    //first remove previous neighbor list:
    
    for (int i = 0; i < n; i++) { 
        neighbor_list[i].clear(); 
    } 
    //then load the new neighbor list
    // Thanks to Spencer's nice reference code for assignment 2, here
    // we apply the same nice trick to get rid of double counting and 
    // cut computational complexity by half
     
    for (int i = 0; i < n; i++) { 
        for (int j = i+1; j < n; j++) { 
            // for this assignment, r = neighCut
            if(get_distance(i,j) <=r){
                neighbor_list[i].push_back(j); // i's new neighbor is j
                neighbor_list[j].push_back(i); // j's new neighbor is i
            }
        }
    } 
}

void block::get_atomic_stress(){ 
    const double omega = pow(rMin,2)*(w-1)*(h-1)/n;

    for(int i = 0; i<n; i++){
        //first remove the previous atomic stress and set it to 0
        atomlist[i][2] = 0;
        atomlist[i][3] = 0;
        atomlist[i][4] = 0;

        for (unsigned int  j = 0; j < neighbor_list[i].size();j++){
            //re calculate the atomic stress
            double r = get_distance(i,neighbor_list[i][j]);
            double rx = get_xdistance(i,neighbor_list[i][j]);
            double ry= get_ydistance(i,neighbor_list[i][j]);
            //sigma_xx
            atomlist[i][2] += dpsi(r)*rx*rx/(r*omega);
 
            //sigma_yy
            atomlist[i][3] += dpsi(r)*ry*ry/(r*omega);

            //sigma_xy
            atomlist[i][4] += dpsi(r)*rx*ry/(r*omega);
        }
    }
}

void block::rdf(string filename){
    stringstream st;
    const double omega = pow(rMin,2)*(w-1)*(h-1)/n;
    vector<double> dislist;
    vector<double> Rlist;
    vector<double> glist;
    get_neighbor_list(rdfMax);
    
    for (int i = 0; i < 400; i++) { 
        for (unsigned int j = 0; j < neighbor_list[i].size(); j++) { 
            dislist.push_back(get_distance(i,neighbor_list[i][j]));
        } 
    } 

    //dislist is sorted in increment order
    sort(dislist.begin(),dislist.end());

    double d = dislist[0];

    // count each distance
    int c = 0;
    for (unsigned int i = 0; i < dislist.size(); i++) { 
        if(abs(dislist[i] - d)<= deltaR){
            c++;
        }
        else{
            //rlist is a list of possible distance of all neighbors 
            //sorted in increment order
            Rlist.push_back(d);
            //glist is a list of occurrence/frequency of ith r for now
            glist.push_back(1.0*c/dislist.size()); //divide dislist to normalize it
            c = 0;
            d = dislist[i];
        }

        if(i == dislist.size()-1){
            //rlist is a list of possible distance of all neighbors 
            //sorted in increment order
            Rlist.push_back(d);
            //glist is a list of occurrence/frequency of ith r for now
            glist.push_back(1.0*c/dislist.size()); //divide dislist.size() to normalize it
            c = 0;
            d = dislist[i];
        }
    } 
    
    // last calculate g(r) and update the glist
    for (unsigned int i = 0; i < glist.size(); i++) { 
        // calculate and update glist into g(r) now, 
        // namely the radial distribution function
        glist[i] = glist[i]*omega/(2.0*PI*deltaR*Rlist[i]);
        st<< Rlist[i]<<"   "<<glist[i]<<"\n"; 
    } 

    // record rdf into a txt file for latter plotting
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();

    //last reset the neighbor list to r = neighCut
    //since we modified the r to be rdfMax every time when 
    //we calculate the rdf
    get_neighbor_list(neighCut);
}


void block::out_config(string filename){

    stringstream st;
    // record xpos, ypos and stresses of each atom into a txt file 
    // file for latter plotting
    for (int i = 0; i < n; i++) { 
        st << atomlist[i][0]<<" "<<atomlist[i][1]<<" "<<atomlist[i][2]<<" ";
        st << atomlist[i][3]<<" "<<atomlist[i][4]<<endl;
    } 
    
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();

}
double block::total_energy(){
    double E = 0;
    for (int i = 0; i < n; i++) { 
        for (unsigned int j = 0; j < neighbor_list[i].size(); j++) { 
            E += psi(get_distance(i,neighbor_list[i][j]));
        } 
    } 
    // Since E_tot = 1/2 sum_{i,j} psi(r_{i,j})
    E = E/2.0;
    return E;
}

double block::calc_force(){

    double F = 0;
    double curFx = 0;
    double curFy = 0;
    double curF = 0;
    // calc Force on each atom and find the F_max
    for (int i = 0; i < n; i++) { 
        curFx = 0;
        curFy = 0;
        for(unsigned int j = 0; j<neighbor_list[i].size();j++){
            double r_ij = get_distance(i,neighbor_list[i][j]);
            double r_ij_x = get_xdistance(i,neighbor_list[i][j]);
            double r_ij_y = get_ydistance(i,neighbor_list[i][j]);
            curFx+= -dpsi(r_ij)*(r_ij_x/r_ij);
            curFy+= -dpsi(r_ij)*(r_ij_y/r_ij);
        }
        // update the atomlist
        atomlist[i][5]=curFx;
        atomlist[i][6]=curFy;
        curF = sqrt(pow(curFx,2) + pow(curFy,2));
        if(F<curF){
            F = curF;
        }
    } 
    return F;
}

void block::calc_single_force(int idx){

    double curFx = 0;
    double curFy = 0;
    // calc Force on each atom and find the F_max
    for(unsigned int j = 0; j<neighbor_list[idx].size();j++){
        double r_ij = get_distance(idx,neighbor_list[idx][j]);
        double r_ij_x = get_xdistance(idx,neighbor_list[idx][j]);
        double r_ij_y = get_ydistance(idx,neighbor_list[idx][j]);
        curFx+= -dpsi(r_ij)*(r_ij_x/r_ij);
        curFy+= -dpsi(r_ij)*(r_ij_y/r_ij);
    }
    // update the atomlist
    atomlist[idx][5]=curFx;
    atomlist[idx][6]=curFy;
}

double block::qcalc_mforce(){
    double maxF = 0;
    double curFx=0;
    double curFy=0;
    double curF=0;
    for (int i = 0; i < n; ++i) { 
        curFx=atomlist[i][5];
        curFy=atomlist[i][6];
        curF = sqrt(pow(curFx,2) + pow(curFy,2));
        if(maxF<curF){
            maxF = curF;
        }
    } 
    return maxF;
}

array<double,4> block::calc_total_stress_pressure(){
    double sxx,syy,sxy,P;
    // P: total hydrostatic pressure:
    // for each atom i: p_i = -1/2*sum(sxx_i +syy_i)
    // P = sum{i}(p_i)
    sxx = syy = sxy = P = 0;
    for (int i = 0; i < n; ++i) { 
        sxx += atomlist[i][2]; 
        syy += atomlist[i][3]; 
        sxy += atomlist[i][4]; 
        P += -1.0/2.0*(atomlist[i][2] + atomlist[i][3]);
    } 

    array<double,4> sum_stress= {{sxx,syy,sxy,P}};
    return sum_stress;
}

void block::SD(double lambda,int it,string filename){
    stringstream st;
    int j = 0;
    int k = 0;
    while(j<it){
        for (int i = 0; i < n; i++) { 
            // get the force on each atom
            calc_single_force(i);
            // Move each atom following the direction
            // of it's force
            atomlist[i][0] += lambda*atomlist[i][5];  
            atomlist[i][1] += lambda*atomlist[i][6];
        } 

        if(it>=1000){
            if(k==it/1000){
                //first update the neighborlist
                get_neighbor_list(neighCut);
                //update the stresses
                get_atomic_stress();
                //calc hydrostatic pressure
                array<double,4> sp = calc_total_stress_pressure();
                //calc total energy
                double E = total_energy();
                //calc maximum force
                double F = qcalc_mforce();
                //record the state/configuration data
                st<<j<<"   "<<E<<"   "<<sp[0]<<"  "<<sp[1]<<"  ";
                st<<sp[2]<<"  "<<sp[3]<<"  "<<F<<endl;
                k = 0;
            }
        }
        else{
            //first update the neighborlist
            get_neighbor_list(neighCut);
            //update the stresses
            get_atomic_stress();
            //calc hydrostatic pressure
            array<double,4> sp = calc_total_stress_pressure();
            //calc total energy
            double E = total_energy();
            //calc maximum force
            double F = qcalc_mforce();
            //record the state/configuration data
            st<<j<<"   "<<E<<"   "<<sp[0]<<"  "<<sp[1]<<"  ";
            st<<sp[2]<<"  "<<sp[3]<<"  "<<F<<endl;
        }

        j++;
        k++;
    }
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();

}
///------------------------------Outside block class ------------------
double psi(double r){
    double result;
    //  r >= rCut
    if(r>=rCut){result = 0;}            
    //  rTail <= r < rCut
    else if(r>=rTail){result = A*pow(r-rCut,3) + B*pow(r-rCut,2);}  
    //  r < rTail
    else{result = 4.0*epsilon*(pow(sigma/r,12) - pow(sigma/r,6));}  
            
    return result;
} 

double dpsi(double r){
    double result;
    if(r>=rCut){result = 0;}
    else if (r>= rTail){result = 3.0*A*pow(r-rCut,2) + 2.0*B*(r-rCut); }
    else{result = (24.0*epsilon/r)*(pow(sigma/r,6)-2.0*pow(sigma/r,12));}
    return result;
}



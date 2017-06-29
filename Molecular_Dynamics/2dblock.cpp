// 2dblock.cpp
// HW4 Molecular dynamics 
// 2-D Block
// Author: Yuding Ai
// Date: 2017.03.21

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
        array<double,9> atom;
        atom[0] = xpos;
        atom[1] = ypos;
        atom[2] = atom[3] =atom[4] = atom[5] = atom[6]=0.0;
        atom[7] = 0; //set initial velocity to be 0;
        atom[8] = 0; //set initial velocity to be 0;
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

array<double,9> &block::operator[](int index){return atomlist[index];}


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
        st<< Rlist[i]<<" "<<glist[i]<<"\n"; 
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
        // st << atomlist[i][3]<<" "<<atomlist[i][4]<<endl;
        st << atomlist[i][3]<<" "<<atomlist[i][4]<<" "<<atomlist[i][5]<<" ";
        st << atomlist[i][6]<<" "<<atomlist[i][7]<<" "<<atomlist[i][8]<<endl;
    } 
    
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();

}
double block::total_penergy(){
    // this total_penergy is the total potential energy
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
                //calc total potential energy
                double E = total_penergy();
                //calc maximum force
                double F = qcalc_mforce();
                //record the state/configuration data
                st<<j<<" "<<E<<" "<<sp[0]<<" "<<sp[1]<<" ";
                st<<sp[2]<<" "<<sp[3]<<" "<<F<<endl;
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
            double E = total_penergy();
            //calc maximum force
            double F = qcalc_mforce();
            //record the state/configuration data
            st<<j<<" "<<E<<" "<<sp[0]<<" "<<sp[1]<<" ";
            st<<sp[2]<<" "<<sp[3]<<" "<<F<<endl;
        }

        j++;
        k++;
    }
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();
}

void block::set_config(string filename){
    string line;
    ifstream myfile(filename);
    if(myfile.is_open()){
        int i = 0;
        while (std::getline(myfile,line)){
            stringstream linestream(line);
            string data;
            std::getline(linestream,data,' ');
            atomlist[i][0] = std::stod(data);
            linestream >>atomlist[i][1]>>atomlist[i][2]>>
                atomlist[i][3]>>atomlist[i][4]
                >>atomlist[i][5]>>atomlist[i][6]
                >>atomlist[i][7]>>atomlist[i][8];
            i++;
        }
    }
}
void block::update_position(){
    /* According to Verlet algorithm         
     * r_i((j+1)*deltat) = r_i(j*deltat) + v_i(j*deltat)*deltat + 
     * (deltat)^2/2m_i *F_i(j*deltat)
     */

    //Recall in my atom's data structure, 
    //atom[0]  --- xposition
    //atom[1]  --- yposition
    //atom[5]  --- F_x
    //atom[6]  --- F_y
    //atom[7]  --- v_x
    //atom[8]  --- v_y
    //And when we calc r, F and v are all corresponding to j*delte_t
    for(unsigned int idx = 0; idx< atomlist.size();idx++){
        //x_position
        atomlist[idx][0] = atomlist[idx][0] + atomlist[idx][7]*delta_t 
            + delta_t*delta_t/(2*mass)* atomlist[idx][5];

        //y_position
        atomlist[idx][1] = atomlist[idx][1] + atomlist[idx][8]*delta_t 
            + delta_t*delta_t/(2*mass)* atomlist[idx][6];
    }
}

void block::update_velosity(array<double,400> F_oldx, array<double,400> F_oldy){
    /* According to Verlet algorithm         
     * v_i((j+1)*deltat) = v_i(j*deltat) + deltat/2m_i*(F_i((j+1)*deltat)
     * + F_i(j*deltat) )
     * 
     * To have F_i((j+1)*delta_t), we should first store the old Force 
     * F_i(j*delta_t) then update the force corresponding to (j+1)*delta_t
     * configuration
     */
    for(unsigned int idx = 0; idx< atomlist.size();idx++){
        //v_x
        atomlist[idx][7] = atomlist[idx][7] + delta_t/(2*mass) *(atomlist[idx][5]+ F_oldx[idx]);

        //v_y
        atomlist[idx][8] = atomlist[idx][8] + delta_t/(2*mass) *(atomlist[idx][6]+ F_oldy[idx]);
    }
}

void block::MD_mini(int it){
    stringstream st;
    stringstream st_temp_stress;
    array<double,4> totstress; //store the total stress and pressure
    array<double,10> pressurearray; //store the pressure of the 10 MD simulations
    array<double,10> temperaturearray; //store the temperature of the 10 MD simulations
    array<double,400> F_oldx;
    array<double,400> F_oldy;

    int counter = 0;
    //I have tried many runs on this model and
    //according to my experiments: 
    //when delta_t = 1E-14:
    //A. At least 30000 iterations is needed for the first MD
    //simulation to reach the equalibirum. For the next few MD simulations
    //, fewer and fewer steps is needed to reach its equalibirum. so 
    //here we choose 50000 iterations for each MD simulation.
    //and after 30000 iterations, we start to collect date as to calculate 
    //the equalibirum average.
    //
    //B. I found that 10 times of sucessive MD simulation with 50000 
    //iterations for each MD simulation would be sufficient enough 
    //to reach the minimun potential energy.(PE_min = -12.3343 eV) which is
    //even better than the result from my previous Molecular static method
    //(-11.7813 eV)
    //
    //Therefore, here we do 10 successive MD simulations to minimize 
    //the potential energy.
    
    for (int n = 0; n<10;n++){
        //For each MD, before looping, we first set vx and vy to be k_fac*v;
        for(unsigned int j = 0; j< atomlist.size();j++){
            atomlist[j][7] = 0;
            atomlist[j][8] = 0;
        }
        double KE= 0;
        double KEsum = 0;
        double Psum=0;

        //then get the initial Force and potential energy
        double PE = total_penergy();
        calc_force();

        int i = 0;
        while(i < it){
            for(unsigned int j = 0; j< atomlist.size();j++){
                F_oldx[j] = atomlist[j][5];
                F_oldy[j] = atomlist[j][6];
            }
            //for eacn iteration, first get the new position
            update_position();
            //then once atoms are moved, we update the neighbor_list
            get_neighbor_list(neighCut);
            // calc the new force
            calc_force();
            //then get the new velocity(the new force is been updated during
            update_velosity(F_oldx, F_oldy);

            //calc the PE and KE corresponding to the new configuration
            PE = total_penergy();
            KE = 0; //recalc the KE
            for (unsigned int k = 0; k<atomlist.size();k++){
                KE += 0.5* mass*(atomlist[k][7]*atomlist[k][7] + 
                        atomlist[k][8]*atomlist[k][8]);
            }
            //store the infomation into stringtream
            st<< KE <<" "<<PE<<" "<<counter<<endl;
            cout<< KE <<" "<<PE<<" "<<counter<<endl;
            i++;
            counter++;

            //According to my result, for delta_t = 1E-14, equalibirum would reached
            //after 30000 iterations, so from there, we start to record the KE value 
            //as to find the average <KE>, then to obtain temperature by N*k_b*T = <KE>
            //FYI, My program takes: 4213.73 sec to compele 10 sets of MD simulation
            //with 50000 itetarions each. so the computation is pretty slow
            
            if(i>30000){
                KEsum +=KE;
                totstress = calc_total_stress_pressure();
                Psum += totstress[3];
            }
        }
        
        KEsum /=(it - 30000); //get the <KE>
        Psum /=(it - 30000); //get the <P>
        double temp = KEsum/(400*k_b);
        get_atomic_stress(); 
        totstress = calc_total_stress_pressure();
        pressurearray[n] = totstress[3];
        temperaturearray[n] = temp;

        st_temp_stress<<temp <<" "<<totstress[0]<<" "<<totstress[1]<<" "<<totstress[2]
            <<" "<<totstress[3]<<" "<<Psum<<endl;
    }

    ofstream myfile1("MD_temp_stress.txt");
    string data1 = st_temp_stress.str();
    myfile1<< data1;
    myfile1.close();

    ofstream myfile("MD_Energy.txt");
    string data = st.str();
    myfile<< data;
    myfile.close();
}


void block::MD(int it,double desire_T,double old_T,string dir){
    // This program is basically the same as MD_mini using the same algorithm
    // but this time is to simulate the temperature dependence MD
    stringstream st;
    stringstream st_temp_stress;
    array<double,4> totstress; //store the total stress and pressure
    array<double,400> F_oldx;
    array<double,400> F_oldy;
    
    int counter = 0;
   
    //before starts the MD simulation, we first calculate k_factor based on the 
    //choice of desire_T and a preknown old_T which is obtained from previous simulation;
    //the very initial state, set old_T to be 1 and desire_T =  0
    //
    //as K = (2Nk_b*T_desire/sum(m_i*<vi^2>))^0.5 = (T_desire/T_old)^0.5
    
    double k_fac = pow(desire_T/old_T,0.5); //set initial k_fac
    // for (int n = 0; n<20;n++){
    for (int n = 0; n<30;n++){
        //For each MD, before looping, we first scale the veloclty
        for(unsigned int j = 0; j< atomlist.size();j++){
            atomlist[j][7] = atomlist[j][7]*k_fac;
            atomlist[j][8] = atomlist[j][8]*k_fac;
        }
        double KE= 0;
        double velocity = 0;
        double KEsum_tempdenpendent = 0;
        double hydroP=0;
        double Tdependent_temp = 0;

        //then get the initial Force and potential energy
        double PE = total_penergy();
        calc_force();

        int i = 0;
        while(i < it){
            //get the old force
            for(unsigned int j = 0; j< atomlist.size();j++){
                F_oldx[j] = atomlist[j][5];
                F_oldy[j] = atomlist[j][6];
            }           
            //for eacn iteration, first get the new position
            update_position();
            //then once atoms are moved, we update the neighbor_list
            get_neighbor_list(neighCut);
            // calc the new force
            calc_force();
            //then get the new velocity(the new force is been updated during
            update_velosity(F_oldx, F_oldy);

            //calc the PE and KE corresponding to the new configuration
            PE = total_penergy();
            KE = 0; //recalc the KE
            for (unsigned int k = 0; k<atomlist.size();k++){
                KE += 0.5* mass*(atomlist[k][7]*atomlist[k][7] + 
                        atomlist[k][8]*atomlist[k][8]);
                velocity += pow(atomlist[k][7]*atomlist[k][7]+
                        atomlist[k][8]*atomlist[k][8],0.5);

            }

            get_atomic_stress(); 
            totstress = calc_total_stress_pressure();
            velocity /=atomlist.size(); 
            cout<< KE <<" "<<PE<<" "<<" "<<velocity<<" "<<counter<<" "<<totstress[3]<<endl;
            st<< KE <<" "<<PE<<" "<<" "<<velocity<<" "<<counter<<" "<<totstress[3]<<endl;
            i++;
            counter++;
            KEsum_tempdenpendent +=KE;
            hydroP +=totstress[3];

        }
        
        //record the final config of each MD simulation
        KEsum_tempdenpendent /= (it); //get <KE>
        hydroP /= (it); //get <P>
        // autocorr /=(it-1000); //get <v(t)v(0)>
        Tdependent_temp = KEsum_tempdenpendent/(400*k_b);
        get_atomic_stress(); 
        totstress = calc_total_stress_pressure();

        cout<<k_fac<<" "<<Tdependent_temp<<" "<<i<<endl;
        st_temp_stress<<Tdependent_temp <<" "<<totstress[0]<<" "<<totstress[1]
            <<" "<<totstress[2]<<" "<<totstress[3]<<" "<<hydroP<<endl;

        //update the k_fac
        k_fac = pow(desire_T/Tdependent_temp,0.5);
    }

    string filename1 = dir +"/part2_MD_temp_stress.txt";
    ofstream myfile1(filename1);
    string data1 = st_temp_stress.str();
    myfile1<< data1;
    myfile1.close();

    string filename = dir +"/part2_MD_Energy.txt";
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();
}

void block::autocorrelation(int it,string dir){
    stringstream st;
    stringstream st_temp_stress;
    array<double,400> F_oldx;
    array<double,400> F_oldy;
    vector<array<double,400> > v_all;

    int counter = 0;

    double KE= 0;
    double KEsum_tempdenpendent = 0;

    //then get the initial Force and potential energy
    double PE = total_penergy();
    calc_force();

    int i = 0;
    counter = 0;

    double velocity = 0;
    while(i < it){
        velocity = 0;
        //get the old force
        for(unsigned int j = 0; j< atomlist.size();j++){
            F_oldx[j] = atomlist[j][5];
            F_oldy[j] = atomlist[j][6];
        }
        //for eacn iteration, first get the new position
        update_position();
        //then once atoms are moved, we update the neighbor_list
        get_neighbor_list(neighCut);
        // calc the new force
        calc_force();
        //then get the new velocity(the new force is been updated during
        update_velosity(F_oldx, F_oldy);

        //calc the PE and KE corresponding to the new configuration
        PE = total_penergy();
        KE = 0; //recalc the KE
        array<double,400> v_run;
        for (unsigned int k = 0; k<atomlist.size();k++){
            KE += 0.5* mass*(atomlist[k][7]*atomlist[k][7] +
                    atomlist[k][8]*atomlist[k][8]);

            velocity = pow(atomlist[k][7]*atomlist[k][7]+
                    atomlist[k][8]*atomlist[k][8],0.5);
            v_run[k] = velocity;
        }
        v_all.push_back(v_run);
        cout << counter<<"velocity "<<KE <<" "<<velocity<<endl;
        
        i++;
        counter++;
        KEsum_tempdenpendent +=KE;
    }

    double autoco = 0;
    double vo = 0;
    double vt = 0;
    for (int k = 0; k<5000;k++){
        autoco = 0;
        for(int i = 0; i<5000;i++ ){
            for(int j = 0; j<400;j++){
                vo = v_all[i][j];
                vt = v_all[i+k][j];
                autoco = autoco + vo*vt;
            }
        }
        autoco = autoco/(400*5000);
        st<<autoco<<" "<<k<<endl;
    }

    string filename =dir + "/auto.txt";
    ofstream myfile(filename);
    string data = st.str();
    myfile<< data;
    myfile.close();
}

double block::diffusion(int it){
    stringstream st;
    stringstream st_temp_stress;
    array<double,400> F_oldx;
    array<double,400> F_oldy;
    vector<array<double,400> > r_x;
    vector<array<double,400> > r_y;

    int counter = 0;

    double KE= 0;
    double KEsum_tempdenpendent = 0;

    //then get the initial Force and potential energy
    double PE = total_penergy();
    calc_force();

    int i = 0;
    counter = 0;

    double velocity = 0;
    while(i < it){
        velocity = 0;
        //get the old force
        for(unsigned int j = 0; j< atomlist.size();j++){
            F_oldx[j] = atomlist[j][5];
            F_oldy[j] = atomlist[j][6];
        }
        //for eacn iteration, first get the new position
        update_position();
        //then once atoms are moved, we update the neighbor_list
        get_neighbor_list(neighCut);
        // calc the new force
        calc_force();
        //then get the new velocity(the new force is been updated during
        update_velosity(F_oldx, F_oldy);

        //calc the PE and KE corresponding to the new configuration
        PE = total_penergy();
        KE = 0; //recalc the KE
        array<double,400> r_runx;
        array<double,400> r_runy;
        for (unsigned int k = 0; k<atomlist.size();k++){
            KE += 0.5* mass*(atomlist[k][7]*atomlist[k][7] +
                    atomlist[k][8]*atomlist[k][8]);

            r_runx[k] = atomlist[k][0];
            r_runy[k] = atomlist[k][1];
        }
        r_x.push_back(r_runx);
        r_y.push_back(r_runy);
        cout << counter<<"velocity "<<KE <<" "<<velocity<<endl;
        
        i++;
        counter++;
        KEsum_tempdenpendent +=KE;
    }

    double rox = 0;
    double roy = 0;
    double rtx = 0;
    double rty = 0;
    double D= 0;
    for (int k = 0; k<5000;k++){
        D = 0;
        for(int i = 0; i<5000;i++ ){
            for(int j = 0; j<400;j++){
                rox = r_x[i][j];
                roy = r_y[i][j];
                rtx = r_x[i+k][j];
                rty = r_y[i+k][j];
                D += (rtx-rox)*(rtx-rox)+(rty-roy)*(rty-roy);
            }
        }
        D = D/(400*5000);
    }
    D /= (6*5000);
    double lnD=0;
    lnD = log(D);
    return lnD;
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



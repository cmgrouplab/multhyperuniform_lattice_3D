//generate hopefully multihyperuniform packings on a given 2D lattice 
// -> this is to understand the competition between local geometrical frustrations and possible long-range interaction
// -> the hard-core exclusion is gauranteed by the singular occuptation nature of the lattice
// -> the long range same-species interaction is imposed by maximizing the sum of nearest neighbor distance


//author: Yang Jiao, yang.jiao.2@asu.edu
//started: 01/20/2024
//adapted 10/25/24 for HEA Carbide with NaCl structure

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define MAXS 20000 //the upper bound for random number generator
#define SD 3 //space dimension

double pi = 3.141592654;

#define MAXX 10 //the linear dimension of the lattice, applying to all lattice, w.r.t. the unit cell, HAS TO BE EVEN
#define N_tot MAXX*MAXX*MAXX/2 //the total number of metal particles, the remaining half are carbons

#define n_s 2 //number of species
int* N_c; //number of "cells" of each specides
int* Type_c;//the type of cell with an index


int config[MAXX][MAXX][MAXX]; //this gives the index of particle occupying each lattice site, -1 for unoccupied sites or vacancies

int** particle_group; //collection of index for all particle of a particluar type, for more efficient search and update
double** coords; //Euclidean coordinates of the cells on lattice, as a function of index
double** dist_matrix; //distance matrix between any particle i and j, regardless of types, 
                     //efficient search can be achieved by only work with index for a particular group
double* dist_nn; // the nearest neighbor same-species distance list, this will be updated if particles are moved;
                 //this will also be used for getting the energy
                 //it is very important to 
int* dist_nn_index; //this saves the index of particle that give the neareast distance                  
                 

double box_length = 1.0; //linear size of the simulation box
                 
double ell = 1.0/(double)MAXX; //this is the length scale to convert the total box length unity

double* Range; //the interaction range of like-species
                 
//***********************
//the following are for annealing and updates

int particle_index_I; //the index of the two switched particles
int particle_index_J;

int indexi; //this is the matrix indices for the sites
int indexj;
int indexk;
int indexm;
int indexn;
int indexq;

double coords_I_old[SD]; //the original coords of two switched particles
double coords_J_old[SD];

//int matrix_index_I[SD]; //the original matrix entries for the two switched pixels 
//int matrix_index_J[SD];

  

//************************
//the anealing paramters

double alpha = 0.9; //the cooling schedule parameter
int TN = 40; //the cooling stages
double T = 0.0015; //the current temperature

int Nevl = 5000; //the number of evolution per stage, this depends on the total number of particles

int N_stage = 2000; //this is for brute-froce steepest decent stages, not used for annealing

//*****************************


void read_parameters()
{
	cout<<"reading parameters now ...."<<endl;
	
	//first read in particle parameters ...
	 cout<<"total number of different species n_s = "<<n_s<<endl;
 	 N_c = new int[n_s];
 
 	int temp_N_tot = 0;

 	 for(int i=0; i<n_s; i++)
    	{
      	cout<<"number of cells for species "<<i<<" N_c["<<i<<"] = "; cin>>N_c[i];
      	temp_N_tot += N_c[i];
    	}
    	
    if(temp_N_tot != N_tot)
	{
		cout<<"total particle number if not right, adding up numbers for all species does not give N_tot!"<<endl;
		exit(1);
	}	

  	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  
  
    Type_c = new int[N_tot];
  
   coords = new double*[N_tot];
  	for(int i=0; i<N_tot; i++)
    {
      coords[i] = new double[SD];
    }

	particle_group = new int*[n_s];  //save the index of particle of species i
	for(int i=0; i<n_s; i++)
	{
		particle_group[i] = new int[N_c[i]];
	}
	
	dist_matrix = new double*[N_tot]; //saves the distance between any pair of particles in the system
	for(int i=0; i<N_tot; i++)
   {
   	   dist_matrix[i] = new double[N_tot]; 
   }


  dist_nn = new double[N_tot]; //saves the same-species nearest neighbor distance, will be used to compute energy
  
  dist_nn_index = new int[N_tot]; //save the index that gives dist_nn for partilce i
  
  Range = new double[n_s];
  
   for(int i=0; i<n_s; i++)
    {
      Range[i] = sqrt(box_length*box_length/(N_c[i]));
      //this is the RSA density 
      cout<<"MAXIMAL interaction range of cells for species "<<i<<" Range["<<i<<"] = "<<Range[i]<<endl;
      cout<<"Rescale factor (preferred to be between 0 and 2, larger values might cause error) for species" <<i<<endl;
      
      //now allow rescaling the soft shell
      double tmp_range_scale;
      cin>>tmp_range_scale;
      
      Range[i] = tmp_range_scale*Range[i];
      cout<<"RESCALED interfaction range of cells for species "<<i<<" Range["<<i<<"] = "<<Range[i]<<endl;
      
    }
 
  cout<<"**************************"<<endl;
  cout<<"the default annealing parameters are below:"<<endl;
  
   //now for the cooling schedule
  cout<<"Starint temp T0 = "<<T<<endl;
  cout<<"Decreasing ratio: alpha = "<<alpha; cout<<endl;
  cout<<"Number of decreasing T stages: TN = "<<TN; cout<<endl;
  cout<<"Number of pxiel move per stage: Nevl = "<<Nevl; cout<<endl;
  
  int temp_flag = 0;
  cout<<"type 1 to reset the parameters, and type 0 to keep the current value..."<<endl;
  cin>>temp_flag;	

 	if(temp_flag == 1)
 	{
	 
  //now for the cooling schedule
  cout<<"Starint temp T0 = "; cin>>T; cout<<endl;
  cout<<"Decreasing ratio: alpha = "; cin>>alpha; cout<<endl;
  cout<<"Number of decreasing T stages: TN = "; cin>>TN; cout<<endl;
  cout<<"Number of pxiel move per stage: Nevl = "; cin>>Nevl; cout<<endl;
 }
  else
  {	
  	cout<<"default annealing parameters are used!"<<endl;
  }
  
}

void init_config()
{
	cout<<"initializing configurations ..."<<endl;
	
	//set all entry to be -1 indicating non-occuptation 
	for(int i=0; i<MAXX; i++)
		for(int j=0; j<MAXX; j++)
			for(int k=0; k<MAXX; k++)
			config[i][j][k] = -1;

	//now put in the Carbons, taking the FCC subsite
    for(int i=0; i<MAXX/2; i++)
		for(int j=0; j<MAXX/2; j++)
			for(int k=0; k<MAXX/2; k++)
			{
				config[2*i][2*j][2*k] = -10;
				config[2*i+1][2*j+1][2*k] = -10;
				config[2*i+1][2*j][2*k+1] = -10;
				config[2*i][2*j+1][2*k+1] = -10;
			}
			
	
	
	int temp_ct = 0; //this is a global counter for all particles 
	
	int temp_sub_ct = 0; //this is for each species
	
	//now put in the first n_s-1 speices by randomly selecting avaialble sites..
	for(int i=0; i<(n_s-1); i++)
	{
		cout<<"initializing species "<<i<<endl;
		
		temp_sub_ct = 0; //this is for each species
		
		for(int j=0; j<N_c[i]; j++)
		{
			int m = rand()%MAXX;
     		int n = rand()%MAXX;
			int q = rand()%MAXX;
    

    	    while(config[m][n][q]!=-1)
			{
	  			m = rand()%MAXX;
	  			n = rand()%MAXX;
				q = rand()%MAXX;
			}

      		config[m][n][q] = temp_ct; //put the index of the particel at the site
      		
      		//map the matrix index to coords
      		coords[temp_ct][0] = ell*m;
      		coords[temp_ct][1] = ell*n;
			coords[temp_ct][2] = ell*q;
      		
      		//save the type of this index
			Type_c[temp_ct] = i; //i is the species index 
      		
      		//add this type to the corresponding group 
      		particle_group[i][temp_sub_ct] = temp_ct;
      		
      		
      		//update the counter
      		temp_ct++;
      		
      		temp_sub_ct++;
      		
      		if(temp_sub_ct > N_c[i])
      		{
      			cout<<"temp_sub_ct = "<<temp_sub_ct<<" > N_c["<<i<<"] = "<<N_c[i]<<endl;
      			exit(1);
			  }
		}
	}
	
	
	//*******************************************************
	cout<<"initializing species "<<(n_s -1 )<<endl;
	
	temp_sub_ct = 0;
	
	//last put in the last species by taking up all available sites
	for(int i=0; i<MAXX; i++)
		for(int j=0; j<MAXX; j++)
		for(int k=0; k<MAXX; k++)
		{
			if(config[i][j][k] == -1)
			{
				config[i][j][k] = temp_ct; //put the index of the particel at the site
				
				//map the matrix index to coords
      			coords[temp_ct][0] = ell*i;
      			coords[temp_ct][1] = ell*j;
				coords[temp_ct][2] = ell*k;
      		
      			//save the type of this index
				Type_c[temp_ct] = n_s-1; // the species index 
      		
      			//add this type to the corresponding group 
      			particle_group[n_s-1][temp_sub_ct] = temp_ct;
				
				//update the counter
      			temp_ct++;
      		
      			temp_sub_ct++;
      		
      			//cout<<"i = "<<i<<endl;
      		
      			if(temp_sub_ct > N_c[n_s-1])
      			{
      				cout<<"temp_sub_ct = "<<temp_sub_ct<<" > N_c["<<(n_s -1)<<"] = "<<N_c[i]<<endl;
      				exit(1);
			  	}
			}
		}
		
	cout<<"finishing initalizing configurations"<<endl;
	cout<<"temp_ct = "<<temp_ct<<endl;
	
	//testing the matrix
	/*
	cout<<"****************"<<endl;
	for(int i=0; i<MAXX; i++)
	{
		for(int j=0; j<MAXX; j++)
			cout<<config[i][j]<<"\t";
		cout<<endl;
	}
	*/
		
	cout<<"****************"<<endl;
	//testing Type_c
	for(int i=0; i<N_tot; i++)
		cout<<"type "<<i<<" = "<<Type_c[i]<<endl;
	
	
	
}

double get_dist(int index1, int index2)
{
  double tempx = fabs(coords[index1][0] - coords[index2][0]);
  if(tempx >= box_length/2.0) tempx = box_length - tempx;

  double tempy = fabs(coords[index1][1] - coords[index2][1]);
  if(tempy >= box_length/2.0) tempy = box_length - tempy;

  double tempz = fabs(coords[index1][2] - coords[index2][2]);
  if(tempz >= box_length/2.0) tempz = box_length - tempz;

  return sqrt(tempx*tempx+tempy*tempy+tempz*tempz);

  
}


double get_dist(double center1[SD], double center2[SD])
{
  double tempx = fabs(center1[0] - center2[0]);
  if(tempx >= box_length/2.0) tempx = box_length - tempx;

  double tempy = fabs(center1[1] - center2[1]);
  if(tempy >= box_length/2.0) tempy = box_length - tempy;

  double tempz = fabs(center1[2] - center2[2]);
  if(tempz >= box_length/2.0) tempz = box_length - tempz;

  return sqrt(tempx*tempx+tempy*tempy+tempz*tempz);

  
}


void init_data()
{
	cout<<"initializing data ..."<<endl;
	
	double temp_dist;
	
	//compute the dist_matrix for all, the full calculation will only be carried out once
	for(int i=0; i<N_tot; i++)
		for(int j=0; j<N_tot; j++)
		{
			temp_dist = get_dist(i, j);
			
			dist_matrix[i][j] = temp_dist;
			}	
	
	//use the dist_matrix to get dist_nn list
	
	double mini_dist;
	int temp_type;
	int temp_index;
	int best_index;
	
	for(int i=0; i<N_tot; i++)
	{
		mini_dist = 100000000.0;
		
		temp_type = Type_c[i]; 
		
		//loop over only the sub-group associated with this species
		for(int j=0; j<N_c[temp_type]; j++)
		{
			temp_index = particle_group[temp_type][j]; //take out the index of the sub gropu list 
			
			//update the mini_dist
			if(i != temp_index && mini_dist>dist_matrix[i][temp_index])
			{
				mini_dist = dist_matrix[i][temp_index];
				
				best_index = temp_index;
			}
		}
		
		//once mini_dist is found by checking all same-type particle, put this in the list
		dist_nn[i] = mini_dist;
		dist_nn_index[i] = best_index;
		
		if(Type_c[i] != Type_c[best_index])
		{
			cout<<"for mini_dist calculation, the particle types are not matching! Recheck!"<<endl;
			exit(1);
		}
		
		cout<<"dist_nn ["<<i<<"] = "<<dist_nn[i]<<endl;
		cout<<"dist_nn_index ["<<i<<"] = "<<dist_nn_index[i]<<endl;
	}
}


void change_config()
{
  int i, j, k, m, n, q;
  int lim = 0;
  
 
  do{   i = rand() % MAXX;
        m = rand() % MAXX;

        j = rand() % MAXX; 
        n = rand() % MAXX;

		k = rand() % MAXX; 
        q = rand() % MAXX;

        lim++; //this is not used

     } while(config[i][j][k]<0 || config[m][n][q]<0 || Type_c[config[i][j][k]] == Type_c[config[m][n][q]]);


  //switch the particle index in these sites
  particle_index_I = config[i][j][k];
  particle_index_J = config[m][n][q];
  
  //cout<<"index_I = "<<particle_index_I<<endl;
  //cout<<"index_J = "<<particle_index_J<<endl;
  
  int temp;
  temp = config[i][j][k];
  config[i][j][k] = config[m][n][q];
  config[m][n][q] = temp;

  indexi = i;
  indexj = j;
  indexk = k;
  
  indexm = m;
  indexn = n;
  indexq = q;
  
  
   //now save the old coords...
  coords_I_old[0] = ell*indexi;
  coords_I_old[1] = ell*indexj;
  coords_I_old[2] = ell*indexk;
  
  coords_J_old[0] = ell*indexm;
  coords_J_old[1] = ell*indexn;
  coords_J_old[2] = ell*indexq;
  
  
  //now updates the coords...
  coords[particle_index_I][0] = ell*indexm;
  coords[particle_index_I][1] = ell*indexn;
  coords[particle_index_I][2] = ell*indexq;
  
  coords[particle_index_J][0] = ell*indexi;
  coords[particle_index_J][1] = ell*indexj;
  coords[particle_index_J][2] = ell*indexk;
  
  //update the distance matrix 
  for(int t=0; t<N_tot; t++)
  {
  	dist_matrix[particle_index_I][t] = dist_matrix[t][particle_index_I] = get_dist(particle_index_I, t);
  	
  	dist_matrix[particle_index_J][t] = dist_matrix[t][particle_index_J] = get_dist(particle_index_J, t);
  }
  
  
  //update the same-species nnl, for new energy calcuulation...
 

}

void resume_config()
{
  int temp;
  //first we resume the config
  temp = config[indexi][indexj][indexk];
  config[indexi][indexj][indexk] = config[indexm][indexn][indexq];
  config[indexm][indexn][indexq]= temp;
  
  
   //now resume the coords...
  coords[particle_index_I][0] = ell*indexi;
  coords[particle_index_I][1] = ell*indexj;
  coords[particle_index_I][2] = ell*indexk;
  
  coords[particle_index_J][0] = ell*indexm;
  coords[particle_index_J][1] = ell*indexn;
  coords[particle_index_J][2] = ell*indexq;
  
  //resume the distance matrix, by re-computing  
  for(int t=0; t<N_tot; t++)
  {
  	dist_matrix[particle_index_I][t] = dist_matrix[t][particle_index_I] = get_dist(particle_index_I, t);
  	
  	dist_matrix[particle_index_J][t] = dist_matrix[t][particle_index_J] = get_dist(particle_index_J, t);
  }
  
  
  //resume the same-species nnl, 

}


//loop over all particles of the same type, directly recompute the distance, use soft shell interaction
double get_energy(int index, double center[SD])
{
	double ener = 0;
	
	int temp_type = Type_c[index];
	int temp_index;
	
	for(int i=0; i<N_c[temp_type]; i++)
    {
    	temp_index = particle_group[temp_type][i];
    	
    	//cout<<"temp_index"<<temp_index<<endl;
    	
      if(temp_index != index)
	{
	  double temp_dist = get_dist(center, coords[temp_index]) - Range[Type_c[index]];

	  if(temp_dist<0)
	    {
	      ener += temp_dist*temp_dist;
	    }
	}
    }

  return ener;
}


//loop over all particles of the same type, use data in distance marix, use soft shell interaction
double get_energyII(int index, double center[SD])
{
	double ener = 0;
	
	int temp_type = Type_c[index];
	int temp_index;
	
	for(int i=0; i<N_c[temp_type]; i++)
    {
    	temp_index = particle_group[temp_type][i];
    	
    	//cout<<"temp_index"<<temp_index<<endl;
    	
      if(temp_index != index)
	{
	  double temp_dist = dist_matrix[index][temp_index] - Range[Type_c[index]];

	  if(temp_dist<0)
	    {
	      ener += temp_dist*temp_dist;
	    }
	}
    }

  return ener;
}


//directly recompute the distance, use soft shell interaction
double get_total_ener()
{
	double en_tot = 0;
	for(int i=0; i<N_tot; i++)
    en_tot +=  get_energy(i, coords[i]);
    
    return en_tot;
}

//use data in distance marix, use soft shell interaction
double get_total_enerII()
{
	double en_tot = 0;
	for(int i=0; i<N_tot; i++)
    en_tot +=  get_energyII(i, coords[i]);
    
    return en_tot;
}


//another way to define energy is to sum over all distance, WITHIN Range < box_length/2, for a given speces
//this way, one can control the range of interaction, this is similar to the idea of soft shell,
//but more flexible for incorporating long-range interactions.
//since the lattice is highly degenerate in terms of near-neighbor distances, this is hard to tune to get correct results
double get_total_enerIII()
{
	double en_tot = 0;
	
	//loop over the species
	for(int i=0; i<n_s; i++)
	{
		//loop the particles in the sub-group
		int temp_index1, temp_index2;
		double temp_dist;
		
		for(int m=0; m<N_c[i]; m++)
			for(int n=0; n<N_c[i]; n++)
			{
				//get two indices of the same species
				temp_index1 = particle_group[i][m];
				temp_index2 = particle_group[i][n];
				
				temp_dist = dist_matrix[temp_index1][temp_index2];
				
				if(temp_index1 != temp_index2 && temp_dist < Range[i])
					en_tot += temp_dist;
			}
	}
    
    
    return -en_tot;
}


/*
//this is not valid any more, because two particles of different kind moved simutaneously, and there will be non-trivial cross-terms
double get_cost(int index, double center[SD], double center_old[SD])
{
  double cost = get_energy(index, center) - get_energy(index, center_old);

  return cost;
}

*/


//will modify this to allow printing out intermediate configurations
void print_config()
{
  ofstream fout;
  fout.open("centers.xls");
  
  cout<<"printing out configurations for all species..."<<endl;
  
  int counter = 0;
  int sum = 0;

  fout<<N_tot<<endl<<endl;

  for(int s=0; s<n_s; s++)
    {
      counter = 0;

      fout<<N_c[s]<<endl;
      //fout<<Rad[s]<<endl;
      fout<<Range[s]<<endl;
      
      for(int i=sum; i<N_tot; i++)
	{
	  if(Type_c[i] == s)
	    {
	      fout<<coords[i][0]<<"\t"<<coords[i][1]<<"\t"<<coords[i][2]<<endl;

	      counter++;
	    }

	  if(counter == N_c[s])
	    {
	      sum += counter;
	      break;
	    }
	}

      fout<<"//###########################"<<endl;
      fout<<endl<<endl;

    }

  fout.close();
}

void print_csv()
{
  ofstream fout;
  fout.open("centers.csv");
  
  cout<<"printing out the CSV files for all species..."<<endl;
  
  int counter = 0;
  int sum = 0;

  //fout<<N_tot<<endl<<endl;

  fout<<"Element, x, y, z"<<endl;


//print out the metals first
  for(int s=0; s<n_s; s++)
    {
      counter = 0;

      //fout<<N_c[s]<<endl;
      //fout<<Rad[s]<<endl;
      //fout<<Range[s]<<endl;
      
      for(int i=sum; i<N_tot; i++)
	{
	  if(Type_c[i] == s)
	    {
		  if(s==0)	
	         fout<<"Mo,"<<coords[i][0]<<","<<coords[i][1]<<","<<coords[i][2]<<endl;
		  else if(s==1)	
	         fout<<"Nb,"<<coords[i][0]<<","<<coords[i][1]<<","<<coords[i][2]<<endl;
		  else if(s==2)	
	         fout<<"Sc,"<<coords[i][0]<<","<<coords[i][1]<<","<<coords[i][2]<<endl;
		  else if(s==3)	
	         fout<<"Ta,"<<coords[i][0]<<","<<coords[i][1]<<","<<coords[i][2]<<endl;
		  else if(s==4)	
	         fout<<"Ti,"<<coords[i][0]<<","<<coords[i][1]<<","<<coords[i][2]<<endl;

	      counter++;
	    }

	  if(counter == N_c[s])
	    {
	      sum += counter;
	      break;
	    }
	}

      //fout<<"//###########################"<<endl;
      //fout<<endl<<endl;

    }
	
  //now print out the carbons, by looping overall elements
  for(int i=0; i<MAXX; i++)
	for(int j=0; j<MAXX; j++)
		for(int k=0; k<MAXX; k++)
		{
			if(config[i][j][k]==-10)
				fout<<"C,"<<ell*i<<","<<ell*j<<","<<ell*k<<endl;
		}

  fout.close();
}


double PE(double dE, double T) // the probability that should compare ...
{
  if(dE > 0) return exp(-dE/T);
  else return 1;
}



//**********************
int main()
{
	srand(time(NULL));
	
	//cout<<"testing testing testing"<<endl;
	cout<<"N_tot = "<<N_tot<<endl;

	cout<<"ell = "<<ell<<endl;	
	cout<<"*****************************"<<endl;
	
	read_parameters();
	
	init_config();
	
	init_data();
	
	
	cout<<"****************"<<endl;
	cout<<"start energy minimization..."<<endl;
	
	/*
	for(int i=0; i<N_tot; i++)
		cout<<"E_old "<<i<<" = "<<get_energy(i, coords[i])<<"="<<get_energyII(i, coords[i])<<endl;
  	cout<<"E_old_tot = "<<get_total_ener()<<endl;
  	cout<<"E_old_tot = "<<get_total_enerII()<<endl;
  	cout<<"E_old_tot = "<<get_total_enerIII()<<endl;
	*/
	
	double E_old, E_new; //save the energy for the optimization 
	
	E_old = get_total_enerII();
	
	int N_acc = 0;
	
	
	//now start the annealing processs...
	for(int t= 0; t<TN; t++)
	{
		N_acc = 0;
		
		T = alpha*T;
		
		cout<<"********************"<<endl;
		cout<<"annealing stage "<<t<<" started ...."<<endl;
		cout<<"current temperature T = "<<T<<endl;
	
	for(int k = 0; k<Nevl; k++)
	{
	
	change_config();
	
	//cout<<"E_new_I = "<<get_energy(particle_index_I, coords[particle_index_I])<<endl;
	//cout<<"E_new_J = "<<get_energy(particle_index_J, coords[particle_index_J])<<endl;
	
	/*
	cout<<"****************"<<endl;
	for(int i=0; i<N_tot; i++)
		cout<<"E_new "<<i<<" = "<<get_energy(i, coords[i])<<"="<<get_energyII(i, coords[i])<<endl;
	
	cout<<"E_new_tot = "<<get_total_ener()<<endl;
	cout<<"E_old_tot = "<<get_total_enerII()<<endl;
  	cout<<"E_old_tot = "<<get_total_enerIII()<<endl;
	*/
	
	E_new = get_total_enerII();
	
	double temp_p = (double) (rand()%MAXS)/(double)MAXS;
	
	if(temp_p < PE((E_new-E_old), T))
	{
		E_old = E_new;
		
		N_acc ++;
	}
	else
	{
		resume_config();
	}
	
	
	//cout<<"E_old_I = "<<get_energy(particle_index_I, coords[particle_index_I])<<endl;
	//cout<<"E_old_J = "<<get_energy(particle_index_J, coords[particle_index_J])<<endl;
	
	//cout<<"****************"<<endl;
	/*
	for(int i=0; i<N_tot; i++)
		cout<<"E_old "<<i<<" = "<<get_energy(i, coords[i])<<"="<<get_energyII(i, coords[i])<<endl;
	
	cout<<"E_old_tot = "<<get_total_ener()<<endl;
	cout<<"E_old_tot = "<<get_total_enerII()<<endl;
  	cout<<"E_old_tot = "<<get_total_enerIII()<<endl;
  	*/
  	
  	//cout<<"stage "<<k<<"\t E = "<<E_old<<endl;
  	
  }
  
  cout<<"P_acc = "<<(double)N_acc/(double)Nevl<<endl;
  cout<<"E = "<<E_old<<endl;
  
 }
  
  print_config(); 

  print_csv(); 
	
	return 0;
}

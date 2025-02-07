//CPP file for obtaining statistics of the packing, including:
//     (1) pair correlation functions
//     (2) number variance
//     (3) structure factor  


//version: 04/10/2020
//author: duyu chen

//modified by Yang Jiao on 07.20.22 
//implement node and edge-based sampling to study effects of correlation on variance scaling

//modified 10/14/22
//read in a hyperuniform system and apply random but tiny pertrubtations ...

//modified 12/14/22
//read in 3D packing, computing the corresponding quantities, starting with Sk first ....

using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <random>

//#define N 2500 
#define N 5000
//#define N 421
#define SD 3


double Center[N][SD];
double L[SD];  //side length of box
//double Ly;
//double Lz;

#define Nbin 100 //for g(r)
#define Bin 0.005 //this needs to be consistent with edge length, e.g., L = 1, then Nbin*Bin <=0.5*L

double G[N][Nbin];//this is the local un-normalized g2, to compute node-based variance 
double BG[3*N][Nbin]; //this is for bond-based statistics 


//need to check the change these numbers for SK calculation
int Nk = 20; //this is the number of mesh point along each direction in K space, the unit of mesh is 2pi/L
double Kbin = 6.0; //this should be on the magnitude of 2pi/L 
#define Kbin_num 50

double nnl; //this is the diameter
double cut_off = 1.01; //this is to determine the connectivity network, 1% of diameter as tolerance, should be sufficient 

/*
int Nk = 300;
double Kbin = 0.15;
#define Kbin_num 100
*/

double pi = 3.14159265358979;
#define Rbin 0.002  //for number variance
#define Nr 1000 //sampling points for number variance at each r
#define MAXY 20000 //the maximal number is 32767 

ofstream fout;
ifstream fin;

void read_packing_3D()
{
	cout<<"Make sure to modify the packing data so it only includes N and D !"<<endl;
	cout<<"if this is confirmed, type 1 "<<endl;
	int temp_a; cin>>temp_a;
	
	
	fin.open("write.dat");
	fin>>temp_a; 
	if(temp_a != N)
	{
		cout<<"the particle number if not correct! check again!"<<endl;
		exit(1);
	}
	
	fin>>nnl; //this is the sphere diameter
	cout<<"sphere diamter D = "<<nnl<<endl;
	
	cout<<"reading the sphere center ..."<<endl;
	for(int i=0; i<N; i++)
		for(int j=0; j<SD; j++)
			fin>>Center[i][j];
			
	cout<<"initalizing the box length to be unit in all directions, if OK type 1"<<endl;
	cin>>temp_a;
	
	for(int i=0; i<SD; i++)
		L[i] = 1.0;
	
	fin.close();
		
	
}

/*
void read_config()
{
	FILE * fp;
//if ((fp = fopen("../../10000/0.04/test2/vertex_graphene_process.txt", "r")) == NULL)
	if ((fp = fopen("vertex_relax.txt", "r")) == NULL)
	{
		perror("Cannot open file!\n");
		exit(1);
	}
	else
	{
		double temp_t;
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &Lx);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &Ly);
		cout << "Lx: " << Lx << endl;
		cout << "Ly: " << Ly << endl;
		//L = 1029.0 / L;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < SD; j++)
			{
				fscanf(fp, "%lf", &Center[i][j]);
			}
		//	for(int i = 0; i < N; i++)
		//		std::cout << Center[i][0] << " " << Center[i][1] << std::endl;
		//Nbin = (int)floor(Ly / 2.0 / Bin);
		fclose(fp);
	}
}
*/



double MinDis(int m, int n)
{
  //find the minimal distance between the centers of two polyhedra in Ecludean space...
  //by checking all the images of Poly[m], while keeping Poly[n] in the central box
  //record the index of the box which Poly[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  double dx = fabs(Center[m][0] - Center[n][0]);
  double dy = fabs(Center[m][1] - Center[n][1]);
  double dz = fabs(Center[m][2] - Center[n][2]);
  //double dh = fabs(Center[m][3] - Center[n][3]);
  
  
  if(dx >= L[0]/2.0) dx = L[0] - dx;
  if(dy >= L[1]/2.0) dy = L[1] - dy;
  if(dz >= L[2]/2.0) dz = L[2] - dz;
  //if(dh >= L[3]/2.0) dh = L[3] - dh;
  
  	
  return sqrt(dx*dx + dy*dy + dz*dz); 
 

 /*
  double dist = 1000000000.0; //just a large number...


  //loop over all possible images of Point m, keep n fixed in the center simulation box....
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
	{
	  double tempd[SD]; 
	  tempd[0] = dx + (double)i*Lx;
	  tempd[1] = dy + (double)j*Ly;
	  
	  double tempdist = tempd[0]*tempd[0]+tempd[1]*tempd[1]; 
	  
	  //printf("TempDist[%d][%d][%d] = %f\n", i, j, k, tempdist);
	  
	  if(tempdist < dist) // store the smallest distance...
	    {
			dist = tempdist; 
	    }
	  
	}
  
  
  return sqrt(dist); //this is center-to-center distance...
  */
}

/*
double MinDis(double x1, double y1, double z1, double h1, double x2, double y2, double z2, double h2)
{
  //this is a simplified version, for a rectanglar box, no need to check all images, just the orthogonal minimal distance should be OK
  
  
  double dx = fabs(x1 - x2);
  double dy = fabs(y1 - y2);
  double dz = fabs(z1 - z2);
  double dh = fabs(h1 - h2);
  
  if(dx >= L[0]/2.0) dx = L[0] - dx;
  if(dy >= L[1]/2.0) dy = L[1] - dy;
  if(dz >= L[2]/2.0) dz = L[2] - dz;
  if(dh >= L[3]/2.0) dh = L[3] - dh;
  
  	
  return sqrt(dx*dx + dy*dy + dz*dz + dh*dh); 
  
 
 
}
*/


double Get_NDensity() // the number density, for computing g2...
{
  double VLambda;
  VLambda = L[0]*L[1]*L[2];

  return (double)N/VLambda;
}




void Get_PairCorr()
{

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //compute the minimum length of the lattice vectors
  double * g;
  
  g = new double [Nbin];
  for(int i = 0; i < Nbin; i++)
  {
		g[i] = 0.0;
   }	
  printf("Computing G2 and g2 now....\n");

  
  //loop over all particles
  for(int i = 0; i < N; i++)
      for(int j = 0; j < N; j++)
	{
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //for the radial g2....
	  double temp_dis = MinDis(j, i)/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(j != i && temp_dis < Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	      
	      g[int_dis] = g[int_dis] + 1.0;
	    }
	  }

	 
  double rho_n = Get_NDensity();
  
  int temp_sum = 0;

//This is still 2D.... need to generalize to corresponding dimenions...
  for(int r = 0; r < Nbin; r++)
    {
    	//cout<<"g_r = " <<g[r]<<endl;
    	
    	temp_sum += g[r];
    	
      g[r] = g[r]/(2*N*rho_n*4.0*pi*((r+0.5)*(r+0.5)*(r+0.5)-r*r*r)*Bin*Bin*Bin/3.0);   //rho*4pi*r^2*dr -> volume of a spherical shell with radius (r+0.5)*Bin
    }

 
 	//cout<<"sum_g = "<<temp_sum<<endl;
  
  FILE* fp = fopen("./g2_relax.xls","w");
  for(int r = 0; r < Nbin; r++)
    fprintf(fp, "%lf\t%lf\n", (r + 0.5) * Bin/nnl, g[r]);
  fclose(fp);
  delete [] g;
}


//this needs to be fixed as well

void num_var()
{
	printf("Computing number variance now....\n");
	int Ns = (int)floor(L[0] / 4.0 / Rbin) - 1;
	cout << "Ns = " << Ns << endl;
	double* sigma = new double [Ns];
	double ave, square_ave;
	int Nc;
	double r = Rbin;
	
	double dx[SD];
	double cx[SD];
	double dist;
	
	for(int t = 0; t < Ns; t++)
	{
		cout << "r = " << r << endl;
		ave = 0.0;
		square_ave = 0.0;
		for(int n = 0; n < Nr; n++)
		{
			Nc = 0;
			//double cx = distribution(gen) * Lx;
			//double cy = distribution(gen) * Ly;
			//cout << cx << " " << cy << endl;
			for(int i=0; i<SD; i++)
			   cx[i] = ((double)(rand() % MAXY) / (double)MAXY) * L[i];
			//double cy = ((double)(rand() % MAXY) / (double)MAXY) * Ly;
			
			for(int i = 0; i < N; i++)
			{
				dist = 0;
				
				for(int j=0; j<SD; j++)
				{
					dx[j] = fabs(Center[i][j] - cx[j]);
					if(dx[j] > L[j]/2.0)
					{
						dx[j] = L[j] - dx[j]; //periodic boundary conditions
					}
					
					dist += dx[j]*dx[j];
				}
				
				
				if(dist < r*r)
				{
					Nc++; 
				}
			}
//			cout << "Nc = " << Nc << endl;
			ave = ave + (double)Nc;
			square_ave = square_ave + (double)(Nc * Nc);
		}
		ave = ave / Nr;
		square_ave = square_ave / Nr;
		cout << "ave = " << ave << endl;
		cout << "square_ave = " << square_ave << endl;
		sigma[t] = square_ave - ave * ave;
		r = r + Rbin;
	}
	FILE* fp = fopen("num_var_relax.xls","w");
	for(int r = 0; r < Ns; r++)
		fprintf(fp, "%lf\t%lf\n", (r + 1) * Rbin, sigma[r]);
	fclose(fp);
	delete [] sigma;
}


double GetInnerProduct(double Vector1[SD], double Vector2[SD])
{
	double sum = 0;

	for(int i = 0; i < SD; i++)
		sum += Vector1[i] * Vector2[i];
	return sum;
}

double Get_Sk(double Vector_k[SD])
{
	double Sk;
	double sum_cos = 0.0;
	double sum_sin = 0.0;
	double temp;
	for(int i = 0; i < N; i++)
	{
			temp = GetInnerProduct(Center[i], Vector_k); 
			sum_cos += cos(temp);
			sum_sin += sin(temp);
	}
	Sk = (sum_cos * sum_cos + sum_sin * sum_sin) / (double)N;
	return Sk;
}

void Print_Sk(double K_Histo[Kbin_num], int K_Counter[Kbin_num])
{
	ofstream histo_out;
	histo_out.open("Sk_relax.txt");
	for(int t = 0; t < Kbin_num; t++)
		if(K_Counter[t] > 0)
			histo_out << (t + 0.5) * Kbin *(nnl/(2.0*pi)) << "\t" << K_Histo[t] << endl; //this should be normalized with a local length scale...
			//histo_out << (t + 0.5) * Kbin * mean_nnb_dist / 2.0 / pi << " " << K_Histo[t] << endl;
	histo_out.close();
	
	histo_out.open("Sk_relax.xls");
	for(int t = 0; t < Kbin_num; t++)
		if(K_Counter[t] > 0)
			histo_out << (t + 0.5) * Kbin *(nnl/(2.0*pi)) << "\t" << K_Histo[t] << endl; //this should be normalized with a local length scale...
			//histo_out << (t + 0.5) * Kbin * mean_nnb_dist / 2.0 / pi << " " << K_Histo[t] << endl;
	histo_out.close();
}


void Get_KHistogram()
{
	printf("Computing Sk now....\n");
	double Sk, k_dis;
	int t;
	double KPoint[SD];
	double K_Histo[Kbin_num];
	int K_Counter[Kbin_num];
	for(int i = 0; i < Kbin_num; i++)
	{
		K_Histo[i] = 0.0;
		K_Counter[i] = 0;
	}
	
	//here we onnly consider a portion of the full Fourier space, i.e., only the 1st quadrant, this should be OK for istropic systems.
	//due to the symmetry constraints, only half 
	
	cout<<"generating the k-space points...."<<endl;
	cout<<"Nk = "<<Nk<<endl;
	
	int temp_counter = 0;
	
	
	//we only look at this 1st quadent 
	for(int i = 1; i <= Nk; i++)
		for(int j = 1; j <= Nk; j++)
			for(int k = 1; k<=Nk; k++)
		{
			KPoint[0] = i * 2.0 * pi / L[0];
			KPoint[1] = j * 2.0 * pi / L[1];
			KPoint[2] = k * 2.0 * pi / L[2];
			//KPoint[3] = m * 2.0 * pi / L[3];
			
			Sk = Get_Sk(KPoint);
			k_dis = 0.0;
			for(int d = 0; d < SD; d++)
				k_dis += KPoint[d] * KPoint[d];
			k_dis = sqrt(k_dis);
			t = floor(k_dis / Kbin);
			if(t < Kbin_num)
			{
				K_Histo[t] += Sk;
				K_Counter[t] ++;
			}
			
			temp_counter++;
			cout<<temp_counter<<endl;

		}
		
	for(t = 0; t < Kbin_num; t++)
		if(K_Counter[t] != 0)		
			K_Histo[t] = K_Histo[t] / K_Counter[t];
	//cout << "first bin: " << K_Counter[0] << " " << K_Histo[0] << endl;


	cout<<"printing out the sk..."<<endl;
	Print_Sk(K_Histo, K_Counter);

}


void print_connectivity_matrix()
{
	//loop for all particles, and directly print out the adjancency matrix..
	
	ofstream fout;
	fout.open("adj_matrix.txt");
	
	double temp_dist;
	
	cout<<"please specify cut_off distance (in unit of D), cut_off = ";
	cin>>cut_off;
	
	cout<<"looping overall particles to generate the connectivity matrix..."<<endl;
	
	for(int i=0; i<N; i++) 
	{
		for(int j=0; j<N; j++)
			{
				temp_dist = MinDis(i,j);
				
				if(temp_dist<cut_off*nnl && i!=j)
					fout<<"1"<<"\t";
				else fout<<"0"<<"\t";
			}
			
		fout<<endl;
	}
	
	fout.close();
}


void get_contact_number()
{
	cout<<"computing the contact number for different gap tolerance..."<<endl; 
	
	int gap_pow = -10; //this is the smallest gap pow 
	
	int n_contact[8]; //the dimension gives the number of gap
	int num_gap = 8; //must be the same as the size of n_contact
		
	for(int i=0; i<num_gap; i++)
		n_contact[i] = 0;
		
	
	int cn_counter[N][8]; //the number of contact number for each particle, for rattler statistics
	for(int i=0; i<N; i++)
		for(int j=0; j<num_gap; j++)
			cn_counter[i][j] = 0;	

	
		
	double temp_dist;
		
	cout<<"looping over all particles to compute the distances now ..."<<endl;
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			temp_dist = MinDis(i, j);
			
			for(int m=0; m<num_gap; m++)
			{
				if(temp_dist < nnl*(1+pow(10, gap_pow+m))&&i!=j)
				{
					n_contact[m]++;
					cn_counter[i][m]++;
				}
					
			}
		}
	
	//count the number of spheres in the contact network
	int N_gap[8];
	
	for(int m=0; m<num_gap; m++)
	{
		N_gap[m] = 0;
		
		for(int i=0; i<N; i++)
			{
				if(cn_counter[i][m]>=(SD+1)) //condition for local jamming
					N_gap[m]++;
			}
	}
	
	
	
	fout.open("n_contact.xls");
	fout<<"the average contact does not exclude rattlers!"<<endl;
	for(int m=0; m<num_gap; m++)
		fout<<pow(10, gap_pow+m)<<"\t"<<((double)n_contact[m])/(double)N<<endl;
	fout.close();
	
	fout.open("fr_rattlers.xls");
	for(int m=0; m<num_gap; m++)
		fout<<pow(10, gap_pow+m)<<"\t"<<((double)(N-N_gap[m]))/(double)N<<endl;
	fout.close();
}


int main()
{
	read_packing_3D();
	
	print_connectivity_matrix();
	
	get_contact_number();
	
	Get_PairCorr();
	Get_KHistogram();
	//srand(time(NULL));
	
	//perturb_config();
	
	//print_config();
	
	num_var();
	
   
    
    //int t_a;
    //cin>>t_a;

	return 1;
}


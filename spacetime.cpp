//---------- This code was written by Sabiha Majumder in 2016--------------------//

#include <iostream>
#include <math.h>
using std::cerr;
using namespace std;
#include <fstream>
using std::ofstream;
#include <cstdlib>       // for exit function
#include <ctime>


// parameter definitions
     
const int N = 128;  // system size
const int p_range = 181 ;  // no. of p values in the transect 
const int T = 2000 ;  // No. of time steps

int Nx=N, Ny=N*p_range+4; //Nx=No. of rows and Ny= no. of columns. Extra 4 columns are alloted for the boundary conditions.

//====== function to initialize the matrix ===============//
// This function creates a binary matrix with Nx rows and Ny cloumns. Here x is the proportion of 1's

void create_random_matrix(int*U,float x) {  

	
	int i1, i2;
	
	for (i1 = 0; i1 < Nx; i1++) {
		for (i2 = 0; i2 < Ny; i2++) {
			double number = rand()/(double)RAND_MAX;
			if (number<x) U[i1*Ny+i2]=1;
			else U[i1*Ny+i2]=0;
			
		}
	}

}

//====================select random neighbouring site==============================//
void select_neighbor_of_site(int i ,int j , int*neighbor){
	int left,right,bottom,top ;
	int in = i ,jn = j;
	
	double test = rand()/(double)RAND_MAX ;
	
	if (i==0) top=Nx-1;
	else top = i-1;
	
	if (i==Nx-1) bottom = 0;
	else bottom = i+1 ;
	
	left = j-1 ;
	
	right = j+1;
				
	if (test <= 0.25) in = top;
	else if ( test <= 0.5) in = bottom;
	else if ( test <=0.75) jn =left;
	else jn = right;
	
	neighbor[0] = in;
	neighbor[1] = jn;
	
}

//========select neighbor of pair=====================================//
int select_neighbor_of_pair(int in ,int jn, int i, int j){
	int left,right,top,bottom,leftn,rightn,topn,bottomn , neighbor_of_pair;
	
	
	if (i==0) top=Nx-1;                 //periodic boundary
	else top = i-1;
	
	if (i==Nx-1) bottom = 0;
	else bottom = i+1 ;
	
	left = j-1 ;
	
	right = j+1;
	
	if (in==0) topn=Nx-1;
	else topn = in-1;
	
	if (in==Nx-1) bottomn = 0;
	else bottomn = in+1 ;
	
	leftn = jn-1 ;
	
	rightn = jn+1;
	
	int nn[6] ,c=0;
	
	if ((top*Ny +j) != (in*Ny+jn)) {
		nn[c]=top*Ny + j;
		c+=1;
	}
	if ((bottom*Ny + j) != (in*Ny+jn)) {
		nn[c]=bottom*Ny  + j;
		c+=1;
	}
	if ((i*Ny +right) != (in*Ny+jn)) {
		nn[c]= i*Ny + right;
		c+=1;
	}
	if ((i*Ny +left) != (in*Ny+jn)) {
		nn[c] = i*Ny + left;
		c+=1;
	}
	if ((topn*Ny+jn) != (i*Ny+j)) {
		nn[c]=topn*Ny + jn;
		c+=1;
	}
	if ((bottomn*Ny +jn) != (i*Ny+j)) {
		nn[c]=bottomn*Ny  + jn;
		c+=1;
	}
	if ((in*Ny +rightn) != (i*Ny+j)) {
		nn[c]= in*Ny + rightn;
		c+=1;
	}
	if ((in*Ny +leftn) != (i*Ny+j)) {
		nn[c] = in*Ny + leftn;
		c+=1;
	}
	
	double test =rand()/(double)RAND_MAX ;
	
	if (test <=(0.1666)) neighbor_of_pair= nn[0];
	else if ( test <= (2*0.1666)) neighbor_of_pair= nn[1];
	else if ( test <= (3*0.1666)) neighbor_of_pair= nn[2];
	else if ( test <= (4*0.1666)) neighbor_of_pair= nn[3];
	else if ( test <= (5*0.1666)) neighbor_of_pair= nn[4];
	else neighbor_of_pair = nn[5];
	

return neighbor_of_pair;

} 

		
////////////// main function //////////////////////////////////
int main(){
	
	srand(time(NULL));
	
	int x,l,t,i,j,z;
	double p[p_range], q=0.95;
	double  mean;
	int* neighbor = new int[2];
	
	p[0]=0.2400; //Value of p at the start of the transect
	for(x=1;x<p_range;++x) {
		p[x]= p[x-1] + 0.0001 ;
	}
	 
	int*A = new int[Nx*Ny];
	
	create_random_matrix(A,0.5);  //Random initial matrix
			
		for(t=0;t<T;++t){
		
			for(z=0;z<Nx*Ny ; ++z){                // so that each site gets selected once on an average
				
				i = rand()%Nx;           // selecting one random site
				j = rand()%Ny;
					
				if (j>1 && j<Ny-2){   // We dont consider the first two and last two columns for fixed boundary
					x = (j-2)/N ;    // p value for that cite is p[x]
					
					double test = rand()/(float)RAND_MAX;
					double test1 = rand()/(float)RAND_MAX;
					
					if (A[i*Ny+j]==1){     //if the site is occupied
						
						select_neighbor_of_site(i, j ,neighbor);    //look for a neighbor
						int in = neighbor[0] , jn = neighbor[1];
						
						if (A[in*Ny +jn]==0) {                     //if neighbor is empty
							if (test < p[x]) 
								A[in*Ny+jn]=1;                 //regular cp
							else A[i*Ny+j]=0;
						} 

						else {							  
							if (test < q){
								
								int neighbor_of_pair=select_neighbor_of_pair (in, jn, i, j);  //look for the neighbor of pair 
								A[neighbor_of_pair]=1;
							}
							else if (test1 < 1-p[x])
								A[i*Ny+j]=0;	
						}
					}
					}

				}	
	 // Now fixing the left boundary based on the density of the left matrix
	 
	 			//calculate mean of the left matrix
				x=0;
				mean=0;
				for (i=0;i<N; ++i){     
					for (j=2;j<N+2;++j){      //we are calculating average of NxN matrix for one p-value
						mean += A[i*Ny+j];						
					}
				}
				
				mean = mean/(N*N);
				
				//update first and second column based on the mean
				
				for (i=0; i<N; ++i){
					for (j=0;j<2; ++j){
						double number = rand()/(double)RAND_MAX;
						if (number<mean) A[i*Ny+j]=1;
						else A[i*Ny+j]=0;
					}
				}
				
				//calculate mean of the right matrix
				x=p_range-1;
				mean=0;
				for (i=0;i<N; ++i){     //since N=Nx
					for (j=x*N+2;j<x*N+N+2;++j){      //we are calculating average of NxN matrix for one p-value
						mean += A[i*Ny+j];
					}
				}
				
				mean = mean/(N*N);
				
				//update last and second last column based on the mean
				
				for (i=0; i<N; ++i){
					for (j=N*p_range+2;j<N*p_range+4; ++j){
						double number = rand()/(double)RAND_MAX;
						if (number<mean) A[i*Ny+j]=1;
						else A[i*Ny+j]=0;
					}
				}
				
						
				
				
			}
			
			//Saving data in a file
			
			ofstream outdata;
			ofstream snapshots;
			
			outdata.open("filename1.txt",ios::app);
			snapshots.open("filename2.txt",ios::app);
			
			
			for (x=0;x<p_range;++x){
				mean=0;
				for (i=0;i<N; ++i){     //since N=Nx
					for (j=x*N+2;j<x*N+N+2;++j){      //we are calculating average of NxN matrix for one p-value
						mean += A[i*Ny+j];
						snapshots<<A[i*Ny+j]<<endl;
						
					}
				}
				
				mean = mean/(N*N);
				
				outdata<<mean<<endl;
				cout<<p[x]<<"\t"<<mean<<endl;
				
			}
		
			outdata.close();
			snapshots.close();
			
			
	
	delete[] A;
	delete[] neighbor;
	return 0;
}

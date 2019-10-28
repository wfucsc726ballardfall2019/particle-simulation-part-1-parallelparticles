#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>
#include <omp.h>
using namespace std;


// void print_distance(){
//     cout << "Applying particle " << i << " to particle " << j;
//     double dx = particles[i].x - particles[j].x;
//     double dy = particles[i].y - particles[j].y;
//     double r2 = dx * dx + dy * dy;
//     cout << ", ditance = " << r2 << " cutoff = " << cutoff * cutoff << endl;
// }


//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0, head_index;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    double cutoff = .05;
    double space_dim = sqrt(n * 0.0005);
    cout << "Actual space is " << space_dim << " by " << space_dim << endl;
    // cutoff = space_dim / 
    double cell_edge = space_dim / floor(space_dim / cutoff);
    // cell_edge = space_dim / 4;
    cout << "Cells are size " << cell_edge << " by " << cell_edge << endl;
    int num_cells = pow((space_dim / cell_edge), 2);
    // cout << num_cells << endl;
    int x_cells = (int)sqrt(num_cells);
    cout << "Using grid size " << x_cells << " by " << x_cells << endl;


    int lscl[n]; //might need to check initial values, paper has them all -1s
    int head[num_cells];

    double r[n][2];
    double mc[2];
    int mc1[2];

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	bool flag = false;
    
    #pragma omp parallel private(dmin) num_threads(4)
	{int numthreads = omp_get_num_threads();
		cout << "numthreads= " <<numthreads << endl;
	    for( int step = 0; step < NSTEPS; step++ )
  		  {
    		    // cout << endl << endl;
      		    navg = 0;
      		    davg = 0.0;
      	            dmin = 1.0;
        // for(int i = 0; i < n; i++){
        //     r[i][0] = particles[i].x;
        //     r[i][1] = particles[i].y;
        // }

       		    for(int c = 0; c < num_cells; c++) head[c] = -1;
       			#pragma omp parallel for private(mc, head_index)
			 for(int i = 0; i < n; i++){
            // cout << "initial point: " << particles[i].x << ", " << particles[i].y << endl;
            // for(int a = 0; a < 2; a++) mc[a] = r[i][a] / cell_edge;
          		  mc[0] = floor(particles[i].x / cell_edge);
           		  mc[1] = floor(particles[i].y / cell_edge);
           		  head_index = mc[0] * x_cells + mc[1];
             if(flag == false){
                 cout << "mc: ";
                 for(int i = 0; i < 2; i++){
                     cout << mc[i] << " ";
                 }
                 cout << " cell index: " << head_index << endl;
             }
            lscl[i] = head[head_index];
            head[head_index] = i;
        }
        // if(flag == false){
            // cout << "Head Array: ";
            // for(int i = 0; i < num_cells; i++){
            //     cout << head[i] << " ";
            // }
            // cout << endl << "LSCL Array: ";
            // for(int i = 0; i < n; i++){
            //     cout << lscl[i] << " ";
            // }
            // cout << endl;
        // }
        // flag = true;
    
       
	 for(mc[0] = 0; mc[0] < x_cells; (mc[0])++){
            for(mc[1] = 0; mc[1] < x_cells; (mc[1])++){
                int c = mc[0] * x_cells + mc[1];
               // cout << "Vector index: " << mc[0] << ", " << mc[1] << " scalar index: " << c << endl;
                for(mc1[0]=mc[0]-1; mc1[0] <= mc[0]+1; (mc1[0])++){
                    for(mc1[1]=mc[1]-1; mc1[1] <= mc[1]+1; (mc1[1])++) {
                        // cout << "\t Neighbor index: " << mc1[0] << ", " << mc1[1] << endl;
                        bool test = false;
                        int j;
                        int i;
                        for(int a=0; a<2; a++) {
                            if(mc1[a] < 0 || mc1[a] >= x_cells){
                                // cout << "\tEDGE CASE FOUND: " << mc1[0] << ", " << mc1[1] << endl;
                                test = true;
                            }
                        }
                        int c1 = ((mc1[0]+x_cells)%x_cells)*x_cells + ((mc1[1]+x_cells)%x_cells);
                         //cout << " scalar index: " << c1 << endl;
                        i = head[c];
                        while(i != -1){
                            if(test == false){
                                j = head[c1];
                            }
                            else{
                                j = -1;
                            }
                          
			 	
                            while(j != -1){
                                if(i < j){
                                    // cout << "Applying particle " << i << " to particle " << j;
                                     //double dx = particles[i].x - particles[j].x;
                                    // double dy = particles[i].y - particles[j].y;
                                    // double r2 = dx * dx + dy * dy;
                                    // cout << ", ditance = " << r2 << " cutoff = " << cutoff * cutoff << endl;
                                    apply_force( particles[i], particles[j],&dmin,&davg,&navg);
                                }
                                j = lscl[j];
                            }
                            i = lscl[i];
                            }
                        }
                         //cout << "Vector index: " << mc1[0] << ", " << mc1[1] << " scalar index: " << endl;
                    }
                }
            }
        


        // navg = 0;
        //     davg = 0.0;
        // dmin = 1.0;
        //
        //  compute forces
        //
        // cout << endl << endl << endl << "NAIVE ALGORITHM" << endl;
        // for( int i = 0; i < n; i++ )
        // {
        //     particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ ){
        //         cout << "Applying particle " << i << " to particle " << j;
        //         double dx = particles[i].x - particles[j].x;
        //         double dy = particles[i].y - particles[j].y;
        //         double r2 = dx * dx + dy * dy;
        //         cout << ", ditance = " << r2 << " cutoff = " << cutoff * cutoff << endl;
        //         apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        //     }
        // }
 
        //
        //  move particles
        //
 	#pragma omp parallel for 
       	 for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

      	  if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          #pragma omp master
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }

	  #pragma omp critical	
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //

	  #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}

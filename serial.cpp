#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>
#include <vector>
using namespace std;


//
//  benchmarking program
//
int main( int argc, char **argv )
{     
    int navg,nabsavg=0;
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

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    

    double cutoff = 0.02;
    double space_dim = sqrt(n * 0.0005);
    // cout << "Actual space is " << space_dim << " by " << space_dim << endl;
    double cell_edge = space_dim / floor(space_dim / cutoff);
    // cell_edge = space_dim / 4;
    // cout << "Cells are size " << cell_edge << " by " << cell_edge << endl;
    int num_cells = pow((space_dim / cell_edge), 2);
    // cout << num_cells << endl;
    int cells_in_row = (int)sqrt(num_cells);
    // cout << "Using grid size " << x_cells << " by " << x_cells << endl;


    vector<vector<int> > cell_vector(num_cells);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        // cout << endl;
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute forces
        //
        for(int i = 0; i < num_cells; i++){
            cell_vector[i].clear();
        }
        for(int i = 0; i < n; i++){
            // cout << particles[i].x << " " << particles[i].y << " ";
            int cell_x = floor(particles[i].x / cell_edge);
            int cell_y = floor(particles[i].y / cell_edge);
            int cell_index = cell_x * cells_in_row + cell_y;
            // cout << "--> (" << cell_x << ", " << cell_y << ") --> " << cell_index << endl;
            cell_vector[cell_index].push_back(i);
        }
        // for(int i = 0; i < num_cells; i++){
        //     cout << "Cell " << i << ": ";
        //     for(int j = 0; j < cell_vector[i].size(); j++){
        //         cout << cell_vector[i][j] << " ";
        //     }
        //     cout << endl;
        // }

        for(int cell_x = 0; cell_x < cells_in_row; cell_x++){
            for(int cell_y = 0; cell_y < cells_in_row; cell_y++){
                int cell_index = cell_x * cells_in_row + cell_y;
                for (int p = 0; p < cell_vector[cell_index].size(); p++){
                    particles[cell_vector[cell_index][p]].ax = particles[cell_vector[cell_index][p]].ay = 0;
                }
                // cout << "Cell index " << cell_x << ", " << cell_y << " (" << cell_index << ")" << endl;

                for(int neighbor_x = cell_x - 1; neighbor_x <= cell_x + 1; neighbor_x++){
                    for(int neighbor_y = cell_y - 1; neighbor_y <= cell_y + 1; neighbor_y++){

                        if((neighbor_x >=0 && neighbor_x <= cells_in_row - 1) && (neighbor_y >=0 && neighbor_y <= cells_in_row - 1)){
                            int neighbor_index = ((neighbor_x + cells_in_row) % cells_in_row) * cells_in_row + ((neighbor_y + cells_in_row) % cells_in_row);
                            // cout << "\tNeighbor " << neighbor_x << ", " << neighbor_y << " (" << neighbor_index << ")" << endl;

                            for(int cell_p = 0; cell_p < cell_vector[cell_index].size(); cell_p++){
                                // particles[cell_vector[cell_index][l]].ax = particles[cell_vector[cell_index][l]].ay = 0;
                                for(int neighbor_p = 0; neighbor_p < cell_vector[neighbor_index].size(); neighbor_p++){
                                    // if(cell_vector[cell_index][l] != cell_vector[neighbor_index][m]){
                                        // cout << "\t\tApplying " << cell_vector[cell_index][l] << " (" << particles[cell_vector[cell_index][l]].ax << ", " << particles[cell_vector[cell_index][l]].ay << ") " << " to " << cell_vector[neighbor_index][m] << " (" << particles[cell_vector[neighbor_index][m]].ax << ", " << particles[cell_vector[neighbor_index][m]].ay << ") " << endl;
                                        apply_force( particles[cell_vector[cell_index][cell_p]], particles[cell_vector[neighbor_index][neighbor_p]],&dmin,&davg,&navg);
                                    // }
                                }
                            }
                        }
                    }
                }
            }
        }
            // for(int j = i - cells_in_row; j <= i + cells_in_row; j += cells_in_row){
            //     for(int k = j - 1; k <= j + 1; k++){
            //         if(k >=0 && k <= num_cells - 1){
            //             cout << "\tNeighbor " << k << endl;
            //             for(int l = 0; l < cell_vector[i].size(); l++){
            //                 particles[cell_vector[i][l]].ax = particles[cell_vector[i][l]].ay = 0;
            //                 for(int m = 0; m < cell_vector[k].size(); m++){
            //                     // cout << l << " " << m << endl;
            //                     if(cell_vector[i][l] < cell_vector[k][m]){
                                    
            //                     }
            //                 }
            //             }
            //         }
            //     }
            // }
        


        // for( int i = 0; i < n; i++ )
        // {
        //     particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ ){
        //         // if(!flag){
        //             cout << "Applying " << i << " (" << particles[i].ax << ", " << particles[i].ay << ") to " << j << " (" << particles[j].ax << ", " << particles[j].ay << ")" << endl;
        //         // }
		// 		apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        //     }
        //     flag = true;
        // }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
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
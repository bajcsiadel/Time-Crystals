//
//  globaldata.h
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef globaldata_h
#define globaldata_h

#include <stdio.h>
#include <time.h>
#include "jsmn.h"

struct global_struct
{
    double SX, SY;
    double halfSX, halfSY;
    
    // number of pinningsites
    int     N_pinningsites;
    // x coordinate of pinningsite
    double  *pinningsite_x;
    // y coordinate of pinningsite
    double  *pinningsite_y;
    // Total force affecting on pinningsite in direction x
    double  *pinningsite_fx;
    // Total force affecting on pinningsite in direction y
    double  *pinningsite_fy;
    // Color of pinningsite
    int     *pinningsite_color;
	// Scale of external force affecting on pinningsite in direction x
    double  *pinningsite_direction_x;
	// Scale of external force affecting on pinningsite in direction y
    double  *pinningsite_direction_y;
	// Travelled distance in direction x since the Verlet list was rebuilt
    double  *pinningsite_dx_so_far;
	// Travelled distance in direction y since the Verlet list was rebuilt
    double  *pinningsite_dy_so_far;
    // Radius of pinningsite
    double  pinningsite_R;
    // Internal force affecting in pinningsite
    double  pinningsite_force;
    // Number of particles in the pinningsite
    unsigned int *particle_in_pinningsite;
    
    // Pinningsite lattice width
    double pinning_lattice_constant;
    // External force affecting on pinningsite
    double pinning_driving_force;
    int pinning_direction_change;
    
    // number of particles
    int     N_particles;
	// x coordinate of particles
    double  *particle_x;
	// y coordinate of particles
    double  *particle_y;
	// Total force affecting on particle in x direction
    double  *particle_fx;
	// Total force affecting on particle in y direction
    double  *particle_fy;
	// Color of particle
    int     *particle_color;
	// Scale of external force affecting on particle in direction x
    double  *particle_direction_x;
	// Scale of external force affecting on particle in direction y
    double  *particle_direction_y;
	// Travelled distance in direction x since the Verlet list was rebuilt
    double  *particle_dx_so_far;
	// Travelled distance in direction y since the Verlet list was rebuilt
    double  *particle_dy_so_far;
	// Travelled distance in the whole simulation in direction x
    double  *particle_all_dx;
	// Travelled distance in the whole simulation in direction y
    double  *particle_all_dy;
	// Travelled distance in the whole simulation
    double  *particle_all_dr2;

    // External force affecting on particle
    double  particle_driving_force;
    double  partile_particle_screening_length;
    double  partile_particle_screening_wavevector;
    
    int N_Verlet;
    int N_Verlet_max; //initial allocation + later, longest allocation
    int *Verletlisti;
    int *Verletlistj;
    int flag_to_rebuild_Verlet;
    
    double pinningsite_grid_dx;
    double pinningsite_grid_dy;
    unsigned int Nx_pinningsite_grid;
    unsigned int Ny_pinningsite_grid;
    int ***pinningsite_grid;
    unsigned int max_pinningsite_per_grid;
    
    double Verlet_cutoff_distance;
    double Verlet_cutoff_distance_squared;
    double Verlet_intershell_squared;
    
    double dt;
    
    unsigned int total_time;
    unsigned int echo_time;
    unsigned int movie_time;
    unsigned int statistics_time;
    unsigned int time;

    FILE *moviefile;
    char *moviefile_name;
    FILE *statisticsfile;
    char *statisticsfile_name;

    jsmntok_t* t;
    char* JSON;

    time_t start_time;
    time_t end_time;
};

struct flag_struct
{
    short int system_size_SX_set;
    short int system_size_SY_set;
};



extern struct global_struct global;
extern struct flag_struct flag;
#endif /* globaldata_h */

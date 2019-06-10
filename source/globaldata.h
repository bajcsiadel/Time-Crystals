//
//  globaldata.h
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef globaldata_h
#    define globaldata_h

#    include <stdio.h>
#    include <time.h>
#    include "jsmn.h"

typedef struct global_struct
{
	double  SX, SY;
	double  halfSX, halfSY;

	// number of pinningsites
	int     N_pinningsites;
	// x coordinate of pinningsite
	double *pinningsite_x;
	// y coordinate of pinningsite
	double *pinningsite_y;
	// Total force affecting on pinningsite in direction x
	double *pinningsite_fx;
	// Total force affecting on pinningsite in direction y
	double *pinningsite_fy;
	// Color of pinningsite
	int    *pinningsite_color;
	// Scale of external force affecting on pinningsite in direction x
	double *pinningsite_direction_x;
	// Scale of external force affecting on pinningsite in direction y
	double *pinningsite_direction_y;
	// Travelled distance in direction x since the Verlet list was rebuilt
	double *pinningsite_dx_so_far;
	// Travelled distance in direction y since the Verlet list was rebuilt
	double *pinningsite_dy_so_far;
	// Radius of pinningsite
	double  pinningsite_R;
	// Internal force affecting in pinningsite
	double  pinningsite_force;
	// Number of particles in the pinningsite
	unsigned int *particle_in_pinningsite;

	// Pinningsite lattice width
	double  pinning_lattice_constant;
	// External force affecting on pinningsite
	double  pinning_driving_force;
	int     pinning_direction_change;

	// number of particles
	int     N_particles;
	// x coordinate of particles
	double *particle_x;
	// y coordinate of particles
	double *particle_y;
	// Total force affecting on particle in x direction
	double *particle_fx;
	// Total force affecting on particle in y direction
	double *particle_fy;
	// Color of particle
	int    *particle_color;
	// Scale of external force affecting on particle in direction x
	double *particle_direction_x;
	// Scale of external force affecting on particle in direction y
	double *particle_direction_y;
	// Travelled distance in direction x since the Verlet list was rebuilt
	double *particle_dx_so_far;
	// Travelled distance in direction y since the Verlet list was rebuilt
	double *particle_dy_so_far;
	// Travelled distance in the whole simulation in direction x
	double *particle_all_dx;
	// Travelled distance in the whole simulation in direction y
	double *particle_all_dy;
	// Travelled distance in the whole simulation
	double *particle_all_dr2;

	// External force affecting on particle
	double  particle_driving_force;
	double  partile_particle_screening_length;
	double  partile_particle_screening_wavevector;

	// Number of particle pairs in the Verlet list
	int     N_Verlet;
	// Initial allocation + later, longest allocation
	int     N_Verlet_max;
	// Particle ID
	int    *Verletlisti;
	// Particle ID
	int    *Verletlistj;
	// Flag if it is needed to rebuild Verlet list or not
	int     flag_to_rebuild_Verlet;

	// Number of grids in direction x
	unsigned int Nx_pinningsite_grid;
	// Number of grids in direction y
	unsigned int Ny_pinningsite_grid;
	// Width of a grid cell
	double  pinningsite_grid_dx;
	// Height of a grid cell
	double  pinningsite_grid_dy;
	// 3-dimensional matrix containing the pinningsite grid. First two dimensions describe the grid and the third dimension contains pinningsite ID in the current cell
	int  ***pinningsite_grid;
	// Maximum number of pinningsite enabled to be in a cell
	unsigned int max_pinningsite_per_grid;

	double  Verlet_cutoff_distance;
	double  Verlet_cutoff_distance_squared;
	double  Verlet_intershell_squared;

	// Time step between iterations in the system
	double  dt;

	// Total iteration number
	unsigned int total_time;
	// Echo data to screen
	unsigned int echo_time;
	// Write data to movie file after each movie_time iterations
	unsigned int movie_time;
	// Write data to statistics file after each statistics_time iterations
	unsigned int statistics_time;
	// Current interation
	unsigned int time;

	// Movie file descriptor
	FILE   *moviefile;
	// Path to movie file
	char   *moviefile_name;
	// Statistics file descriptor
	FILE   *statisticsfile;
	// Path to statistics file
	char   *statisticsfile_name;

	// Token array from JSON object
	jsmntok_t *t;
	// String containing the JSON object
	char   *JSON;

	// Beginnin of the simulation in miliseconds
	time_t  start_time;
	// End of the simulation in miliseconds
	time_t  end_time;
} global_struct;

typedef struct flag_struct
{
	short int system_size_SX_set;
	short int system_size_SY_set;
} flag_struct;



extern global_struct global;
extern flag_struct flag;
#endif /* globaldata_h */

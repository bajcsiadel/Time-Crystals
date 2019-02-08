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
    double halfSX,halfSY;
    
    int N_pinningsites;
    double *pinningsite_x;
    double *pinningsite_y;
    double *pinningsite_fx;
    double *pinningsite_fy;
    int *pinningsite_color;
    double *pinningsite_direction_x;
    double *pinningsite_direction_y;
    double *pinningsite_dx_so_far;
    double *pinningsite_dy_so_far;
    double pinningsite_R;
    double pinningsite_force;
    unsigned int *particle_in_pinningsite;
    
    double pinning_lattice_constant;
    double pinning_driving_force;
    int pinning_direction_change;
    
    int N_particles;
    double *particle_x;
    double *particle_y;
    double *particle_fx;
    double *particle_fy;
    int *particle_color;
    double *particle_direction_x;
    double *particle_direction_y;
    double *particle_dx_so_far;
    double *particle_dy_so_far;
    double *particle_all_dx;
    double *particle_all_dy;

    double particle_driving_force;
    double partile_particle_screening_length;
    double partile_particle_screening_wavevector;
    
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

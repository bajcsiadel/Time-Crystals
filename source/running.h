//
//  running.h
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef running_h
#define running_h

#include <stdio.h>

void run_simulation(void);

void calculate_external_forces_on_pinningsites(void);
void move_pinningsites(void);
void fold_pinningsite_back_PBC(int i);
void rebuild_pinning_grid(void);

void calculate_external_forces_on_particles(void);
void calculate_pairwise_forces(void);
void calculate_pinning_force_on_particles(void);
void check_Verlet_rebuild_condition_and_set_flag(void);
void rebuild_Verlet_list(void);
void move_particles(void);
void fold_particle_back_PBC(int i);

void distance_squared_folded_PBC(double x0,double y0,double x1,double y1,
        double *r2_return, double *dx_return,double *dy_return);

void write_cmovie_frame(void);

void adjust_pinningsite_directions(void);
void rotate_pinningsite_directions(void);

void calculate_statistics(void);
void write_statistics(void);

//functions for testing the program
void test_program_by_coloring(void);

void delete_arrays();

#endif /* running_h */

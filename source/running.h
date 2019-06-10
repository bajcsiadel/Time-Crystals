//
//  running.h
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef running_h
#    define running_h

#    include <stdio.h>
#    include <time.h>

/**
 * Running the simulation with the previously setted parameters.
*/
void    run_simulation();

/**
 * Adding external force affecting on pinningsites to total forces that affects it.
*/
void    calculate_external_forces_on_pinningsites();
/**
 * Moves the pinning sites one time step.
*/
void    move_pinningsites();
/**
 * Fold back an element from system (pinningsite or particle) into the PBC simulation box. Assumes it did not jump more than a box length. If it did the simulation is already broken anyhow
 * (x, y)	-coordinates of the element
*/
void    fold_back_PBC(double *, double *);
/**
 * Rebuild pinningsite grid. At first call allocating a 3-dimensional matrix, represention the pinning grid.
 * 1-D	- rows of pinningsite grid
 * 2-D	- columns of pinningsite grid
 * 3-D	- pinningsite ID which is in the given grid
 * With this it will be easier and faster to calculate pinningsite force on particles, because for a particle is enough to look at pinninsites in the neighbour cells.
*/
void    rebuild_pinning_grid();

/**
 * Adding external force affecting on particles to total forces that affects on it.
*/
void    calculate_external_forces_on_particles();
/**
 * Calculate pairwise force on particles, because they have a tossing force between them. For calculating this value we use the Verlet list and by using it we can accelerete the calculations.
*/
void    calculate_pairwise_forces();
/**
 * Pinningsites have an internal force that points to its center. This function calculates and adds the pinningsites force on particle using pinningsite grid.
*/
void    calculate_pinning_force_on_particles();
/**
 * Pinningsites have an internal force that points to its center. This function calculates and adds the pinningsites force on particle by iterating through each pinningsites and checking its distance from every particle, therefore this method costs N_pinningstie * N_partile iterations.
*/
void    calclulate_pinning_force_on_particles_without_grid();
/**
 * Check if any particle moved (potentially) enough to enter the inner Verlet shell coming from the outside, if one did, we have to rebuild the Verlet lists.
*/
void    check_Verlet_rebuild_condition_and_set_flag();
/**
 * Rebuilding the Verlet list. When it is called for the first time it allocates the arrays.
*/
void    rebuild_Verlet_list();
/**
 * Moves the particles one time step.
*/
void    move_particles();

/**
 * Calculates the shortest distance squared between 2 points in a PBC configuration. This is squared because we want to save on sqrt with the lookup table. Also used by the Verlet rebuild flag check where I check how much a particle moved.
 * (x0, y0)	- cootdinates of first point
 * (x1, y1)	- coordinates of second point
 * r2_return- distance square between the points
 * dx_return- distance between the points on axis x
 * dy_return- distance between the points on axis y
 * 
*/
void    distance_squared_folded_PBC(double, double, double,
									double, double *, double *, double *);

/**
 * Writing current state of the simulation into movie file as its frame.
*/
void    write_cmovie_frame();

/**
 * Change the direction of the pinning sites they will move on a triangle
*/
void    adjust_pinningsite_directions();
/**
 * Change the direction of the pinning sites they will move on a circle
*/
void    rotate_pinningsite_directions();

/**
 * Calculate values needed to write to statistics file.
*/
void    calculate_statistics();
/**
 * Write data into statistics file.
*/
void    write_statistics();

/**
 * Functions for testing the program
*/
void    test_program_by_coloring();

/**
 * At the end of the simulation writes to the screent the starting, ending, and total time (difference between end and start) the simulation ran.
*/
void    write_time();
/**
 * Free array used during the simulation.
*/
void    delete_arrays();

#endif /* running_h */

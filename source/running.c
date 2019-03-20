//
//  running.c
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#include "running.h"
#include "globaldata.h"
#include "color.h"

#include <stdlib.h>
#include <math.h>

#define PI 3.1415926535

//runs the simulation for the required steps
void run_simulation()
{
    rebuild_Verlet_list(); //build for the first time
    rebuild_pinning_grid(); //build for the first time

    time(&global.start_time);
    for (global.time = 0; global.time <= global.total_time; global.time++)
    {
        //adjust_pinningsite_directions();
        //rotate_pinningsite_directions();
        calculate_external_forces_on_pinningsites();
        move_pinningsites();
        // rebuild_pinning_grid();
        
        // calculate_external_forces_on_particles();
        calculate_pairwise_forces();
        // calculate_pinning_force_on_particles();
        calclulate_pinning_force_on_particles_without_grid();

        //this is where we calculate statistics
        calculate_statistics();
        move_particles();
        
        check_Verlet_rebuild_condition_and_set_flag();
        
        //one particle moved enough to rebuild the Verlet list
        if (global.flag_to_rebuild_Verlet == 1)
            rebuild_Verlet_list();
        
        //echo time
        if (global.time % global.echo_time == 0)
        {
            printf("%d / %d\n", global.time, global.total_time);
            fflush(stdout);
        }

        //statistics time
        if (global.time % global.statistics_time == 0 && global.time!=0)
            write_statistics();

        //movie write time
        if (global.time % global.movie_time == 0)
            write_cmovie_frame();  
    }
    delete_arrays();
    
    fclose(global.moviefile);
    fclose(global.statisticsfile);

    time(&global.end_time);
    write_time();
}

//change the direction of the pinning sites
//they will move on a circle
void rotate_pinningsite_directions()
{
    double theta, sint, cost;
    double newx, newy;
    int i;

    //theta  = 3.1415/180.0*0.005;
    theta = global.pinning_driving_force * global.dt * 0.2857;
    sint = sin(theta);
    cost = cos(theta);

    //printf("%lf %lf\n", sint, cost);

    //rotation matrix
    //cos(theta)*x -sin(theta)*y
    //sin(theta)*x +cos(theta)*y

    for(i=0;i<global.N_particles;i++)
    {
        newx = global.particle_direction_x[i] * cost - global.particle_direction_y[i] * sint;
        newy = global.particle_direction_x[i] * sint + global.particle_direction_y[i] * cost;
        
        global.particle_direction_x[i] = newx;
        global.particle_direction_y[i] = newy;
    }

    printf("%lf %lf\n", global.particle_direction_x[i], global.particle_direction_y[i]);

}

//change the direction of the pinning sites
//they will move on a triangle
void adjust_pinningsite_directions()
{
    int i, k;

    if ( (global.pinningsite_dx_so_far[0]*global.pinningsite_dx_so_far[0] +
        global.pinningsite_dy_so_far[0]*global.pinningsite_dy_so_far[0])
        >= global.pinning_lattice_constant * global.pinning_lattice_constant * 1.0)
    {
        //change directions now
        
        for (k = 0; k < global.N_pinningsites; k++)
            switch ( (global.pinningsite_color[k] - 4 + global.pinning_direction_change) % 3)
            {
                case 0: {
                        global.pinningsite_direction_x[k] = - 0.5;
                        global.pinningsite_direction_y[k] = - sqrt(3)/2.0;
                        break;
                        }
                case 1: {
                        global.pinningsite_direction_x[k] = 1.0;
                        global.pinningsite_direction_y[k] = 0.0;
                        break;
                        }
                case 2: {
                        global.pinningsite_direction_x[k] = - 0.5;
                        global.pinningsite_direction_y[k] = + sqrt(3)/2.0;
                        break;
                        }
            }
        global.pinning_direction_change++;
        for (i = 0; i < global.N_pinningsites; i++)
        {
            global.pinningsite_dx_so_far[i] = 0.0;
            global.pinningsite_dy_so_far[i] = 0.0;
        }
    }
}

void calculate_external_forces_on_pinningsites()
{
    int i;

    for(i = 0; i < global.N_pinningsites; i++)
    {
        global.pinningsite_fx[i] += global.pinningsite_direction_x[i] * global. pinning_driving_force;
        global.pinningsite_fy[i] += global.pinningsite_direction_y[i] * global. pinning_driving_force;
    }
}



//moves the pinning sites one time step
void move_pinningsites()
{

    int i;
    double dx, dy;

    for(i = 0; i < global.N_pinningsites; i++)
    {
        dx = global.pinningsite_fx[i] * global.dt;
        dy = global.pinningsite_fy[i] * global.dt;
        
        global.pinningsite_x[i] += dx;
        global.pinningsite_y[i] += dy;
        
        global.pinningsite_dx_so_far[i] += dx;
        global.pinningsite_dy_so_far[i] += dy;
        
        fold_pinningsite_back_PBC(i);
        
        global.pinningsite_fx[i] = 0.0;
        global.pinningsite_fy[i] = 0.0;
    }
}

void fold_pinningsite_back_PBC(int i)
{
    //fold back the pinningsite into the pBC simulation box
    //assumes it did not jump more thana  box length
    //if it did the simulation is already broken anyhow

    if (global.pinningsite_x[i] < 0) global.pinningsite_x[i] += global.SX;
    if (global.pinningsite_y[i] < 0) global.pinningsite_y[i] += global.SY;
    if (global.pinningsite_x[i] >= global.SX) global.pinningsite_x[i] -= global.SX;
    if (global.pinningsite_y[i] >= global.SY) global.pinningsite_y[i] -= global.SY;
}

void calculate_external_forces_on_particles()
{
    int i;

    for(i = 0; i < global.N_particles; i++)
    {
        global.particle_fx[i] += global.particle_direction_x[i] * global.particle_driving_force;
        global.particle_fy[i] += global.particle_direction_y[i] * global.particle_driving_force;
        //printf("ext_fx= %lf\n", global.particle_direction_x[i] * global. particle_driving_force);
    }
}

void calculate_pairwise_forces()
{
    int ii, i, j;
    double r, r2, f;
    double dx, dy;

    for(ii = 0; ii < global.N_Verlet; ii++)
    {
        //obtain the i, j from the Verlet list
        i = global.Verletlisti[ii];
        j = global.Verletlistj[ii];
        
        //perform the pairwise force calculation
        distance_squared_folded_PBC(global.particle_x[i], global.particle_y[i], 
                global.particle_x[j], global.particle_y[j], &r2, &dx, &dy);
        
        //non-tabulated version
        //try to not divide just multiply division is costly
        r = sqrt(r2);
        if (r < 0.2)
        {
            //printf("WARNING:PARTICLES TOO CLOSE. LOWER FORCE USED\n");
            f = 100.0;
        }
        else
        {
            f = 1 / r2 * exp(-r * global.partile_particle_screening_wavevector);
            //if (r<1.0) printf("r=%lf f=%lf\n", r, f);
        }
            
        //division and multiplication for projection to the x, y axes
        
        f = f / r;
        //if (r<1.0) printf("%f/r=lf fx=%lf fy=%lf\n\n", f, f*dx, f*dy);
        
        global.particle_fx[i] -= f * dx;
        global.particle_fy[i] -= f * dy;
        global.particle_fx[j] += f * dx;
        global.particle_fy[j] += f * dy;
    }
}

void calculate_pinning_force_on_particles()
{
    unsigned int i, j, k, n, m;
    int gi, gj, gj2;
    float dr2, dx, dy;
    float r, x, y;

    r = global.pinningsite_R;
    for (i = 0; i < global.N_particles; i++)
    {
        // getting particle's position
        x = global.particle_x[i];
        y = global.particle_y[i];

        gi = (int) (x / global.pinningsite_grid_dx) - 1;
        gj = (int) (y / global.pinningsite_grid_dy) - 1;

        for (m = 0; m < 3; m++)
        {
            if (gi == -1) gi = global.Nx_pinningsite_grid - 1;
            else if (gi == global.Nx_pinningsite_grid - 1) gi = 0;
            else gi ++;

            gj2 = gj;
            for (n = 0; n < 3; n++)
            {
                if (gj2 == -1) gj2 = global.Ny_pinningsite_grid - 1;
                else if (gj2 == global.Ny_pinningsite_grid - 1) gj2 = 0;
                else gj2 ++;

                j = 0;
                while (j < global.max_pinningsite_per_grid && global.pinningsite_grid[gi][gj2][j] != -1)
                {
                    k = global.pinningsite_grid[gi][gj2][j];
                    distance_squared_folded_PBC(x, y, global.pinningsite_x[k], global.pinningsite_y[k], &dr2, &dx, &dy);
                    if (dr2 < r * r)
                    {
                        global.particle_fx[i] += dx / r * global.pinningsite_force;
                        global.particle_fy[i] += dy / r * global.pinningsite_force;
                        global.particle_in_pinningsite[k] ++;
                    }
                    j ++;
                }
            }
        }
    }
}

void calclulate_pinning_force_on_particles_without_grid()
{
    unsigned int i, j, k, n, m;
    double gi, gj, gj2;
    double dr2, dx, dy;
    double r, x, y;

    r = global.pinningsite_R;
    for (j = 0; j < global.N_pinningsites; j++)
    {
        gi = global.pinningsite_x[j];
        gj = global.pinningsite_y[j];
        for (i = 0; i < global.N_particles; i++)
        {
            // getting particle's position
            x = global.particle_x[i];
            y = global.particle_y[i];

            distance_squared_folded_PBC(gi, gj, x, y, &dr2, &dx, &dy);
            if (dr2 < r * r)
            {
                global.particle_fx[i] -= dx / r * global.pinningsite_force;
                global.particle_fy[i] -= dy / r * global.pinningsite_force;
                global.particle_in_pinningsite[j] ++;
            }
        }
    }
}

void check_Verlet_rebuild_condition_and_set_flag()
{
    int i;
    double dr2;

    //check if any particle moved (potentially) enough to enter the inner Verlet shell
    //coming from the outside, if one did, rebuild the Verlet lists
    global.flag_to_rebuild_Verlet = 0;

    for (i = 0; i < global.N_particles; i++)
    {
        dr2 = global.particle_dx_so_far[i] * global.particle_dx_so_far[i] +
              global.particle_dy_so_far[i] * global.particle_dy_so_far[i];
            
        //if (i==0) printf("%d %lf\n", global.time, dr2);fflush(stdout);
            
        if (dr2 >= global.Verlet_intershell_squared)
        {
            global.flag_to_rebuild_Verlet = 1;
            break; //exit the for cycle
        }
    }
}

//rebuilds the Verlet list
void rebuild_Verlet_list()
{
    int i, j;
    double dr2, dx, dy;
    double estimation;

    //initialize the Verlet list for the first time
    if (global.N_Verlet_max == 0)
    {
        //we are building the Verlet list for the first time in the simulation
        printf("Verlet list will be built for the first time\n");
        
        estimation = global.N_particles / (double)global.SX / (double)global.SY;
        printf("System density is %.3lf\n", estimation);
        
        estimation *= PI * global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
        printf("Particles in a R = %.2lf shell = %lf\n", 
            global.Verlet_cutoff_distance, estimation);
        
        global.N_Verlet_max = (int) estimation * global.N_particles / 2;
        printf("Estimated N_Verlet_max = %d\n", global.N_Verlet_max);
        
        global.Verletlisti = (int *) malloc(global.N_Verlet_max * sizeof(int));
        global.Verletlistj = (int *) malloc(global.N_Verlet_max * sizeof(int));
    }

    //build the Verlet list
    global.N_Verlet = 0;

    for (i = 0; i < global.N_particles; i++)
    {
        for (j = i + 1; j < global.N_particles; j++)
        {
            distance_squared_folded_PBC(global.particle_x[i], global.particle_y[i], 
                global.particle_x[j], global.particle_y[j], &dr2, &dx, &dy);
                
            //if (global.time==0)
            //    if ((i==150)||(j==150)) printf("(%d %d %lf)\n", i, j, dr2);
                
            if (dr2 < 36.0)
            {
                global.Verletlisti[global.N_Verlet] = i;
                global.Verletlistj[global.N_Verlet] = j;
                
                //if (global.time==0)
                //if ((i==150)||(j==150)) printf("(%d %d)\n", i, j);
                
                global.N_Verlet++;
                if (global.N_Verlet >= global.N_Verlet_max)
                {
                    printf("Verlet list reallocated from %d\n", global.N_Verlet_max);
                    global.N_Verlet_max = (int)(1.1*global.N_Verlet);
                    global.Verletlisti = (int *) realloc(global.Verletlisti , global.N_Verlet_max * sizeof(int));
                    global.Verletlistj = (int *) realloc(global.Verletlistj , global.N_Verlet_max * sizeof(int));
                    printf("New Verlet list max size = %d\n", global.N_Verlet_max);
                }
            }
        }
    }

    //printf("Counted Verlet list lenght = %d\n", global.N_Verlet);
    //fflush(stdout);

    global.flag_to_rebuild_Verlet = 0;
    for (i = 0; i < global.N_particles; i++)
    {
        global.particle_dx_so_far[i] = 0.0;
        global.particle_dy_so_far[i] = 0.0;
    }
}

void rebuild_pinning_grid()
{
    unsigned int i, j, k;
    unsigned int gi, gj;

    if (global.pinningsite_grid == NULL)
    {
        //build the pinningsite grid for the first time;
        global.Nx_pinningsite_grid = (int) (global.SX / (2 * global.pinningsite_R)) + 1;
        global.Ny_pinningsite_grid = (int) (global.SY / (2 * global.pinningsite_R)) + 1;
        // global.Nx_pinningsite_grid = 3;
        // global.Ny_pinningsite_grid = 3;
        global.max_pinningsite_per_grid = 500;
        printf("Pinning sites grid is %d x %d\n", global.Nx_pinningsite_grid, global.Ny_pinningsite_grid);
        global.pinningsite_grid_dx = global.SX / global.Nx_pinningsite_grid;
        global.pinningsite_grid_dy = global.SY / global.Ny_pinningsite_grid;
        printf("Pinning sites cell is %.2lf x %.2lf\n", global.pinningsite_grid_dx, global.pinningsite_grid_dy);
        
        //initialize the grid
        global.pinningsite_grid = (int ***) malloc(global.Nx_pinningsite_grid * sizeof(int **));
        for(i = 0; i < global.Nx_pinningsite_grid; i++) {
            global.pinningsite_grid[i] = (int **) malloc(global.Ny_pinningsite_grid * sizeof(int *));
            for (j = 0; j < global.Ny_pinningsite_grid; j++)
                global.pinningsite_grid[i][j] = (int *) malloc(global.max_pinningsite_per_grid * sizeof(int)); 
        }
    }
        
    //always do this - zero the values
    for (i = 0; i < global.Nx_pinningsite_grid; i++)
        for (j = 0; j < global.Ny_pinningsite_grid; j++)
            for (k = 0; k < global.max_pinningsite_per_grid; k++)
                global.pinningsite_grid[i][j][k] = -1;
        
    //always do this - fill up the values
    for (i = 0; i < global.N_pinningsites; i++)
    {
        gi = (unsigned int) (global.pinningsite_x[i] / global.pinningsite_grid_dx);
        gj = (unsigned int) (global.pinningsite_y[i] / global.pinningsite_grid_dy);
        
        /*if ((gi>=global.Nx_pinningsite_grid)||(gj>=global.Ny_pinningsite_grid))
            {
            printf("Out of bounds pinningsite\n");
            exit(1);
            }*/
        //this cannot fail (hopefully)
        k = 0;
        while (k < global.max_pinningsite_per_grid && global.pinningsite_grid[gi][gj][k] != -1) k++;
        if (k == global.max_pinningsite_per_grid)
        {
            COLOR_WARNING;
            printf("Too much pinningsite in one grid!\n");
            COLOR_DEFAULT;
        } else {
            global.pinningsite_grid[gi][gj][k] = i;
        }
    }
}


//calculates the shortest distance squared between 2 points in a PBC configuration
//this is squared because I want to save on sqrt with the lookup table
//also used by the Verlet rebuild flag check where I check how much a particle moved
void distance_squared_folded_PBC(double x0, double y0, double x1, double y1,  double *r2_return,  double *dx_return, double *dy_return)
{
    double dr2;
    double dx, dy;

    dx = x1 - x0;
    dy = y1 - y0;

    //PBC fold back
    //if any distance is larger than half the box
    //the copy in the neighboring box is closer
    if (dx > global.halfSX) dx -= global.SX;
    if (dx <= -global.halfSX) dx += global.SX;
    if (dy > global.halfSY) dy -= global.SY;
    if (dy <= -global.halfSY) dy += global.SY;

    dr2 = dx*dx + dy*dy;

    //if ((global.time==0)&&(sqrt(dr2)<6.0))
    //    {
    //    printf("%lf %lf %lf %lf -> %lf %lf = %lf\n", x0, y0, x1, y1, dx, dy,  sqrt(dr2));
    //    }


    *r2_return = dr2;
    *dx_return = dx;
    *dy_return = dy;
}


//moves the particles one time step
void move_particles()
{
    int i;
    double dx, dy;

    for(i = 0; i < global.N_particles; i++)
    {
        dx = global.particle_fx[i] * global.dt;
        dy = global.particle_fy[i] * global.dt;
        
        global.particle_x[i] += dx;
        global.particle_y[i] += dy;
        
        global.particle_dx_so_far[i] += dx;
        global.particle_dy_so_far[i] += dy;

        global.particle_all_dx[i] += dx;
        global.particle_all_dy[i] += dy;
        
        fold_particle_back_PBC(i);
        
        global.particle_fx[i] = 0.0;
        global.particle_fy[i] = 0.0;
    }
}

void fold_particle_back_PBC(int i)
{
    //fold back the particle into the pBC simulation box
    //assumes it did not jump more thana  box length
    //if it did the simulation is already broken anyhow

    if (global.particle_x[i] < 0) global.particle_x[i] += global.SX;
    if (global.particle_y[i] < 0) global.particle_y[i] += global.SY;
    if (global.particle_x[i] >= global.SX) global.particle_x[i] -= global.SX;
    if (global.particle_y[i] >= global.SY) global.particle_y[i] -= global.SY;
}

void write_cmovie_frame()
{
    int i;
    float floatholder;
    int intholder;

    //implement the color testing
    // test_program_by_coloring();

    //legacy cmovie format for del-plot

    // intholder = global.N_pinningsites + global.N_particles;
    // fwrite(&intholder, sizeof(int), 1, global.moviefile);

    // intholder = global.time;
    // fwrite(&intholder, sizeof(int), 1, global.moviefile);

    // for (i = 0; i < global.N_pinningsites; i++)
    // {
    //     intholder = global.pinningsite_color[i];
    //     fwrite(&intholder, sizeof(int), 1, global.moviefile);
    //     intholder = i;//ID
    //     fwrite(&intholder, sizeof(int), 1, global.moviefile);
    //     floatholder = (float)global.pinningsite_x[i];
    //     fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
    //     floatholder = (float)global.pinningsite_y[i];
    //     fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
    //     floatholder = global.pinningsite_R;//cum_disp, cmovie format
    //     fwrite(&floatholder, sizeof(float), 1, global.moviefile);
    // }


    // for (i = 0; i < global.N_particles; i++)
    // {
    //     intholder = global.particle_color[i];
    //     fwrite(&intholder, sizeof(int), 1, global.moviefile);
    //     intholder = i;//ID
    //     fwrite(&intholder, sizeof(int), 1, global.moviefile);
    //     floatholder = (float)global.particle_x[i];
    //     fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
    //     floatholder = (float)global.particle_y[i];
    //     fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
    //     floatholder = global.pinningsite_R / 3.0;//cum_disp, cmovie format
    //     fwrite(&floatholder, sizeof(float), 1, global.moviefile);
    // }

    float dr2, dx, dy;
    float r, x, y;
    int i, j, gi, gj, gj2, n, m, k, count;

    r = global.pinningsite_R;

    x = global.particle_x[50];
    y = global.particle_y[50];

    count = 1;

    gi = (int) (x / global.pinningsite_grid_dx) - 1;
    gj = (int) (y / global.pinningsite_grid_dy) - 1;

    for (m = 0; m < 3; m++)
    {
        if (gi == -1) gi = global.Nx_pinningsite_grid - 1;
        else if (gi == global.Nx_pinningsite_grid - 1) gi = 0;
        else gi ++;

        gj2 = gj;
        for (n = 0; n < 3; n++)
        {
            if (gj2 == -1) gj2 = global.Ny_pinningsite_grid - 1;
            else if (gj2 == global.Ny_pinningsite_grid - 1) gj2 = 0;
            else gj2 ++;

            j = 0;
            while (j < global.max_pinningsite_per_grid && global.pinningsite_grid[gi][gj2][j] != -1)
            {
                k = global.pinningsite_grid[gi][gj2][j];
                count ++;
                distance_squared_folded_PBC(x, y, global.pinningsite_x[k], global.pinningsite_y[k], &dr2, &dx, &dy);
                if (dr2 < r * r)
                {
                    global.particle_fx[i] += dx / r * global.pinningsite_force;
                    global.particle_fy[i] += dy / r * global.pinningsite_force;
                    global.particle_in_pinningsite[k] ++;
                }
                j ++;
            }
        }
    }

    intholder = count;
    fwrite(&intholder, sizeof(int), 1, global.moviefile);

    intholder = global.time;
    fwrite(&intholder, sizeof(int), 1, global.moviefile);

    for (i = 0; i < global.N_pinningsites; i++)
    {
        intholder = global.pinningsite_color[i];
        fwrite(&intholder, sizeof(int), 1, global.moviefile);
        intholder = i;//ID
        fwrite(&intholder, sizeof(int), 1, global.moviefile);
        floatholder = (float)global.pinningsite_x[i];
        fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
        floatholder = (float)global.pinningsite_y[i];
        fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
        floatholder = global.pinningsite_R;//cum_disp, cmovie format
        fwrite(&floatholder, sizeof(float), 1, global.moviefile);
    }


    for (i = 0; i < global.N_particles; i++)
    {
        intholder = global.particle_color[i];
        fwrite(&intholder, sizeof(int), 1, global.moviefile);
        intholder = i;//ID
        fwrite(&intholder, sizeof(int), 1, global.moviefile);
        floatholder = (float)global.particle_x[i];
        fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
        floatholder = (float)global.particle_y[i];
        fwrite(&floatholder, sizeof(float), 1,  global.moviefile);
        floatholder = global.pinningsite_R / 3.0;//cum_disp, cmovie format
        fwrite(&floatholder, sizeof(float), 1, global.moviefile);
    }
}

void calculate_statistics()
{
}

void write_statistics()
{
    int j;

    fprintf(global.statisticsfile, "%d ", global.time);
    double dx, dy, avg_particle_per_pinningsite;
    double avg_particle_per_horizontal_pinningsite, avg_particle_per_left_up_pinningsite, avg_particle_per_left_down_pinningsite;

    dx = dy = 0.0;
    for (j = 0; j < global.N_particles; j++) {
        dx += global.particle_all_dx[j];
        dy += global.particle_all_dy[j];
    }

    avg_particle_per_pinningsite = 0.0;
    avg_particle_per_horizontal_pinningsite = 0.0;
    avg_particle_per_left_up_pinningsite = 0.0;
    avg_particle_per_left_down_pinningsite = 0.0;
    for (j = 0; j < global.N_pinningsites; j++)
    {
        unsigned int par_per_pin = global.particle_in_pinningsite[j];
        avg_particle_per_pinningsite += par_per_pin;
        if (global.pinningsite_direction_x[j] == 1.0 && global.pinningsite_direction_y[j] == 0.0) avg_particle_per_horizontal_pinningsite += par_per_pin;
        else if (global.pinningsite_direction_x[j] == -0.5 && global.pinningsite_direction_y[j] == -sqrt(3) / 2) avg_particle_per_left_down_pinningsite += par_per_pin;
        else if (global.pinningsite_direction_x[j] == -0.5 && global.pinningsite_direction_y[j] == +sqrt(3) / 2) avg_particle_per_left_up_pinningsite += par_per_pin;
        global.particle_in_pinningsite[j] = 0;
    }

    fprintf(global.statisticsfile, "%lf ", dx / global.N_particles);
    fprintf(global.statisticsfile, "%lf ", dy / global.N_particles);
    fprintf(global.statisticsfile, "%lf ", avg_particle_per_horizontal_pinningsite / global.N_pinningsites / global.statistics_time);
    fprintf(global.statisticsfile, "%lf ", avg_particle_per_left_down_pinningsite / global.N_pinningsites / global.statistics_time);
    fprintf(global.statisticsfile, "%lf ", avg_particle_per_left_up_pinningsite / global.N_pinningsites / global.statistics_time);
    fprintf(global.statisticsfile, "%lf ", avg_particle_per_pinningsite / global.N_pinningsites / global.statistics_time);
    fprintf(global.statisticsfile, "\n");
    fflush(global.statisticsfile);
}



//TESTS


void test_program_by_coloring()
{
    int ii;
    int i, j, k;

    //testing the Verlet list by coloring one particle's neightbors one color
    for(i=0;i<global.N_particles;i++)
        global.particle_color[i] = 2;

    for(ii=0;ii<global.N_Verlet;ii++)
        {
        i = global.Verletlisti[ii];
        j = global.Verletlistj[ii];
        if (i==150) global.particle_color[j] = 0;
        if (j==150) global.particle_color[i] = 0;
        }
    //printf("\n");

    //testing the pinning grid by coloring the pins
    // for(i = 0; i < global.Nx_pinningsite_grid; i++)
    //     for(j = 0; j < global.Ny_pinningsite_grid; j++)
    //         for (k = 0; k < global.max_pinningsite_per_grid; k++)
    //             if (global.pinningsite_grid[i][j][k] == -1)
    //                 break;
    //             else
    //             {
    //                 int d = global.pinningsite_grid[i][j][k];
    //                 global.pinningsite_color[d] = (i + j) % 2 + 4;
    //                 // if (d < 24)
    //                 //     printf("%d: (%lf, %lf) -> (%d, %d) -> %d\n", d, global.pinningsite_x[d], global.pinningsite_y[d], i, j, global.pinningsite_color[d]);
    //             }
}

void write_time()
{
	printf("Simulation\n\t- started at:\t%s\t", 
        ctime(&global.start_time));
    printf("- ended at:\t%s\t- total runtime: ", 
        ctime(&global.end_time));
    double time_diff = difftime(global.end_time, global.start_time);
    int day = time_diff / (60 * 60 * 24);
    time_diff -= day * 60 * 60 * 24;
    int hour = time_diff / (60 * 60);
    time_diff -= hour * 60 * 60;
    int minute = time_diff / 60;
    time_diff -= minute * 60;
    int sec = time_diff;
	COLOR_NOTE;
    if (day != 0) 
        if (day == 1) printf("1 day ");
        else printf("%d days ", day);
    if (hour != 0)
        if (hour == 1) printf("1 hour ");
        else printf("%d hours ", hour);
    if (minute != 0) 
        if (minute == 1) printf("1 minute ");
        else printf("%d minutes ", minute);
    printf("%d", sec);
    if (sec <= 1) printf(" second\n");
        else printf(" seconds\n");
	COLOR_DEFAULT;
}

void delete_arrays() 
{
    free(global.moviefile_name);
    free(global.statisticsfile_name);
    free(global.pinningsite_x);
    free(global.pinningsite_y);
    free(global.pinningsite_fx);
    free(global.pinningsite_fy);
    free(global.pinningsite_color);
    free(global.pinningsite_direction_x);
    free(global.pinningsite_direction_y);
    free(global.pinningsite_dx_so_far);
    free(global.pinningsite_dy_so_far);
    free(global.particle_in_pinningsite);
    
    free(global.particle_x);
    free(global.particle_y);
    free(global.particle_fx);
    free(global.particle_fy);
    free(global.particle_color);
    free(global.particle_direction_x);
    free(global.particle_direction_y);
    free(global.particle_dx_so_far);
    free(global.particle_dy_so_far);
    free(global.particle_all_dx);
    free(global.particle_all_dy);

    free(global.Verletlisti);
    free(global.Verletlistj);

    unsigned int i, j;
    for (i = 0; i < global.Nx_pinningsite_grid; i++) {
        for (j = 0; j < global.Ny_pinningsite_grid; j++)
            free(global.pinningsite_grid[i][j]);
        free(global.pinningsite_grid[i]);
    }
    free(global.pinningsite_grid);
}
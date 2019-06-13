//
//  running.c
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#include "running.h"
#include "globaldata.h"
#include "utils.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.1415926535

// runs the simulation for the required steps
void run_simulation()
{
	int console_width;
	// get console width and substracting start and end character from progressbar and leaving space for writing the precentage 
	console_width = get_console_width() - 11;
	rebuild_Verlet_list();		// build Verlet list for the first time
	print_log(stdout, NOTE, __FILE__, __LINE__, "Simulation started", 0);

	// save time when the simulation started
	time(&global.start_time);
	for (global.time = 0; global.time <= global.total_time; global.time++)
	{
		// adjust_pinningsite_directions();
		// rotate_pinningsite_directions();
		calculate_external_forces_on_pinningsites();
		move_pinningsites();
		rebuild_pinning_grid();

		// calculate_external_forces_on_particles();
		calculate_pairwise_forces();
		calculate_pinning_force_on_particles();
		// calclulate_pinning_force_on_particles_without_grid();

		// this is where we calculate statistics
		calculate_statistics();
		move_particles();

		check_Verlet_rebuild_condition_and_set_flag();

		// one particle moved enough to rebuild the Verlet list
		if (global.flag_to_rebuild_Verlet == 1)
			rebuild_Verlet_list();

		// echo time
		if (global.time % global.echo_time == 0)
		{
			console_width = get_console_width() - 12;
			progress_bar((double) global.time / (double) global.total_time,
						 console_width, "%d / %d", global.time,
						 global.total_time);
		}

		// statistics time
		if (global.time % global.statistics_time == 0 && global.time != 0)
		 write_statistics();

		// movie write time
		if (global.time % global.movie_time == 0)
		 write_cmovie_frame();
	}

	delete_arrays();
	// fclose(global.moviefile);
	// fclose(global.statisticsfile);

	// saving time when the simulation finished
	time(&global.end_time);
	// printing out the execution time of the simulation
	write_time();
}

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

	for (i = 0; i < global.N_particles; i++)
	{
		newx =
			global.particle_direction_x[i] * cost -
			global.particle_direction_y[i] * sint;
		newy =
			global.particle_direction_x[i] * sint +
			global.particle_direction_y[i] * cost;

		global.particle_direction_x[i] = newx;
		global.particle_direction_y[i] = newy;
	}

	printf("%lf %lf\n", global.particle_direction_x[i],
		   global.particle_direction_y[i]);

}

void adjust_pinningsite_directions()
{
	int i, k;

	if ((global.pinningsite_dx_so_far[0] * global.pinningsite_dx_so_far[0] +
		 global.pinningsite_dy_so_far[0] * global.pinningsite_dy_so_far[0])
		>=global.pinning_lattice_constant * global.pinning_lattice_constant *
		1.0)
	{
		//change directions now
		for (k = 0; k < global.N_pinningsites; k++)
			switch ((global.pinningsite_color[k] - 4 +
					 global.pinning_direction_change) %3)
			{
				case 0:
				{
					global.pinningsite_direction_x[k] = -0.5;
					global.pinningsite_direction_y[k] = -sqrt(3) / 2.0;
					break;
				}
				case 1:
				{
					global.pinningsite_direction_x[k] = 1.0;
					global.pinningsite_direction_y[k] = 0.0;
					break;
				}
				case 2:
				{
					global.pinningsite_direction_x[k] = -0.5;
					global.pinningsite_direction_y[k] = +sqrt(3) / 2.0;
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

	for (i = 0; i < global.N_pinningsites; i++)
	{
		// for every pinningsite calculate the external force on it by multiplying its magnitude with its direction for each (x and y) directions
		global.pinningsite_fx[i] +=
			global.pinningsite_direction_x[i] * global.pinning_driving_force;
		global.pinningsite_fy[i] +=
			global.pinningsite_direction_y[i] * global.pinning_driving_force;
	}
}

void move_pinningsites()
{
	int i;
	double dx, dy;

	for (i = 0; i < global.N_pinningsites; i++)
	{
		// calculating travelled distence in a time step depending on the force affecting it
		dx = global.pinningsite_fx[i] * global.dt;
		dy = global.pinningsite_fy[i] * global.dt;

		// addig the calculated distaces to the current position
		global.pinningsite_x[i] += dx;
		global.pinningsite_y[i] += dy;

		// and add also to the total travelling so far
		global.pinningsite_dx_so_far[i] += dx;
		global.pinningsite_dy_so_far[i] += dy;

		fold_back_PBC(&global.pinningsite_x[i], &global.pinningsite_y[i]);

		// set forces to zero, prepare to the next iteration
		global.pinningsite_fx[i] = 0.0;
		global.pinningsite_fy[i] = 0.0;
	}
}

void fold_back_PBC(double *x, double *y)
{
	if (*x < 0)
		*x += global.SX;
	if (*y < 0)
		*y += global.SY;
	if (*x >= global.SX)
		*x -= global.SX;
	if (*y >= global.SY)
		*y -= global.SY;
}

void calculate_external_forces_on_particles()
{
	int i;

	for (i = 0; i < global.N_particles; i++)
	{
		// for each particle add external for to total force by multiplying the force's magnitude with its direction
		global.particle_fx[i] +=
			global.particle_direction_x[i] * global.particle_driving_force;
		global.particle_fy[i] +=
			global.particle_direction_y[i] * global.particle_driving_force;
	}
}

void calculate_pairwise_forces()
{
	int ii, i, j;
	double r, r2, f;
	double dx, dy;

	// iterating through the elements in the Verlet list
	for (ii = 0; ii < global.N_Verlet; ii++)
	{
		// obtain the i, j particle IDs from the Verlet list
		i = global.Verletlisti[ii];
		j = global.Verletlistj[ii];

		// calculating distance between them
		distance_squared_folded_PBC(global.particle_x[i],
									global.particle_y[i],
									global.particle_x[j],
									global.particle_y[j], &r2, &dx, &dy);

		// non-tabulated version
		// try to not divide just multiply division is costly
		r = sqrt(r2);
		if (r < 0.2)
		{
			// particles are too close so there is a big force between them
			f = 100.0;
		}
		else
		{
			// calculating force depending on the distance between the particles
			f = 1 / r2 * exp(-r *
							 global.partile_particle_screening_wavevector);
		}

		// division and multiplication for projection to the x, y axes
		f = f / r;

		// adding and substracting the force to the total forces
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
	double dr2, dx, dy;
	double r, x, y;

	r = global.pinningsite_R;
	for (i = 0; i < global.N_particles; i++)
	{
		// get current particle position
		x = global.particle_x[i];
		y = global.particle_y[i];

		// we have to check for each particle the surrounding (9) squares in the grid, because in those could
		// be pinnningsites close enough to effect it

		// we subtract the particle's cell position by one, to start with the left upper square the
		// iteration. Because this subtraction we can get -1 if the particle is in the first row or column
		// therefore gi and gj have to be int.
		gi = (int) (x / global.pinningsite_grid_dx) -1;
		gj = (int) (y / global.pinningsite_grid_dy) -1;

		// iterate through rows in pinningsite grid
		for (m = 0; m < 3; m++)
		{
			// fold pinning grid row back
			if (gi == -1)
				gi = global.Nx_pinningsite_grid - 1;
			else if (gi == global.Nx_pinningsite_grid - 1)
				gi = 0;
			else if (m != 0)
				gi++;

			// remember original column index
			gj2 = gj;
			// iterate through columns in pinningsite grid
			for (n = 0; n < 3; n++)
			{
				// fold pinning grid column back
				if (gj2 == -1)
					gj2 = global.Ny_pinningsite_grid - 1;
				else if (gj2 == global.Ny_pinningsite_grid - 1)
					gj2 = 0;
				else if (n != 0)
					gj2++;

				j = 0;
				// iterate through pinningsites in the current grid cell
				while (j < global.max_pinningsite_per_grid
					   && global.pinningsite_grid[gi][gj2][j] != -1)
				{
					// getting pinningsite ID
					k = global.pinningsite_grid[gi][gj2][j];
					// calculate distance between pinningsite center and particle
					distance_squared_folded_PBC(x, y, global.pinningsite_x[k],
												global.pinningsite_y[k], &dr2,
												&dx, &dy);
					if (dr2 < r * r)
					{
						// calculate force depending on their distance
						global.particle_fx[i] +=
							dx / r * global.pinningsite_force;
						global.particle_fy[i] +=
							dy / r * global.pinningsite_force;
						// calculate number of particles in the current pinningsite
						global.particle_in_pinningsite[k]++;
					}
					j++;
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
	// iterate through each pinningsite
	for (j = 0; j < global.N_pinningsites; j++)
	{
		// get pinningsite's position
		gi = global.pinningsite_x[j];
		gj = global.pinningsite_y[j];

		// iterate through each particle
		for (i = 0; i < global.N_particles; i++)
		{
			// get particle's position
			x = global.particle_x[i];
			y = global.particle_y[i];

			// calculate distance between pinningsite and particle
			distance_squared_folded_PBC(gi, gj, x, y, &dr2, &dx, &dy);
			if (dr2 < r * r)
			{
				// caluclate force depending on their distance
				global.particle_fx[i] -= dx / r * global.pinningsite_force;
				global.particle_fy[i] -= dy / r * global.pinningsite_force;
				// calculate number of particles in the current pinningsite
				global.particle_in_pinningsite[j]++;
			}
		}
	}
}

void check_Verlet_rebuild_condition_and_set_flag()
{
	int i;
	double dr2;

	global.flag_to_rebuild_Verlet = 0;

	// iterate through partilces
	for (i = 0; i < global.N_particles; i++)
	{
		// calculate total travelled distance square from the position when it was placed to the Verlet list
		dr2 = global.particle_dx_so_far[i] * global.particle_dx_so_far[i] +
			global.particle_dy_so_far[i] * global.particle_dy_so_far[i];;

		if (dr2 >= global.Verlet_intershell_squared)
		{
			// just one match is enough to rebuild the lists
			global.flag_to_rebuild_Verlet = 1;
			break;				// exit the for cycle
		}
	}
}

void rebuild_Verlet_list()
{
	int i, j;
	double dr2, dx, dy;
	double density, estimation;

	if (global.N_Verlet_max == 0)
	{
		// initialize the Verlet list for the first time in the simulation
		// calculate density of the system
		density =
			global.N_particles / (double) global.SX / (double) global.SY;

		estimation = density *
			PI * global.Verlet_cutoff_distance *
			global.Verlet_cutoff_distance;

		// estimating number of particle pairs in Verlet list
		global.N_Verlet_max = (int) estimation * global.N_particles;

		global.Verletlisti =
			(int *) malloc(global.N_Verlet_max * sizeof(int));
		global.Verletlistj =
			(int *) malloc(global.N_Verlet_max * sizeof(int));

		print_log(stdout, NOTE, __FILE__, __LINE__,
				  "Build Verlet list for the first time", 3,
				  "System density is %.3lf\n", density,
				  "Particles in a R = %.2lf shell = %lf\n",
				  global.Verlet_cutoff_distance, estimation,
				  "Estimated N_Verlet_max = %d\n", global.N_Verlet_max);
	}

	// build the Verlet list - filling it with particle pairs
	global.N_Verlet = 0;

	for (i = 0; i < global.N_particles; i++)
	{
		for (j = i + 1; j < global.N_particles; j++)
		{
			// calculating distance between two particle
			distance_squared_folded_PBC(global.particle_x[i],
										global.particle_y[i],
										global.particle_x[j],
										global.particle_y[j], &dr2, &dx, &dy);

			// if the two particle is close enough to affect each other place them in the Verlet list
			if (dr2 < 36.0)
			{
				global.Verletlisti[global.N_Verlet] = i;
				global.Verletlistj[global.N_Verlet] = j;

				global.N_Verlet++;
				// if end of list is reached, then reallocate arrays
				if (global.N_Verlet >= global.N_Verlet_max)
				{
					int old = global.N_Verlet_max;
					// calculate new max size of Verlet list
					global.N_Verlet_max = global.N_Verlet + 10;
					// reallocate arrays
					global.Verletlisti = (int *) realloc(global.Verletlisti,
														 global.N_Verlet_max *
														 sizeof(int));
					global.Verletlistj = (int *) realloc(global.Verletlistj,
														 global.N_Verlet_max *
														 sizeof(int));
					print_log(stdout, WARNING, __FILE__, __LINE__,
							  "Reallocate Verlet list", 2,
							  "Old max size = %d\n", old,
							  "New max size = %d\n", global.N_Verlet_max);
				}
			}
		}
	}

	// Verlet list was rebiult => setting rebuild flag to False and travelled distance since rebuilt to 0
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
		// build the pinningsite grid for the first time;
		// calculate pinningsite grid size depending on the diameter of a pinningsite
		global.Nx_pinningsite_grid =
			(int) (global.SX / (2 * global.pinningsite_R)) +1;
		global.Ny_pinningsite_grid =
			(int) (global.SY / (2 * global.pinningsite_R)) +1;
		// set maximum number of pinningsites allowed to be in a grid at the same time
		global.max_pinningsite_per_grid = 5;
		// calculate cell width and height in the system's measure
		global.pinningsite_grid_dx = global.SX / global.Nx_pinningsite_grid;
		global.pinningsite_grid_dy = global.SY / global.Ny_pinningsite_grid;
		print_log(stdout, NOTE, __FILE__, __LINE__, "Build pinningsite grid",
				  2, "Pinning sites grid is %d x %d\n",
				  global.Nx_pinningsite_grid, global.Ny_pinningsite_grid,
				  "Pinning sites cell is %.2lf x %.2lf\n",
				  global.pinningsite_grid_dx, global.pinningsite_grid_dy);

		// allocate the grid
		global.pinningsite_grid =
			(int ***) malloc(global.Nx_pinningsite_grid * sizeof(int **));
		for (i = 0; i < global.Nx_pinningsite_grid; i++)
		{
			// allocate columns
			global.pinningsite_grid[i] =
				(int **) malloc(global.Ny_pinningsite_grid * sizeof(int *));
			for (j = 0; j < global.Ny_pinningsite_grid; j++)
			{
				// allocate elements in a cell
				global.pinningsite_grid[i][j] =
					(int *) malloc(global.max_pinningsite_per_grid *
								   sizeof(int));
			}
		}
	}

	// always do this - zero the values
	for (i = 0; i < global.Nx_pinningsite_grid; i++)
		for (j = 0; j < global.Ny_pinningsite_grid; j++)
			for (k = 0; k < global.max_pinningsite_per_grid; k++)
				global.pinningsite_grid[i][j][k] = -1;

	// always do this - fill up the values
	for (i = 0; i < global.N_pinningsites; i++)
	{
		// calculate cell's row and column from pinningsite's position
		gi = (unsigned int) (global.pinningsite_x[i] /
							 global.pinningsite_grid_dx);
		gj = (unsigned int) (global.pinningsite_y[i] /
							 global.pinningsite_grid_dy);
		k = 0;
		// find empty (-1) spot in pinningsite list in the cell
		while (k < global.max_pinningsite_per_grid
			   && global.pinningsite_grid[gi][gj][k] != -1)
			k++;
		// checking pinningsit number in the cell
		if (k == global.max_pinningsite_per_grid)
			print_log(stdout, WARNING, __FILE__, __LINE__,
					  "Too much pinningsite in one grid!", 0);
		else
			global.pinningsite_grid[gi][gj][k] = i;
	}
}


void distance_squared_folded_PBC(double x0, double y0, double x1, double y1,
								 double *r2_return, double *dx_return,
								 double *dy_return)
{
	double dr2;
	double dx, dy;

	dx = x1 - x0;
	dy = y1 - y0;

	// PBC fold back
	// if any distance is larger than half the box
	// the copy in the neighboring box is closer
	if (dx > global.halfSX)
		dx -= global.SX;
	if (dx <= -global.halfSX)
		dx += global.SX;
	if (dy > global.halfSY)
		dy -= global.SY;
	if (dy <= -global.halfSY)
		dy += global.SY;

	dr2 = dx * dx + dy * dy;

	*r2_return = dr2;
	*dx_return = dx;
	*dy_return = dy;
}

void move_particles()
{
	int i;
	double dx, dy;

	for (i = 0; i < global.N_particles; i++)
	{
		// calculate travelled distance in the time step
		dx = global.particle_fx[i] * global.dt;
		dy = global.particle_fy[i] * global.dt;

		// move particle
		global.particle_x[i] += dx;
		global.particle_y[i] += dy;

		// increase total travelled distance since last rebuilt of Verlet list
		global.particle_dx_so_far[i] += dx;
		global.particle_dy_so_far[i] += dy;

		// increase total travelled distance
		global.particle_all_dx[i] += fabs(dx);
		global.particle_all_dy[i] += fabs(dy);
		global.particle_all_dr2[i] += dx * dx + dy * dy;

		fold_back_PBC(&global.particle_x[i], &global.particle_y[i]);

		// reseting force of particle to 0
		global.particle_fx[i] = 0.0;
		global.particle_fy[i] = 0.0;
	}
}

void write_cmovie_frame()
{
	int i;
	float floatholder;
	int intholder;

	// implement the color testing
	test_program_by_coloring();

	// legacy cmovie format for del-plot

	// write total number of object on the frame
	intholder = global.N_pinningsites + global.N_particles;
	fwrite(&intholder, sizeof(int), 1, global.moviefile);

	// write time step
	intholder = global.time;
	fwrite(&intholder, sizeof(int), 1, global.moviefile);

	// write pinningsites' data
	for (i = 0; i < global.N_pinningsites; i++)
	{
		// color
		intholder = global.pinningsite_color[i];
		fwrite(&intholder, sizeof(int), 1, global.moviefile);
		// ID
		intholder = i;
		fwrite(&intholder, sizeof(int), 1, global.moviefile);
		// x coordinate
		floatholder = (float) global.pinningsite_x[i];
		fwrite(&floatholder, sizeof(float), 1, global.moviefile);
		// y coordinate
		floatholder = (float) global.pinningsite_y[i];
		fwrite(&floatholder, sizeof(float), 1, global.moviefile);
		// cum_disp, cmovie format
		floatholder = global.pinningsite_R;
		fwrite(&floatholder, sizeof(float), 1, global.moviefile);
	}

	// write particles' data
	for (i = 0; i < global.N_particles; i++)
	{
		// color
		intholder = global.particle_color[i];
		fwrite(&intholder, sizeof(int), 1, global.moviefile);
		// ID
		intholder = i;
		fwrite(&intholder, sizeof(int), 1, global.moviefile);
		// x coordinate
		floatholder = (float) global.particle_x[i];
		fwrite(&floatholder, sizeof(float), 1, global.moviefile);
		// y coordinate
		floatholder = (float) global.particle_y[i];
		fwrite(&floatholder, sizeof(float), 1, global.moviefile);
		//cum_disp, cmovie format
		floatholder = global.pinningsite_R / 3.0;
		fwrite(&floatholder, sizeof(float), 1, global.moviefile);
	}
}

void calculate_statistics()
{
}

void write_statistics()
{
	int j;

	// wrute time step
	fprintf(global.statisticsfile, "%d ", global.time);
	double dx, dy, avg_particle_per_pinningsite, dr2;
	double avg_particle_per_horizontal_pinningsite,
		avg_particle_per_left_up_pinningsite,
		avg_particle_per_left_down_pinningsite;

	dr2 = dx = dy = 0.0;
	// calculate travelled distance average
	for (j = 0; j < global.N_particles; j++)
	{
		dx += global.particle_all_dx[j];
		dy += global.particle_all_dy[j];
		dr2 += global.particle_all_dr2[j];
	}

	// calculate particle number per pinningsite in total and for the 3 type of pinningsite
	avg_particle_per_pinningsite = 0.0;
	avg_particle_per_horizontal_pinningsite = 0.0;
	avg_particle_per_left_up_pinningsite = 0.0;
	avg_particle_per_left_down_pinningsite = 0.0;
	for (j = 0; j < global.N_pinningsites; j++)
	{
		unsigned int par_per_pin = global.particle_in_pinningsite[j];
		avg_particle_per_pinningsite += par_per_pin;
		if (global.pinningsite_direction_x[j] == 1.0
			&& global.pinningsite_direction_y[j] == 0.0)
			avg_particle_per_horizontal_pinningsite += par_per_pin;
		else if (global.pinningsite_direction_x[j] == -0.5
				 && global.pinningsite_direction_y[j] == -sqrt(3) / 2)
			avg_particle_per_left_down_pinningsite += par_per_pin;
		else if (global.pinningsite_direction_x[j] == -0.5
				 && global.pinningsite_direction_y[j] == +sqrt(3) / 2)
			avg_particle_per_left_up_pinningsite += par_per_pin;
		global.particle_in_pinningsite[j] = 0;
	}

	// wirte calculated values into statistics file
	fprintf(global.statisticsfile, "%lf ", dx / global.N_particles);
	fprintf(global.statisticsfile, "%lf ", dy / global.N_particles);
	// fprintf(global.statisticsfile, "%lf ", dr2 / global.N_particles);
	fprintf(global.statisticsfile, "%lf ",
			avg_particle_per_horizontal_pinningsite / global.N_pinningsites /
			global.statistics_time);
	fprintf(global.statisticsfile, "%lf ",
			avg_particle_per_left_down_pinningsite / global.N_pinningsites /
			global.statistics_time);
	fprintf(global.statisticsfile, "%lf ",
			avg_particle_per_left_up_pinningsite / global.N_pinningsites /
			global.statistics_time);
	fprintf(global.statisticsfile, "%lf ",
			avg_particle_per_pinningsite / global.N_pinningsites /
			global.statistics_time);
	fprintf(global.statisticsfile, "\n");
	fflush(global.statisticsfile);
}



//TESTS


void test_program_by_coloring()
{
	int ii;
	int i, j, k;

	//testing the Verlet list by coloring one particle's neightbors one color
	for (i = 0; i < global.N_particles; i++)
		global.particle_color[i] = 2;

	global.particle_color[150] = 3;

	for (ii = 0; ii < global.N_Verlet; ii++)
	{
		i = global.Verletlisti[ii];
		j = global.Verletlistj[ii];
		if (i == 150)
			global.particle_color[j] = 0;
		if (j == 150)
			global.particle_color[i] = 0;
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
	// print starting and ending time
	printf("Simulation\n\t- started at:\t%s\t", ctime(&global.start_time));
	printf("- ended at:\t%s\t- total runtime: ", ctime(&global.end_time));
	// calculate difference between end and start.
	// it will result in milliseconds
	double time_diff = difftime(global.end_time, global.start_time);
	// calculate number of dasys
	int day = time_diff / (60 * 60 * 24);
	time_diff -= day * 60 * 60 * 24;
	// calculate number of hours
	int hour = time_diff / (60 * 60);
	time_diff -= hour * 60 * 60;
	// calculate minutes
	int minute = time_diff / 60;
	time_diff -= minute * 60;
	// remaining will be the seconds
	int sec = time_diff;
	COLOR_NOTE;
	if (day != 0)
		if (day == 1)
			printf("1 day ");
		else
			printf("%d days ", day);
	if (hour != 0)
		if (hour == 1)
			printf("1 hour ");
		else
			printf("%d hours ", hour);
	if (minute != 0)
		if (minute == 1)
			printf("1 minute ");
		else
			printf("%d minutes ", minute);
	printf("%d", sec);
	if (sec <= 1)
		printf(" second\n");
	else
		printf(" seconds\n");
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
	for (i = 0; i < global.Nx_pinningsite_grid; i++)
	{
		for (j = 0; j < global.Ny_pinningsite_grid; j++)
			free(global.pinningsite_grid[i][j]);
		free(global.pinningsite_grid[i]);
	}
	free(global.pinningsite_grid);
}

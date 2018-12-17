//
//  initializer.c
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#include <math.h>
#include <regex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "initializer.h"
#include "globaldata.h"

//reads the parameter file
//adjusts the parameters accordingly
int read_init_file(const char* filename)
{
	FILE *f = fopen(filename, "rt");
	char tmp[255];
    int n;
    size_t len_tmp, len_json;

    n = 3;
    global.JSON = (char*)malloc(n);
    len_json = 1;
	if (f != NULL)
	{  
        printf("Initializing data from the parameterfile\n");
        printf("Opened parameter file %s\n", filename);
		while (!feof(f))
		{
			fgets(tmp, 255, f);
            len_tmp = strlen(tmp);
			len_json += len_tmp;
            if (len_json > n)
            {
                n += 255;
                global.JSON = (char*)realloc(global.JSON, n);
            }
            strncat(global.JSON, tmp, len_tmp);
            global.JSON[len_json] = '\0';
		}
        
		if (global.JSON[strlen(global.JSON) - 6] == '}' && global.JSON[strlen(global.JSON) - 3] == '}')
			global.JSON[strlen(global.JSON) - 5] = '\0';
	} else {
		global.JSON = "{}";
        printf("\033[1;31mFailed to open parameter file %s\033[0m\n", filename);
        fflush(stdout);
        exit(1);
	}

	fclose(f);
    return n;
}

int check_file_name(const char* regex_str, char* toTest)
{
	regex_t regexCompiled;
	int result;
	regmatch_t pmatch[1];

	if (regcomp(&regexCompiled, regex_str, REG_EXTENDED))
	{
		printf("\033[1;31mCould not create regex!\033[0m\n");
		return 0;
	}

	result = !regexec(&regexCompiled, toTest, 1, pmatch, 0) ? ((pmatch[0].rm_so == 0) && (pmatch[0].rm_eo == strlen(toTest))) : 0;
	regfree(&regexCompiled);
	return result;
}

int get_token(int start, int end, const char* token)
{
    int i;
	int l2 = strlen(token);
	for (i = start; i < end; i++)
	{
		if (global.t[i].type == JSMN_STRING)
		{
			int l1 = global.t[i].end - global.t[i].start;
			if (l1 == l2)
			{
				char* tmp = (char*) malloc(l1);
                memcpy(tmp, &global.JSON[global.t[i].start], l1);
				if (strncmp(tmp, token, l1) == 0)
				{
					free(tmp);
					return i;
				}
				free(tmp);
			}
		}
	}
	return -1;
}

void init_data()
{
	//T es decimated_number kimaradt
	jsmn_parser p;
	int n, r, x, y, object_end;
	size_t len, regex_len, s, len_file_path;
	char* regex_str, *toTest;

	regex_len = 36;
    len_file_path = 50;

	len = (size_t) strlen(global.JSON);
	jsmn_init(&p);
	n = jsmn_parse(&p, global.JSON, len, global.t, 0);
	if (n < 0) {
		printf("\033[1;31mSome error occured!\033[0m\n");
		return;
	} else {
		global.t = (jsmntok_t*) malloc(sizeof(jsmntok_t) * n);
	}
	jsmn_init(&p);
	r = jsmn_parse(&p, global.JSON, len, global.t, n);
	if (r < 1 && !global.t[0].type == JSMN_OBJECT) {
		printf("\033[1;31mJSON object not found!\033[0m\n");
		return;
	}
	
	if ((x = get_token(1, r, "simulation_box")) >= 0) {
		if (global.t[x + 1].type == JSMN_OBJECT) {
			object_end = x + 1 + 2*global.t[x + 1].size;
			if ((y = get_token(x + 2, object_end, "SX")) >= 0)
				global.SX = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.SX = 72.0;

			if ((y = get_token(x + 2, object_end, "SY")) >= 0)
				global.SY = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.SY = 10;
		}
	} else {
		global.SX = 72.0;
		global.SY = 72.0;
	}
    printf("Initializing the simulation box\n");
	printf("\tSX: %f\n", global.SX);
	printf("\tSY: %f\n", global.SY);
    global.halfSX = global.SX / 2.0;
    global.halfSY = global.SY / 2.0;

	if ((x = get_token(1, r, "pinningsite")) >= 0) {
		if (global.t[x + 1].type == JSMN_OBJECT) {
			object_end = x + 1 + 2*global.t[x + 1].size;
			if ((y = get_token(x + 2, object_end, "pinningsite_R")) >= 0)
				global.pinningsite_R = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.pinningsite_R = 0.2;

			if ((y = get_token(x + 2, object_end, "pinningsite_force")) >= 0)
				global.pinningsite_force = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.pinningsite_force = 0.6;

			if ((y = get_token(x + 2, object_end, "pinning_driving_force")) >= 0)
				global.pinning_driving_force = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.pinning_driving_force = 0.5;

			if ((y = get_token(x + 2, object_end, "N_pinningsites")) >= 0)
				global.N_pinningsites = atoi(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.N_pinningsites = 567;

            if ((y = get_token(x + 2, object_end, "pinning_direction_change")) >= 0)
				global.pinning_direction_change = atoi(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.pinning_direction_change = 1;
		}
	} else {
		global.pinningsite_R = 0.2;
		global.pinningsite_force = 0.6;
		global.pinning_driving_force = 0.5;
		global.N_pinningsites = 567;
        global.pinning_direction_change = 1;
	}
    printf("Pinningsite:\n");
    printf("\tNumber: %d\n", global.N_pinningsites);
	printf("\tRadius: %f\n", global.pinningsite_R);
    printf("\tForce: %f\n", global.pinningsite_force);
    printf("\tDriving force: %f\n", global.pinning_driving_force);
    printf("\tDirection change: %d\n", global.pinning_direction_change);

    if ((x = get_token(1, r, "particle")) >= 0) {
		if (global.t[x + 1].type == JSMN_OBJECT) {
			object_end = x + 1 + 2*global.t[x + 1].size;
            if ((y = get_token(x + 2, object_end, "N_particles")) >= 0)
				global.N_particles = atoi(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.N_particles = 1000;

			if ((y = get_token(x + 2, object_end, "particle_driving_force")) >= 0)
				global.particle_driving_force = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.particle_driving_force = 0.2;

			if ((y = get_token(x + 2, object_end, "partile_particle_screening_length")) >= 0)
				global.partile_particle_screening_length = atof(substr(global.JSON, global.t[y + 1].start, global.t[y + 1].end - global.t[y + 1].start));
			else
				global.partile_particle_screening_length = 0.6;
		}
	} else {
		global.particle_driving_force = 0.2;
		global.partile_particle_screening_length = 0.6;
	}

    global.partile_particle_screening_wavevector =  1.0 / global.partile_particle_screening_length;

    printf("Particles:\n");
    printf("\tNumber: %d\n", global.N_particles);
    printf("\tDriving force: %f\n", global.particle_driving_force);
    printf("\tParticle particle screening length: %f\n", global.partile_particle_screening_length);
    printf("\tParticle particle screening wavevector: %lf\n", global.partile_particle_screening_wavevector);

    if ((x = get_token(1, r, "dt")) >= 0) 
        global.dt = atof(substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start));
    else
        global.dt = 0.001;

    if ((x = get_token(1, r, "statistics_time")) >= 0) 
        global.statistics_time = atoi(substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start));
    else
        global.statistics_time = 1000;

    if ((x = get_token(1, r, "movie_time")) >= 0) 
        global.movie_time = atoi(substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start));
    else
        global.movie_time = 100;

    if ((x = get_token(1, r, "total_time")) >= 0) 
        global.total_time = atoi(substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start));
    else
        global.total_time = 100000;
    
    if ((x = get_token(1, r, "echo_time")) >= 0) 
        global.echo_time = atoi(substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start));
    else
        global.echo_time = 1000;

    printf("Timestep: %f\n", global.dt);
    printf("Echo times:\n");
    printf("\tTo screen: %d\n", global.echo_time);
    printf("\tTo statistics file: %d\n", global.statistics_time);
    printf("\tTo movie file: %d\n", global.movie_time);
    printf("Total running time: %d\n", global.total_time);

    global.moviefile_name = (char *)malloc(len_file_path);
    strncpy(global.moviefile_name, "./results/movies/",  18);
    global.moviefile_name[18] = '\0';
	if ((x = get_token(1, r, "moviefile")) >= 0) {
		regex_str = (char*) malloc(regex_len);
		strncpy(regex_str, "^[a-z0-9A-Z_]+$", (size_t) 15);
        regex_str[15] = '\0';
        toTest = substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start);
		if (check_file_name(regex_str, toTest)) 
            strncat(toTest, ".mvi", 4);
        s = strlen(regex_str) - 1;
		regex_str[s] = '\0';
		strncat(regex_str, "(\\.mvi)$", (size_t) 7);
        regex_str[s + 7] = '\0';
        len = strlen(global.moviefile_name);
		if (!check_file_name(regex_str, toTest)) {
            strncat(global.moviefile_name, "result.mvi", 10);
            len += 10;
        } else {
            strncat(global.moviefile_name, toTest, strlen(toTest));
            len +=strlen(toTest);
        }
        global.moviefile_name[len] = '\0'; 
		free(regex_str);
        free(toTest);
	} else {
        len = snprintf(NULL, 0, "./results/movies/result%2.2f.mvi", global.pinningsite_force);
        global.moviefile_name = (char *)realloc(global.moviefile_name, len + 1);
		snprintf(global.moviefile_name, len + 1, "./results/movies/result%2.2f.mvi", global.pinningsite_force);
    }
    printf("File names:\n");
	printf("\tMoviefile: %s\n", global.moviefile_name);
    
    global.statisticsfile_name = (char *)malloc(len_file_path);
    strncpy(global.statisticsfile_name, "./results/stats/",  18);
    global.statisticsfile_name[18] = '\0';
	if ((x = get_token(1, r, "statfile")) >= 0) {
		regex_str = (char*) malloc(regex_len);
		strncpy(regex_str, "^[a-z0-9A-Z_]+$", (size_t) 15);
        regex_str[15] = '\0';
		toTest = substr(global.JSON, global.t[x + 1].start, global.t[x + 1].end - global.t[x + 1].start);
		if (check_file_name(regex_str, toTest)) 
            strncat(toTest, ".txt", 4);
        s = strlen(regex_str) - 1;
		regex_str[s] = '\0';
		strncat(regex_str, "(\\.txt)$", (size_t) 7);
        regex_str[s + 7] = '\0';
        len = strlen(global.statisticsfile_name);
		if (!check_file_name(regex_str, toTest)) {
            strncat(global.statisticsfile_name, "stat.txt", 8);
            len += 9;
        } else {
            strncat(global.statisticsfile_name, toTest, strlen(toTest));
            len += strlen(toTest);
        }
        global.statisticsfile_name[len] = '\0';
        free(regex_str);
        free(toTest);
	} else {
        len = snprintf(NULL, 0, "./results/stats/stat%2.2f.txt", global.pinningsite_force);
        global.statisticsfile_name = (char *)realloc(global.statisticsfile_name, len + 1);
		snprintf(global.statisticsfile_name, len + 1, "./results/stats/stat%2.2f.txt", global.pinningsite_force);
    }
	printf("\tStatfile: %s\n", global.statisticsfile_name);

	free(global.t);
}

//init the common variables in the simulation
void init_simulation()
{
    printf("Post Initializing data\n");
    
    global.Verlet_cutoff_distance = 1.5 * global.partile_particle_screening_length;
    global.Verlet_cutoff_distance_squared = global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
    global.Verlet_intershell_squared = global.Verlet_cutoff_distance - global.partile_particle_screening_length;
    global.Verlet_intershell_squared = global.Verlet_intershell_squared / 2.0;
    global.Verlet_intershell_squared *= global.Verlet_intershell_squared;

    printf("Verlet cutoff distance = %.2lf\n", global.Verlet_cutoff_distance);
    printf("Verlet cutoff distance squared = %.2lf\n", global.Verlet_cutoff_distance_squared);
    printf("Half of Verlet intershell distance = %.2lf\n", sqrt(global.Verlet_intershell_squared));
    printf("Half of Verlet intershell distance squared = %.2lf\n", global.Verlet_intershell_squared);

    //zero everything so rebuild Verlet can find this the first time
    global.Verletlisti = NULL;
    global.Verletlistj = NULL;
    global.N_Verlet = 0;
    global.N_Verlet_max = 0;

    //zero everythiong so rebuild_pinning_grid can find this the first time
    global.pinningsite_grid = NULL;
    global.Nx_pinningsite_grid = 0;
    global.Ny_pinningsite_grid = 0;

    //init the unit vectors for the 3 directions
    //left and down
    global.ex[0] =  -0.5;
    global.ey[0] =  -sqrt(3)/2.0;
    //right
    global.ex[1] =  1.0;
    global.ey[1] =  0.0;
    //left and up
    global.ex[2] =  -0.5;
    global.ey[2] =  sqrt(3)/2.0;
        
    global.sumforce[0] = 0.0;
    global.sumforce[1] = 0.0;
    global.sumforce[2] = 0.0;
    
    global.N_sum[0] = 0;
    global.N_sum[1] = 0;
    global.N_sum[2] = 0;
    
}

void init_pinningsites()
{
    global.pinningsite_x = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_y = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_fx = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_fy = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_color = (int *)malloc(global.N_pinningsites*sizeof(int));
    global.pinningsite_direction_x = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_direction_y = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_dx_so_far = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_dy_so_far = (double *)malloc(global.N_pinningsites*sizeof(double));

    //init_pinningsites_square_lattice();
    init_pinningsites_triangular_lattice();
}

void init_pinningsites_square_lattice()
{
    int i, j, k;
    int N_rows, N_columns;
    double lattice_const;

    N_rows = (int)sqrt((double)global.N_pinningsites);
    N_columns = N_rows;

    lattice_const = global.SX / (double) N_rows;
    global.pinning_lattice_constant = lattice_const;

    k = 0;

    for(i = 0; i < N_rows; i++)
        for(j = 0; j < N_columns; j++)
        {
            global.pinningsite_x[k] = (i + 0.5) * lattice_const;
            global.pinningsite_y[k] = (j + 0.5) * lattice_const;
            global.pinningsite_dx_so_far[k] = 0.0;
            global.pinningsite_dy_so_far[k] = 0.0;
            global.pinningsite_fx[k] = 0.0;
            global.pinningsite_fy[k] = 0.0;
            k++;
        }

    printf("pinningsites initialized\n");
    printf("N_pinningsites = %d\n", k);
    printf("Square lattice\nRows = %d\nColumns = %d\nLattice constant = %.2lf\n", N_rows, N_columns, lattice_const);
}

void init_pinningsites_triangular_lattice()
{
    int i, j, k;
    int N_rows, N_columns;
    double lattice_const;

    N_columns = (int)sqrt((double)global.N_pinningsites);
    lattice_const = global.SX/(double)N_columns;
    global.pinning_lattice_constant = lattice_const;

    //keeping same number of rows and columns
    N_rows = N_columns;
    //if (N_rows%2==1) N_rows++;

    /*
    //rescaling rows rto nearest number for SY
    N_rows = (int) floor(global.SY / (lattice_const * (sqrt(3)/2.0)));
    printf("Rescaled row number to = %d\n", N_rows);
    */

    global.SY = N_rows * lattice_const * (sqrt(3) / 2.0);

    k = 0;

    for (i = 0; i < N_rows; i++)
        for (j = 0; j < N_columns; j++)
        {
            global.pinningsite_x[k] = (i + 0.25 + 0.5*(j%2) ) * lattice_const;
            global.pinningsite_y[k] = (j + 0.25) * lattice_const * (sqrt(3)/2.0);
            global.pinningsite_dx_so_far[k] = 0.0;
            global.pinningsite_dy_so_far[k] = 0.0;
            global.pinningsite_fx[k] = 0.0;
            global.pinningsite_fy[k] = 0.0;
            
            //set the color and direction of the pinningsite
            
            if (j%2==0) global.pinningsite_color[k] = 4 + i % 3;
            else global.pinningsite_color[k] = 4 + (i+2) % 3;
            
            //printf("%d %d %d\n", i, j, global.pinningsite_color[k]);
            
            //driection based on color
            switch (global.pinningsite_color[k] - 4)
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
            k++;
        }

    printf("Pinning sites initialized\n");
    printf("Triangular lattice\n\tRows: %d\n\tColumns: %d\n\tLattice constant: %.2lf\n", N_rows, N_columns, lattice_const);
    printf("\tRescaled SY to: %.2lf\n", global.SY);
}


void init_particles()
{
    global.particle_x = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_y = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_fx = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_fy = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_color = (int *)malloc(global.N_particles*sizeof(int));
    global.particle_direction_x = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_direction_y = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_dx_so_far = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_dy_so_far = (double *)malloc(global.N_particles*sizeof(double));

    //init_particles_square_lattice();
    //init_particles_triangular_lattice();
    init_particles_randomly();

    printf("Particles initialized\n");
}

//calculates the shortest distance between 2 points in a PBC configuration
//only used in the random deposition chekc which happens once at initialization
//so this is not so time crucial left the square root inside
double distance_folded_PBC(double x0, double y0, double x1, double y1)
{
    double r;
    double dx, dy;

    dx = x1 - x0;
    dy = y1 - y0;

    //PBC fold back
    //if any distance is lrger than half the box
    //the copy in the neighboring box is closer
    if (dx > global.halfSX) dx -=global.SX;
    if (dx <= -global.halfSX) dx +=global.SX;
    if (dy > global.halfSY) dx -=global.SY;
    if (dy <= -global.halfSY) dx +=global.SY;

    r = sqrt( dx*dx + dy*dy );

    return r;
}


void init_particles_randomly()
{
    int i, j;
    double x_try, y_try;
    double r_min;
    double dr;
    int overlap;
    int N_trials;

    r_min = 0.2;

    for(i = 0; i < global.N_particles; i++)
    {
        x_try = 0.0;
        y_try = 0.0;
        
        //check overlap with previous particles
        //assume there is overlap to get into the cycle
        overlap = 1;
        N_trials = 0;
        
        while ((overlap == 1) && (N_trials < global.N_particles))
        {    
            //attempt to place the particle
            x_try = global.SX * rand() / (RAND_MAX + 1.0);
            y_try = global.SY * rand() / (RAND_MAX + 1.0);
    
            //assume this was good
            overlap = 0;
            
            for(j = 0; j < i; j++)
            {
                //calculate distance
                dr = distance_folded_PBC(x_try, y_try, global.particle_x[j], global.particle_y[j]);
                if (dr < r_min)
                {
                    overlap = 1; //found overlap
                    N_trials++;
                    break; //no need to check with other particles
                }
            }
        }
            
        if (N_trials==global.N_particles)
        {
            printf("/033[1;31mCan't place particles randomly, quitting\n\033[0m");
            exit(1);
        }
        
        global.particle_x[i] = x_try;
        global.particle_y[i] = y_try;
        global.particle_dx_so_far[i] = 0.0;
        global.particle_dy_so_far[i] = 0.0;
        global.particle_fx[i] = 0.0;
        global.particle_fy[i] = 0.0;
        
        if (rand() / (RAND_MAX + 1.0) < 0.5)
        {
            global.particle_direction_x[i] = - 1.0;
            global.particle_direction_y[i] = 0.0;
            global.particle_color[i] = 2;
        } 
        else
        {
            global.particle_direction_x[i] = + 1.0;
            global.particle_direction_y[i] = 0.0;
            global.particle_color[i] = 3;
        }    
    }

    printf("Random arrangement of particles initialized\n");
}

void init_particles_square_lattice()
{
    int i, j, k;
    int N_rows, N_columns;
    double lattice_const;

    N_rows = (int)sqrt((double)global.N_particles);
    N_columns = N_rows;

    lattice_const = global.SX/(double)N_rows;
    global.pinning_lattice_constant = lattice_const;

    k = 0;

    for(i = 0; i < N_rows; i++)
        for(j = 0; j < N_columns; j++)
        {
            global.particle_x[k] = (i + 0.5) * lattice_const;
            global.particle_y[k] = (j + 0.5) * lattice_const;
            global.particle_dx_so_far[k] = 0.0;
            global.particle_dy_so_far[k] = 0.0;
            global.particle_fx[k] = 0.0;
            global.particle_fy[k] = 0.0;
            k++;
        }
    
    printf("Square lattice of particles initialized\n");
    printf("N_particles = %d\n", k);
    printf("Square lattice\nRows = %d\nColumns = %d\nLattice constant = %.2lf\n", N_rows, N_columns, lattice_const);
}

void init_particles_triangular_lattice()
{
    int i, j, k;
    int N_rows, N_columns;
    double lattice_const;

    N_columns = (int)sqrt((double)global.N_particles);
    lattice_const = global.SX/(double)N_columns;
    global.pinning_lattice_constant = lattice_const;

    //keeping same number of rows and columns
    N_rows = N_columns;
    //if (N_rows%2==1) N_rows++;

    /*
    //rescaling rows rto nearest number for SY
    N_rows = (int) floor(global.SY / (lattice_const * (sqrt(3)/2.0)));
    printf("Rescaled row number to = %d\n", N_rows);
    */

    global.SY = N_rows * lattice_const * (sqrt(3)/2.0);
    global.halfSY = global.SY/2.0;
    printf("Rescaled SY to = %.2lf\n", global.SY);

    k=0;

    for(i=0;i<N_rows;i++)
        for(j=0;j<N_columns;j++)
            {
            global.particle_x[k] = (i + 0.25 + 0.5*(j%2) ) * lattice_const;
            global.particle_y[k] = (j + 0.25) * lattice_const * (sqrt(3)/2.0);
            global.particle_dx_so_far[k] = 0.0;
            global.particle_dy_so_far[k] = 0.0;
            global.particle_fx[k] = 0.0;
            global.particle_fy[k] = 0.0;
            
            //set the color and direction of the particle
            
            if (j%2==0) global.particle_color[k] = 2 + i % 3;
            else global.particle_color[k] = 2 + (i+2) % 3;
            
            //printf("%d %d %d\n", i, j, global.particle_color[k]);
            
            //driection based on color
            switch (global.particle_color[k]-2)
                {
                case 0: {
                        global.particle_direction_x[k] = - 0.5;
                        global.particle_direction_y[k] = - sqrt(3)/2.0;
                        break;
                        }
                case 1: {
                        global.particle_direction_x[k] = 1.0;
                        global.particle_direction_y[k] = 0.0;
                        break;
                        }
                case 2: {
                        global.particle_direction_x[k] = - 0.5;
                        global.particle_direction_y[k] = + sqrt(3)/2.0;
                        break;
                        }
                }
            k++;
            }

    printf("Particles initialized\n");
    printf("N_particles = %d\n", k);
    printf("Triangular lattice\nRows = %d\nColumns = %d\nLattice constant = %.2lf\n", N_rows, N_columns, lattice_const);

    FILE *test;

    test = fopen("testfile1.txt", "wt");
    for(i=0;i<global.N_particles;i++)
        if (global.particle_color[i]==0) fprintf(test, "%lf %lf\n", global.particle_x[i], global.particle_y[i]);
    fclose(test);
    test = fopen("testfile2.txt", "wt");
    for(i=0;i<global.N_particles;i++)
        if (global.particle_color[i]==1) fprintf(test, "%lf %lf\n", global.particle_x[i], global.particle_y[i]);
    fclose(test);
    test = fopen("testfile3.txt", "wt");
    for(i=0;i<global.N_particles;i++)
        if (global.particle_color[i]==2) fprintf(test, "%lf %lf\n", global.particle_x[i], global.particle_y[i]);
    fclose(test);
}

void init_files()
{
    global.moviefile = fopen(global.moviefile_name, "wb");
    if (global.moviefile == NULL)
    {
        printf("\033[1;31mCould not create/open movie file\033[0m\n");
        exit(2);
    }

    global.statisticsfile = fopen(global.statisticsfile_name, "wt");
    if (global.statisticsfile == NULL)
    {
        printf("\033[1;31mCould not create/open statistics file\033[0m\n");
        exit(2);
    }
}

char* substr(const char* from, int start, int count)
{
    char* result = (char *)malloc(count + 1);
    int i, n;
    for (i = start, n = 0; i < start + count; i++, n++)
    {
        result[n] = from[i];
    }
    result[n] = '\0';
    return result;
}
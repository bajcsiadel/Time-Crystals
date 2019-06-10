//
//  initializer.h
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef initializer_h
#    define initializer_h

#    include <stdio.h>

#    define PROJECT_NAME "Time-Crystals"

// inital values for the simulation

#    define __SX 72.0
#    define __SY 72.0

#    define __PINNINGSITE_R				0.2
#    define __PINNINGSITE_FORCE			0.6
#    define __PINNINGSITE_DRIVING_FORCE	0.5
#    define __N_PINNINGSITES			567
#    define __PINNING_DIRECTION_CHANGE	1

#    define __N_PARTICLES				1000
#    define __PARTICLE_DRIVING_FORCE	0.2
#    define __PARTICLE_PARTICLE_SCREENING_LENGTH 4.0

#    define __DT			0.001
#    define __STAT_TIME		100
#    define __MOVIE_TIME	100
#    define __TOTAL_TIME	100000
#    define __ECHO_TIME		1000

/**
 * Checking if the given filename correspods to the given regular expression.
 * 
 * regex_str - regular expression
 * to_test   - string to test if it matches regex
 * 
 * Return value:
 * 1 - to_test matches the regular expression
 * 0 - otherwise
*/
int     check_file_name(const char *, const char *);
/**
 * Calculates the shortest distance between 2 points in a PBC configuration. 
 * Only used in the random deposition check which happens once at initialization
 * so this is not so time crucial left the square root inside
 * 
 * (x0, y0) - first point's coordinates
 * (x1, y1) - second point's coordinates
 * 
 * Return value: the distance between the two given point
*/
double  distance_folded_PBC(double x0, double y0, double x1, double y1);
/**
 * Checking if file given as parameter exists.
 * 
 * Return value:
 * 1 - file exists
 * 0 - file does not exist
*/
int     file_exists(const char *filename);
/**
 * Getting current working directory.
 * The value will be returned in the parameter (current_working_dir).
*/
void    get_current_woring_dir(char *current_working_dir);
/**
 * Defining the path to the result directory. PROJECT_NAME macro should be set to the name of the root directory of the project.
*/
void    get_path_to_result_folder(char *path, size_t path_size);
/**
 * Get token from JSON dictionary.
 * start    - start position in the dictionary
 * end      - end position in the dictionary. The token will be searched between these two values.
 * token    - searched token name
 * 
 * Return value: start position of the token if it is in the JSON object, otherwise -1.
*/
int     get_token(int, int, const char *);
/**
 * Reading data from JSON object and initializing simulation system data.
*/
void    init_data();
/**
 * Allocating arrays needed to describe the behavior of particles.
*/
void    init_particles();
/**
 * Generating particle's properties (coordinates, color - defines direction forces as well) at random.
*/
void    init_particles_randomly();
/**
 * Arrange particles in a regular square lattice.
*/
void    init_particles_square_lattice();
/**
 * Arrange particles in a triangular lattice.
*/
void    init_particles_triangular_lattice();
/**
 * Allocating arrays needed to describe the behavior of pinningsites.
*/
void    init_pinningsites();
/**
 * Arrange pinningsites in a regular square lattice.
*/
void    init_pinningsites_square_lattice();
/**
 * Arrange pinningsites in a triangular lattice.
*/
void    init_pinningsites_triangular_lattice();
/**
 * Initializing Verlet list and pinningsite grid.
*/
void    init_simulation();
/**
 * Open movie and statistics files.
*/
void    open_files();
/** 
 * Reads the parameter file adjusts the parameters accordingly
 * filename - parameter file
 * 
 * Return value: length of the JSON object string
*/
int     read_init_file(const char *);
/**
 * Setting seed to random. If is set in the console then that value else time(NULL).
 * seed - console parameter
*/
void    set_seed(const char *);
/**
 * The substring is the portion of the object that starts at character position start and spans count characters (or until the end of the string, whichever comes first).
 * 
 * Return value: newly constructed string object with its value initialized to a copy of a substring of this object.
*/
char   *substr(const char *from, int start, int count);

#endif /* initializer_h */

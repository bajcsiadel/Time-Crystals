//
//  initializer.c
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#include <math.h>				// sqrt
#include <regex.h>
#include <string.h>				// strlen, srtstr, strncpy, strncat, strncmp, strchr
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>				// isspace
#include <time.h>				// time
#include <stdarg.h>				// for void something(int arg, ...)

#include "initializer.h"
#include "globaldata.h"
#include "utils.h"


// depending on operating system defining the function file_exists and get_current_dir
// http://www.codebind.com/cpp-tutorial/c-get-current-directory-linuxwindows/
#ifdef WINDOWS
#include <direct.h>
#include <windows.h>
#define get_current_dir _getcwd
int file_exists(const char *filename)
{
	WIN32_FIND_DATA FindFileData;
	HANDLE handle = FindFirstFile(filename, &FindFileData);
	int found = handle != INVALID_HANDLE_VALUE;
	if (found)
	{
		FindClose(handle);
		return 1;
	}
	else
		return 0;
}
#else // LINUX or MAC
#include <unistd.h>
#define get_current_dir getcwd
int file_exists(const char *filename)
{
	if (access(filename, F_OK) != -1)
		return 1;
	else
		return 0;
}
#endif

void get_current_woring_dir(char *current_working_dir)
{
	char buff[FILENAME_MAX];
	// get_current_dir calls getcwd or _getcwd depending on the operating system the program runs
	get_current_dir(buff, FILENAME_MAX);
	// copying the result to the parameter
	memcpy(current_working_dir, buff, strlen(buff));
}

void get_path_to_result_folder(char *path, size_t path_size)
{
	char *current_working_dir, *time_crystals;

	current_working_dir = (char *) malloc(255);
	get_current_woring_dir(current_working_dir);
	// find project root dir in current working path
	time_crystals = strstr(current_working_dir, PROJECT_NAME);
	// setting path to empty string
	strncpy(path, "", 1);
	// end -  address of last position is path variable
	char *buff, *const end = path + path_size;

	// at the beginning 
	buff = path;
	for (; *time_crystals != '\0'; time_crystals++)
	{
		if (*time_crystals == '\\' || *time_crystals == '/')
		{
			buff += snprintf(buff, end - buff, "../");
		}
	}

	snprintf(buff, end - buff, "results/");
	free(current_working_dir);
}

void set_seed(const char *seed)
{
	unsigned int s;
	if (seed)
		// convert char* to int
		s = atoi(seed);
	else
		s = time(NULL);

	// set random seed
	srand(s);
	print_log(stdout, NOTE, __FILE__, __LINE__, NULL, 1,
			  "Seed value: %u\n", s);
}

int read_init_file(const char *filename)
{
	FILE *f = fopen(filename, "rt");
	char tmp[255];
	int i, n;
	size_t len_tmp, len_json;

	n = 3;
	global.JSON = (char *) malloc(n);
	global.JSON[0] = '\0';
	len_json = 1;
	if (f != NULL)
	{
		print_log(stdout, NOTE, __FILE__, __LINE__,
				  "Initializing data from the parameterfile", 1,
				  "Opened parameter file %s\n", filename);
		while (!feof(f))
		{
			// get 255 character from f
			fgets(tmp, 255, f);
			len_tmp = strlen(tmp);
			len_json += len_tmp;
			// widening array if exceeds the current allocated length
			if (len_json > n)
			{
				n += 255;
				global.JSON = (char *) realloc(global.JSON, n);
			}
			// concatenating the newly read string to the JSON object string
			strncat(global.JSON, tmp, len_tmp);
			global.JSON[len_json] = '\0';
		}

		char first = '\0';
		// depending on operating system it can happen that the last character of the object, usually '}',
		// will be duplicated (ATTENTION do NOT write an object as the last element in JSON), therefore 
		// iterating from the end of the string toward is beginning
		for (i = strlen(global.JSON); i >= 0; i--)
		{
			// checking if current character is not white character
			if (!isspace(global.JSON[i]))
				// remembering the first not-white character at the end of the string
				if (first == '\0')
					first = global.JSON[i];
				else
				{
					// if the next not-white character is the same then cut the tail after this character.
					if (first == global.JSON[i])
						global.JSON[i + 1] = '\0';
					break;
				}
		}
	}
	else
	{
		// if file could not be open then initialize JSON object string with a string describing an empty object
		global.JSON = "{}";
		print_log(stdout, ERROR, __FILE__, __LINE__, "Could not open file.",
				  1, "Failed to open parameter file: %s\n", filename);
	}
	fclose(f);
	return n;
}

int check_file_name(const char *regex_str, const char *to_test)
{
	regex_t regexCompiled;
	regmatch_t pmatch[1];
	int result;

	// creating regular expression from string
	if (regcomp(&regexCompiled, regex_str, REG_EXTENDED))
	{
		print_log(stdout, ERROR, __FILE__, __LINE__,
				  "Could not create regex!", 0);
		return 0;
	}

	// running the regular expression on to_test, and checking if it is totally matches
	result = !regexec(&regexCompiled, to_test, 1, pmatch, 0)
		? ((pmatch[0].rm_so == 0) && (pmatch[0].rm_eo == strlen(to_test)))
		: 0;
	// deallocating the regular expression
	regfree(&regexCompiled);
	return result;
}

int get_token(int start, int end, const char *token)
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
				char *tmp = (char *) malloc(l1);
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
	jsmn_parser p;
	int n, r, x, y, object_end;
	size_t len, regex_len, s;
	char *regex_str, *to_test, *path, *dot;

	// maximum length of a regulat expression
	regex_len = 36;

	// length of the object read from the parameter file
	len = (size_t) strlen(global.JSON);
	// initializing the JSON object
	jsmn_init(&p);
	// try (last parameter is 0) to parse a JSON data string into and array of tokens
	// n - will be the length of this array (now we know the amount of space needed)
	n = jsmn_parse(&p, global.JSON, len, global.t, 0);
	if (n < 0)
	{
		print_log(stdout, ERROR, __FILE__, __LINE__, "Some error occured!",
				  0);
		return;
	}
	else
	{
		global.t = (jsmntok_t *) malloc(n * sizeof(jsmntok_t));
	}
	jsmn_init(&p);
	// run JSON parser. It parses a JSON data string into and array of tokens, each
	// describing a single JSON object.
	r = jsmn_parse(&p, global.JSON, len, global.t, n);
	// assume the top-level element is an object
	if (r < 1 || global.t[0].type != JSMN_OBJECT)
	{
		print_log(stdout, ERROR, __FILE__, __LINE__, "JSON object not found!",
				  0);
		return;
	}

	// find token "simulation_box"
	if ((x = get_token(1, r, "simulation_box")) >= 0)
	{
		// assume the value corresponding to this key is an object
		if (global.t[x + 1].type == JSMN_OBJECT)
		{
			// calculate the end of the object
			// token.size gives the number of tokens in the object so it has to be multipied
			// by 2 because a token consists in a key and a value 
			object_end = x + 1 + 2 * global.t[x + 1].size;
			// find token "SX", the simulation box's width
			if ((y = get_token(x + 2, object_end, "SX")) >= 0)
			{
				global.SX = atof(substr(global.JSON, global.t[y + 1].start,
										global.t[y + 1].end - global.t[y +
																	   1].
										start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'SX' not found in the parameter file", 1,
						  "Default value setted.\n");
				global.SX = __SX;
			}

			// find token "SY", the simulation box's height
			if ((y = get_token(x + 2, object_end, "SY")) >= 0)
			{
				global.SY = atof(substr(global.JSON, global.t[y + 1].start,
										global.t[y + 1].end - global.t[y +
																	   1].
										start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'SY' not found in the parameter file", 1,
						  "Default value setted.\n");
				global.SY = __SY;
			}
		}
	}
	else
	{
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'simulation_box' not found in the parameter file", 1,
				  "Default values setted.");
		// if "simulation_box" not found then set default values 
		global.SX = __SX;
		global.SY = __SY;
	}
	print_log(stdout, NOTE, __FILE__, __LINE__, "Simulation box", 2,
			  "SX: %f\n", global.SX, "SY: %f\n", global.SY);
	// calculate mid-point of the simulation box
	global.halfSX = global.SX / 2.0;
	global.halfSY = global.SY / 2.0;

	// find token "pinningsite"
	if ((x = get_token(1, r, "pinningsite")) >= 0)
	{
		// assume the value corresponding to this key is an object
		if (global.t[x + 1].type == JSMN_OBJECT)
		{
			// calculate the end of the object
			object_end = x + 1 + 2 * global.t[x + 1].size;
			// find token "pinningsite_R", radius of the pinningsites
			if ((y = get_token(x + 2, object_end, "pinningsite_R")) >= 0)
			{
				global.pinningsite_R =
					atof(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'pinningsite_R' not found in the parameter file",
						  1, "Default value setted.\n");
				global.pinningsite_R = __PINNINGSITE_R;
			}

			// find token "pinningsite_force", the internal force of the pinningsites
			if ((y = get_token(x + 2, object_end, "pinningsite_force")) >= 0)
			{
				global.pinningsite_force =
					atof(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'pinningsite_force' not found in the parameter file",
						  1, "Default value setted.\n");
				global.pinningsite_force = __PINNINGSITE_FORCE;
			}

			// find token "pinningsite_driving_force", the external force affecting pinningsites
			if ((y =
				 get_token(x + 2, object_end, "pinning_driving_force")) >= 0)
			{
				global.pinning_driving_force =
					atof(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'pinning_driving_force' not found in the parameter file",
						  1, "Default value setted.\n");
				global.pinning_driving_force = __PINNINGSITE_DRIVING_FORCE;
			}

			// find token "N_pinningsites", number of pinningsites in the system
			if ((y = get_token(x + 2, object_end, "N_pinningsites")) >= 0)
			{
				global.N_pinningsites =
					atoi(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'N_pinningsites' not found in the parameter file",
						  1, "Default value setted.\n");
				global.N_pinningsites = __N_PINNINGSITES;
			}

			// find token "pinning_direction change"
			if ((y =
				 get_token(x + 2, object_end,
						   "pinning_direction_change")) >= 0)
			{
				global.pinning_direction_change =
					atoi(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'pinning_direction_change' not found in the parameter file",
						  1, "Default value setted.\n");
				global.pinning_direction_change = __PINNING_DIRECTION_CHANGE;
			}
		}
	}
	else
	{
		// if "pinningsite" not found then set default values
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'pinningsite' not found in the parameter file", 1,
				  "Default values setted.");
		global.pinningsite_R = __PINNINGSITE_R;
		global.pinningsite_force = __PINNINGSITE_FORCE;
		global.pinning_driving_force = __PINNINGSITE_DRIVING_FORCE;
		global.N_pinningsites = __N_PINNINGSITES;
		global.pinning_direction_change = __PINNING_DIRECTION_CHANGE;
	}
	print_log(stdout, NOTE, __FILE__, __LINE__, "Pinningsites", 5,
			  "Number: %d\n", global.N_pinningsites,
			  "Radius: %f\n", global.pinningsite_R,
			  "Force: %f\n", global.pinningsite_force,
			  "Driving force: %f\n", global.pinning_driving_force,
			  "Direction change: %d\n", global.pinning_direction_change);

	// find token "particle"
	if ((x = get_token(1, r, "particle")) >= 0)
	{
		// assume the value corresponding to this key is an object
		if (global.t[x + 1].type == JSMN_OBJECT)
		{
			object_end = x + 1 + 2 * global.t[x + 1].size;
			// find token "N_particles", number of particles in the system
			if ((y = get_token(x + 2, object_end, "N_particles")) >= 0)
			{
				global.N_particles =
					atoi(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'N_particles' not found in the parameter file", 1,
						  "Default value setted.\n");
				global.N_particles = __N_PARTICLES;
			}

			// find token "particle_driving_force", external force affecting the particles
			if ((y =
				 get_token(x + 2, object_end, "particle_driving_force")) >= 0)
			{
				global.particle_driving_force =
					atof(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'particle_driving_force' not found in the parameter file",
						  1, "Default value setted.\n");
				global.particle_driving_force = __PARTICLE_DRIVING_FORCE;
			}

			// find token "particle_particle_screening_length", range square in which particles affect each other
			if ((y =
				 get_token(x + 2, object_end,
						   "partile_particle_screening_length")) >= 0)
			{
				global.partile_particle_screening_length =
					atof(substr(global.JSON, global.t[y + 1].start,
								global.t[y + 1].end - global.t[y + 1].start));
			}
			else
			{
				// if not defined set default value
				print_log(stdout, WARNING, __FILE__, __LINE__,
						  "'partile_particle_screening_length' not found in the parameter file",
						  1, "Default value setted.\n");
				global.partile_particle_screening_length =
					__PARTICLE_PARTICLE_SCREENING_LENGTH;
			}
		}
	}
	else
	{
		// if "particle" not found then set default values
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'particle' not found in the parameter file", 1,
				  "Default values setted.");
		global.N_particles = __N_PARTICLES;
		global.particle_driving_force = __PARTICLE_DRIVING_FORCE;
		global.partile_particle_screening_length =
			__PARTICLE_PARTICLE_SCREENING_LENGTH;
	}

	// calculating reciprocal of particle particle screening length (it will be useful when calculating
	// force between two particles)
	global.partile_particle_screening_wavevector =
		1.0 / global.partile_particle_screening_length;

	print_log(stdout, NOTE, __FILE__, __LINE__, "Particles", 4,
			  "Number: %d\n", global.N_particles,
			  "Driving force: %f\n", global.particle_driving_force,
			  "Particle particle screening length: %f\n",
			  global.partile_particle_screening_length,
			  "Particle particle screening wavevector: %lf\n",
			  global.partile_particle_screening_wavevector);

	// find token "dt", timestep between iterations in the system
	if ((x = get_token(1, r, "dt")) >= 0)
	{
		global.dt = atof(substr(global.JSON, global.t[x + 1].start,
								global.t[x + 1].end - global.t[x + 1].start));
	}
	else
	{
		// if not defined set default value
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'dt' not found in the parameter file", 1,
				  "Default value setted.\n");
		global.dt = __DT;
	}

	// find token "statistics_time", interval of writing statistics file
	if ((x = get_token(1, r, "statistics_time")) >= 0)
	{
		global.statistics_time =
			atoi(substr(global.JSON, global.t[x + 1].start,
						global.t[x + 1].end - global.t[x + 1].start));
	}
	else
	{
		// if not defined set default value
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'statistics_time' not found in the parameter file", 1,
				  "Default value setted.\n");
		global.statistics_time = __STAT_TIME;
	}

	// find token "movie_time", interval of writing movie
	if ((x = get_token(1, r, "movie_time")) >= 0)
	{
		global.movie_time = atoi(substr(global.JSON, global.t[x + 1].start,
										global.t[x + 1].end - global.t[x +
																	   1].
										start));
	}
	else
	{
		// if not defined set default value
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'movie_time' not found in the parameter file", 1,
				  "Default value setted.\n");
		global.movie_time = __MOVIE_TIME;
	}

	// find token "total_time", total iteration on the system
	if ((x = get_token(1, r, "total_time")) >= 0)
	{
		global.total_time = atoi(substr(global.JSON, global.t[x + 1].start,
										global.t[x + 1].end - global.t[x +
																	   1].
										start));
	}
	else
	{
		// if not defined set default value
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'total_time' not found in the parameter file", 1,
				  "Default value setted.\n");
		global.total_time = __TOTAL_TIME;
	}

	// find token "echo_time", interval of writing where the simulation is to the screen
	if ((x = get_token(1, r, "echo_time")) >= 0)
	{
		global.echo_time = atoi(substr(global.JSON, global.t[x + 1].start,
									   global.t[x + 1].end - global.t[x +
																	  1].
									   start));
	}
	else
	{
		// if not defined set default value
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'echo_time' not found in the parameter file", 1,
				  "Default value setted.\n");
		global.echo_time = __ECHO_TIME;
	}

	print_log(stdout, NOTE, __FILE__, __LINE__, NULL, 1,
			  "Timestep: %f\n", global.dt);
	print_log(stdout, NOTE, __FILE__, __LINE__, "Echo times", 3,
			  "To screen: %d\n", global.echo_time,
			  "To statistics file: %d\n", global.statistics_time,
			  "To movie file: %d\n", global.movie_time);
	print_log(stdout, NOTE, __FILE__, __LINE__, NULL, 1,
			  "Total running time: %d\n", global.total_time);

	// get current time
	time_t now = time(NULL);
	char buff[16];
	// creating a string (buff) containing the current time in a given format
	// <year><month><day>_<hour><minute><second>
	strftime(buff, 16, "%Y%m%d_%H%M%S", localtime(&now));

	// find token "filename", name of output files
	// just name without extension
	// this will be the name for the statistics and movie file as well
	if ((x = get_token(1, r, "filename")) >= 0)
	{
		regex_str = (char *) malloc(regex_len);
		// create a regulare expression which matches words containing letters (small and big ones),
		// numbers or underscore characters
		strncpy(regex_str, "^[a-z0-9A-Z_]+$", (size_t) 15);
		regex_str[15] = '\0';
		// get value for the key "filename"
		to_test =
			substr(global.JSON, global.t[x + 1].start,
				   global.t[x + 1].end - global.t[x + 1].start);
		if (!check_file_name(regex_str, to_test))
		{
			// if the filename contains any other character than in the list above
			// then find for dot and replace it with null character (cutting the remaining after dot)
			dot = strchr(to_test, '.');
			if (dot != NULL)
				*dot = '\0';
		}
		if (!check_file_name(regex_str, to_test))
		{
			// if filename still not matches the regular expression then setting default value to it
			// default: <timestamp>_result_<pinningsite_force_integer_part><pinningsite_force_decimal_part>
			print_log(stdout, WARNING, __FILE__, __LINE__,
					  "'filename' contains characters other than a-z, A-Z, 0-9 or _",
					  1, "Default value setted.\n");
			snprintf(to_test, 30, "%s_result_%d%d", buff,
					 (int) global.pinningsite_force,
					 (int) (global.pinningsite_force -
							(int)global.pinningsite_force) *100);
		}
		free(regex_str);
	}
	else
	{
		// if not defined set default value
		// default: <timestamp>_result_<pinningsite_force_integer_part><pinningsite_force_decimal_part>
		print_log(stdout, WARNING, __FILE__, __LINE__,
				  "'filename' not found in the parameter file", 1,
				  "Default value setted.\n");
		to_test = (char *) malloc(30);
		snprintf(to_test, 30, "%s_result_%d_%d", buff,
				 (int) global.pinningsite_force,
				 (int) (global.pinningsite_force -
						(int)global.pinningsite_force) *100);
		to_test[6] = '\0';
	}

	path = (char *) malloc(255);
	get_path_to_result_folder(path, 255);

	// try to make the whole path to movie file
	// len - length of the concatenated string
	len = snprintf(NULL, 0, "%smovies/%s.mvi", path, to_test);
	global.moviefile_name = (char *) malloc(len + 1);
	// make the concatenation
	snprintf(global.moviefile_name, len + 1, "%smovies/%s.mvi", path,
			 to_test);

	// try to make the whole path to statistics file
	len = snprintf(NULL, 0, "%sstats/%s.txt", path, to_test);
	global.statisticsfile_name = (char *) malloc(len + 1);
	snprintf(global.statisticsfile_name, len + 1, "%sstats/%s.txt", path,
			 to_test);

	// if one of the result files already exists then creating a unique filename by adding a timestemp to it
	/*
	if (file_exists(global.moviefile_name)
		|| file_exists(global.statisticsfile_name))
	{
		// try to create concatenated string by adding timestamp to the base filename
		len = snprintf(NULL, 0, "%s_%s", to_test, buff);
		to_test = (char *) realloc(to_test, len + 1);
		snprintf(&to_test[strlen(to_test)], 16, "_%s", buff);

		// try to concatenate the whole path with the new filename
		len = snprintf(NULL, 0, "%smovies/%s.mvi", path, to_test);
		global.moviefile_name =
			(char *) realloc(global.moviefile_name, len + 1);
		snprintf(global.moviefile_name, len + 1, "%smovies/%s.mvi", path,
				 to_test);

		// try to concatenate the whole path with the new filename
		len = snprintf(NULL, 0, "%sstats/%s.txt", path, to_test);
		global.statisticsfile_name =
			(char *) realloc(global.statisticsfile_name, len + 1);
		snprintf(global.statisticsfile_name, len + 1, "%sstats/%s.txt", path,
				 to_test);
	}*/

	free(to_test);

	print_log(stdout, NOTE, __FILE__, __LINE__, "File names", 2,
			  "Moviefile: %s\n", global.moviefile_name,
			  "Statfile:  %s\n", global.statisticsfile_name);

	free(path);
	free(global.t);
}

//init the common variables in the simulation
void init_simulation()
{
	global.Verlet_cutoff_distance =
		1.5 * global.partile_particle_screening_length;
	global.Verlet_cutoff_distance_squared =
		global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
	global.Verlet_intershell_squared =
		global.Verlet_cutoff_distance -
		global.partile_particle_screening_length;
	global.Verlet_intershell_squared = global.Verlet_intershell_squared / 2.0;
	global.Verlet_intershell_squared *= global.Verlet_intershell_squared;

	print_log(stdout, NOTE, __FILE__, __LINE__, "Verlet properties", 4,
			  "Verlet cutoff distance = %.2lf\n",
			  global.Verlet_cutoff_distance,
			  "Verlet cutoff distance squared = %.2lf\n",
			  global.Verlet_cutoff_distance_squared,
			  "Half of Verlet intershell distance = %.2lf\n",
			  sqrt(global.Verlet_intershell_squared),
			  "Half of Verlet intershell distance squared = %.2lf\n",
			  global.Verlet_intershell_squared);

	// zero everything so rebuild Verlet can find this the first time
	global.Verletlisti = NULL;
	global.Verletlistj = NULL;
	global.N_Verlet = 0;
	global.N_Verlet_max = 0;

	// zero everythiong so rebuild_pinning_grid can find this the first time
	global.pinningsite_grid = NULL;
	global.Nx_pinningsite_grid = 0;
	global.Ny_pinningsite_grid = 0;
	global.max_pinningsite_per_grid = 0;
}

void init_pinningsites()
{
	global.pinningsite_x =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_y =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_fx =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_fy =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_color =
		(int *) malloc(global.N_pinningsites * sizeof(int));
	global.pinningsite_direction_x =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_direction_y =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_dx_so_far =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.pinningsite_dy_so_far =
		(double *) malloc(global.N_pinningsites * sizeof(double));
	global.particle_in_pinningsite =
		(unsigned int *) malloc(global.N_pinningsites * sizeof(unsigned int));

	//init_pinningsites_square_lattice();
	init_pinningsites_triangular_lattice();
}

void init_pinningsites_square_lattice()
{
	int i, j, k;
	int N_rows, N_columns;
	double lattice_const;

	// calculate lattice size (row x column)
	N_rows = N_columns = (int) sqrt((double) global.N_pinningsites);

	// calculate row/column width
	global.pinning_lattice_constant = lattice_const =
		global.SX / (double) N_rows;

	k = 0;
	for (i = 0; i < N_rows; i++)
		for (j = 0; j < N_columns; j++)
		{
			// calculating pinningsite's coordinates
			global.pinningsite_x[k] = (i + 0.5) * lattice_const;
			global.pinningsite_y[k] = (j + 0.5) * lattice_const;
			// set other properties to 0
			global.pinningsite_dx_so_far[k] = 0.0;
			global.pinningsite_dy_so_far[k] = 0.0;
			global.pinningsite_fx[k] = 0.0;
			global.pinningsite_fy[k] = 0.0;
			k++;
		}

	global.N_pinningsites = k;
	print_log(stdout, NOTE, __FILE__, __LINE__,
			  "Square lattice of pinningsites initialized", 4,
			  "N_pinningsites = %d\n", k, "Rows = %d\n", N_rows,
			  "Columns = %d\n", N_columns, "Lattice constant = %.2lf\n",
			  lattice_const);
}

void init_pinningsites_triangular_lattice()
{
	int i, j, k;
	int N_rows, N_columns;
	double lattice_const;

	// calculate column numbers in triangular lattice
	N_columns = (int) sqrt((double) global.N_pinningsites);
	// calculate row width
	global.pinning_lattice_constant = lattice_const =
		global.SX / (double) N_columns;

	// keeping same number of rows and columns
	N_rows = N_columns;

	// rescale system size depending on triangular lattice size
	global.SY = N_rows * lattice_const * (sqrt(3) / 2.0);

	k = 0;
	for (i = 0; i < N_rows; i++)
		for (j = 0; j < N_columns; j++)
		{
			// calculate pinningsite coordinates
			global.pinningsite_x[k] =
				(i + 0.25 + 0.5 * (j % 2)) * lattice_const;
			global.pinningsite_y[k] =
				(j + 0.25) * lattice_const * (sqrt(3) / 2.0);
			// set other properties to 0
			global.pinningsite_dx_so_far[k] = 0.0;
			global.pinningsite_dy_so_far[k] = 0.0;
			global.pinningsite_fx[k] = 0.0;
			global.pinningsite_fy[k] = 0.0;
			global.particle_in_pinningsite[k] = 0;

			//set the color and direction of the pinningsite
			if (j % 2 == 0)
				global.pinningsite_color[k] = 4 + i % 3;
			else
				global.pinningsite_color[k] = 4 + (i + 2) % 3;

			//driection based on color
			switch (global.pinningsite_color[k] - 4)
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
			k++;
		}

	global.N_pinningsites = k;
	print_log(stdout, NOTE, __FILE__, __LINE__,
			  "Triangular lattice of pinningsites initialized", 5,
			  "N_pinningsites = %d\n", k, "Rows: %d\n", N_rows,
			  "Columns: %d\n", N_columns, "Lattice constant: %.2lf\n",
			  lattice_const, "Rescaled SY to: %.2lf\n", global.SY);
}


void init_particles()
{
	global.particle_x =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_y =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_fx =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_fy =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_color = (int *) malloc(global.N_particles * sizeof(int));
	global.particle_direction_x =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_direction_y =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_dx_so_far =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_dy_so_far =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_all_dx =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_all_dy =
		(double *) malloc(global.N_particles * sizeof(double));
	global.particle_all_dr2 =
		(double *) malloc(global.N_particles * sizeof(double));

	// init_particles_square_lattice();
	// init_particles_triangular_lattice();
	init_particles_randomly();

	print_log(stdout, NOTE, __FILE__, __LINE__, "Particles initialized", 0);
}

double distance_folded_PBC(double x0, double y0, double x1, double y1)
{
	double r;
	double dx, dy;

	dx = x1 - x0;
	dy = y1 - y0;

	//PBC fold back
	//if any distance is lrger than half the box
	//the copy in the neighboring box is closer
	if (dx > global.halfSX)
		dx -= global.SX;
	if (dx <= -global.halfSX)
		dx += global.SX;
	if (dy > global.halfSY)
		dx -= global.SY;
	if (dy <= -global.halfSY)
		dx += global.SY;

	r = sqrt(dx * dx + dy * dy);

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

	for (i = 0; i < global.N_particles; i++)
	{
		// check overlap with previous particles
		// assume there is overlap to get into the loop
		overlap = 1;
		// there was not any try, to set particle, yet
		N_trials = 0;

		while ((overlap == 1) && (N_trials < global.N_particles))
		{
			// attempt to place the particle
			// rand() / (RAND_MAX + 1.0) -> generates number between [0, 1]
			x_try = global.SX * rand() / (RAND_MAX + 1.0);
			y_try = global.SY * rand() / (RAND_MAX + 1.0);

			// assume this was good
			overlap = 0;

			for (j = 0; j < i; j++)
			{
				// calculate distance
				dr = distance_folded_PBC(x_try, y_try, global.particle_x[j],
										 global.particle_y[j]);
				if (dr < r_min)
				{
					// found overlap
					overlap = 1;
					// inceasing trialings number
					N_trials++;
					// no need to check with other particles, breaks from the for loop
					break;
				}
			}
		}

		// if could not find a spot in N_particle tries then quit, particles can not be setted with random values
		if (N_trials == global.N_particles)
		{
			print_log(stdout, ERROR, __FILE__, __LINE__,
					  "Can't place particles randomly", 0);
			exit(1);
		}

		// set particles coordinate to the ones tested above
		global.particle_x[i] = x_try;
		global.particle_y[i] = y_try;
		// this is the system states before the simulation starts, thus particles have no travelled distance
		global.particle_dx_so_far[i] = 0.0f;
		global.particle_dy_so_far[i] = 0.0f;
		global.particle_all_dx[i] = 0.0f;
		global.particle_all_dy[i] = 0.0f;
		global.particle_all_dr2[i] = 0.0f;
		// and there is no force affecting them
		global.particle_fx[i] = 0.0f;
		global.particle_fy[i] = 0.0f;

		// generating another numer (50%-50% is the appearance of both colors)
		if (rand() / (RAND_MAX + 1.0) < 0.5)
		{
			// setting color to 2 and for this particles there will be an exernal force which pulls them to the left
			global.particle_direction_x[i] = -1.0f;
			global.particle_direction_y[i] = 0.0f;
			global.particle_color[i] = 2;
		}
		else
		{
			// setting color to 3 and for this particles there will be an exernal force which pulls them to the right
			global.particle_direction_x[i] = +1.0f;
			global.particle_direction_y[i] = 0.0f;
			global.particle_color[i] = 3;
		}
	}

	print_log(stdout, NOTE, __FILE__, __LINE__,
			  "Random arrangement of particles initialized", 0);
}

void init_particles_square_lattice()
{
	int i, j, k;
	int N_rows, N_columns;
	double lattice_const;

	// calculate row and column number in the lattice
	N_rows = N_columns = (int) sqrt((double) global.N_particles);

	// calculate width of a row/column
	global.pinning_lattice_constant = lattice_const =
		global.SX / (double) N_rows;

	k = 0;
	for (i = 0; i < N_rows; i++)
		for (j = 0; j < N_columns; j++)
		{
			// calculate coordinates of particles
			global.particle_x[k] = (i + 0.5) * lattice_const;
			global.particle_y[k] = (j + 0.5) * lattice_const;
			global.particle_dx_so_far[k] = 0.0;
			global.particle_dy_so_far[k] = 0.0;
			global.particle_all_dx[k] = 0.0;
			global.particle_all_dy[k] = 0.0;
			global.particle_all_dr2[i] = 0.0;
			global.particle_fx[k] = 0.0;
			global.particle_fy[k] = 0.0;
			k++;
		}

	global.N_particles = k;
	print_log(stdout, NOTE, __FILE__, __LINE__,
			  "Square lattice of particles initialized", 2,
			  "N_particles = %d\n", k,
			  "Rows = %d\nColumns = %d\nLattice constant = %.2lf\n", N_rows,
			  N_columns, lattice_const);
}

void init_particles_triangular_lattice()
{
	int i, j, k;
	int N_rows, N_columns;
	double lattice_const;

	N_columns = (int) sqrt((double) global.N_particles);
	// calculate width of the row/column
	global.pinning_lattice_constant = lattice_const =
		global.SX / (double) N_columns;

	// keeping same number of rows and columns
	N_rows = N_columns;

	// recalculating system's height depending on the triangular lattice
	global.SY = N_rows * lattice_const * (sqrt(3) / 2.0);
	global.halfSY = global.SY / 2.0;

	k = 0;
	for (i = 0; i < N_rows; i++)
		for (j = 0; j < N_columns; j++)
		{
			// calculate coordinates
			global.particle_x[k] = (i + 0.25 + 0.5 * (j % 2)) * lattice_const;
			global.particle_y[k] =
				(j + 0.25) * lattice_const * (sqrt(3) / 2.0);
			// setting other properties to 0
			global.particle_dx_so_far[k] = 0.0;
			global.particle_dy_so_far[k] = 0.0;
			global.particle_all_dx[k] = 0.0;
			global.particle_all_dy[k] = 0.0;
			global.particle_all_dr2[i] = 0.0;
			global.particle_fx[k] = 0.0;
			global.particle_fy[k] = 0.0;

			// set the color and direction of the particle
			if (j % 2 == 0)
				global.particle_color[k] = 2 + i % 3;
			else
				global.particle_color[k] = 2 + (i + 2) % 3;

			//driection based on color
			switch (global.particle_color[k] - 2)
			{
				case 0:
					global.particle_direction_x[k] = -0.5;
					global.particle_direction_y[k] = -sqrt(3) / 2.0;
					break;
				case 1:
					global.particle_direction_x[k] = 1.0;
					global.particle_direction_y[k] = 0.0;
					break;
				case 2:
					global.particle_direction_x[k] = -0.5;
					global.particle_direction_y[k] = +sqrt(3) / 2.0;
					break;
			}
			k++;
		}

	global.N_particles = k;
	print_log(stdout, NOTE, __FILE__, __LINE__,
			  "Triangular lattice of particles initialized", 3,
			  "N_particles = %d\n", k, "Rescaled SY to = %.2lf\n", global.SY,
			  "Rows = %d\nColumns = %d\nLattice constant = %.2lf\n", N_rows,
			  N_columns, lattice_const);
}

void open_files()
{
	// open file
	global.moviefile = fopen(global.moviefile_name, "wb");
	if (global.moviefile == NULL)
	{
		// exit if does not managed to open
		print_log(stdout, ERROR, __FILE__, __LINE__,
				  "Could not create/open movie file", 0);
		exit(2);
	}

	// open file
	global.statisticsfile = fopen(global.statisticsfile_name, "wt");
	if (global.statisticsfile == NULL)
	{
		// exit if does not managed to open
		print_log(stdout, ERROR, __FILE__, __LINE__,
				  "Could not create/open statistics file", 0);
		exit(2);
	}
	fprintf(global.statisticsfile,
			"time dx dy avg_particle_per_horizontal_pinningsite avg_particle_per_left_down_pinningsite avg_particle_per_left_up_pinningsite avg_particles_per_pinningsite\n");
}

char *substr(const char *from, int start, int count)
{
	char *result = (char *) malloc(count + 1);
	int i, n;
	for (i = start, n = 0; i < start + count && from[i] != '\0'; i++, n++)
	{
		result[n] = from[i];
	}
	result[n] = '\0';
	return result;
}

// int print_log(FILE * stream, LogTypes type, const char *filename,
//            const int line, const char *title, const size_t format_number,
//            ...)
// {
//  int done = 1;
//  va_list args;
//  size_t i;
//  time_t now = time(NULL);
//  char buff[20];
//  // get current time
//  strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
//  COLOR_LOG;
//  printf("LOG[%s]: ", buff);
//  COLOR_DEFAULT;

//  switch (type)
//  {
//      case ERROR:
//          COLOR_ERROR;
//          printf("ERROR ");
//          break;

//      case WARNING:
//          COLOR_WARNING;
//          printf("WARNING ");
//          break;

//      case NOTE:
//          COLOR_NOTE;
//          printf("NOTE ");
//          break;

//      default:
//          COLOR_WARNING;
//          printf("WARNING: LogType not defined!");
//          break;
//  }
//  // print place of occurance
//  printf("(%s: line %d) ", filename, line);
//  // print title
//  if (title)
//      printf("%s", title);
//  printf("\n");
//  COLOR_DEFAULT;

//  // initializing args with the parameters after format_number
//  va_start(args, format_number);
//  for (i = 0; i < format_number; i++)
//  {
//      // reading one element for args as char*
//      char *format = va_arg(args, char *);
//      // print the current note
//      done &= vfprintf(stream, format, args);
//  }
//  va_end(args);

//  return done;
// }

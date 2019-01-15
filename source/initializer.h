//
//  initializer.h
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef initializer_h
#define initializer_h

#include <stdio.h>
#include <stdlib.h>

int read_init_file(const char*);
int check_file_name(const char*, char*);
int get_token(int, int, const char*);
void init_data();

void init_simulation(void);

void init_particles(void);
void init_particles_square_lattice(void);
void init_particles_triangular_lattice(void);
void init_particles_randomly(void);

void init_pinningsites(void);
void init_pinningsites_square_lattice(void);
void init_pinningsites_triangular_lattice(void);

void init_files(void);

double distance_folded_PBC(double x0,double y0,double x1,double y1);

char* substr(const char* from, int start, int count);

#endif /* initializer_h */

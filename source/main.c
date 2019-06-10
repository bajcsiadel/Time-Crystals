//
//  main.c
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//
// first problem: time crystal simulation

#include <stdio.h>
#include "initializer.h"
#include "running.h"
#include "utils.h"

int main(int argc, const char *argv[])
{
	if (argc >= 2)
	{
		if (argc == 3)
			set_seed(argv[2]);
		else
			set_seed(NULL);

		read_init_file(argv[1]);
		init_data();

		printf("Soft Matter Simulator\n");

		init_simulation();
		init_pinningsites();
		init_particles();
		open_files();
		hide_cursor();
		run_simulation();
		show_cursor();
	}
	else
	{
		print_log(stdout, ERROR, __FILE__, __LINE__, "Not enough parameters!",
				  2, "Parameter file missing!", "Usage: %s param.json\n",
				  argv[0]);
	}
	return 0;
}

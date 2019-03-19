//
//  main.c
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//
// first problem: tim,e crystal simulation

#include <stdio.h>
#include "initializer.h"
#include "running.h"
#include "color.h"

int main(int argc, const char * argv[])
{    
    if (argc >= 2)
    {
        if (argc == 3) set_seed(argv[2]);
        else set_seed(NULL);
        
        read_init_file(argv[1]);
        init_data();

        printf("Soft Matter Simulator\n");

        init_simulation();
        init_pinningsites();
        init_particles();
        init_files();
        
        run_simulation();
    } else {
        COLOR_ERROR;
        printf("Not enough parameters! Parameter file missing!\nUsage: %s param.json\n", argv[0]);
        COLOR_DEFAULT;
    }
    return 0;
}

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

int main(int argc, const char * argv[]) {
    
    if (argc == 2)
    {
        read_init_file(argv[1]);
        init_data();

        printf("Soft Matter Simulator\n");

        init_simulation();
        init_pinningsites();
        init_particles();
        init_files();
        
        run_simulation();
    } else {
        printf("\033[1;31mNot enough parameters! Parameter file missing!\nUsage: %s param.json\n\033[0m", argv[0]);
    }
    return 0;
}

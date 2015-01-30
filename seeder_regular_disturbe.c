#include <time.h>

#include "bluebottle.h"
#include "domain.h"
#include "particle.h"

void seeder_a1(int Nx, int Ny, int Nz, real a, real rho, real E, real sigma, int o, int t, int r) {
  printf("Running bluebottle seeder for %d particles...\n\n", Nx*Ny*Nz);
  fflush(stdout);
  int fail = 0;
  
  // read domain input
  domain_read_input();
  domain_init();

  nparts = Nx*Ny*Nz;
	
  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);


  real dx = Dom.xl/Nx; //dx in the distance between centers of two nearby particles in x direction
  real dy = Dom.yl/Ny;
  real dz = Dom.zl/Nz;
  
  
  if(dx < 2*a){
	printf(" Too many particles in x direction\n");
	fail = !fail;
  }
  if(dy < 2*a){
	printf("Too many particles in y direction\n");
	fail = !fail; 
  }
  if(dz < 2*a){
	printf("Too many particles in z direction\n");
	fail = !fail; 
 } 

 if(fail) {
    printf("...bluebottle seeder done.\n\n");
    exit(EXIT_FAILURE);
  }
 // real bias_x;
 // real bias_y;
 // real_bias_z; 
  real bias_x = dx - a;
  real bias_y = dy - a;
  real bias_z = dz - a;

  
  
  for(int k = 0; k < Nz; k++)
  {
  	for (int j = 0; j < Ny; j++)
  	{
  		for (int i = 0; i < Nx; i++)
  		{		
  			//srand(time(NULL));
			parts[i + j*Nx + k*(Nx*Ny)].x = Dom.xs + a + i*dx + bias_x*((double) rand() /(double) (RAND_MAX));
  			parts[i + j*Nx + k*(Nx*Ny)].y = Dom.ys + a + j*dy + bias_y*((double) rand() /(double) (RAND_MAX));
  			parts[i + j*Nx + k*(Nx*Ny)].z = Dom.zs + a + k*dz + bias_z*((double) rand() /(double) (RAND_MAX));
  			parts[i + j*Nx + k*(Nx*Ny)].r = a;
  			parts[i + j*Nx + k*(Nx*Ny)].u = 0;
     			parts[i + j*Nx + k*(Nx*Ny)].v = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].w = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].aFx = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].aFy = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].aFz = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].aLx = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].aLy = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].aLz = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].rho = rho;
      			parts[i + j*Nx + k*(Nx*Ny)].E = E;
      			parts[i + j*Nx + k*(Nx*Ny)].sigma = sigma;
      			parts[i + j*Nx + k*(Nx*Ny)].order = o;
      			parts[i + j*Nx + k*(Nx*Ny)].ncoeff = 0;
      			parts[i + j*Nx + k*(Nx*Ny)].translating = t;
      			parts[i + j*Nx + k*(Nx*Ny)].rotating = r;
  		}
  	}
  }
  	
  printf("Writing part_seeder.input...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE];
  // open file for writing
  sprintf(fname, "%spart_seeder.input", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }
  	
  	
  			
    // write the number of particles
  fprintf(ofile, "n %d\n", nparts);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy, parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy, parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "spring_k %f\n", 0.);//parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n", 0., 0., 0.);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean();			
}  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
 

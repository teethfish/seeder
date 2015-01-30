#include <time.h>

#include "bluebottle.h"
#include "domain.h"
#include "particle.h"

void seeder_hex(int Nx, int Ny, int Nz, real a, real rho, real E, real sigma, int o, int t, int r) {
//Nx is the number of particles in the x-direction in the first layer
//Ny is the number of particles in the y-direction in the first layer
//Nz is the number of layer you would like to have in z direction
//ddz is the distance between the two particle centers in z direction, tan(theta) = dz/(dx/2), shows the hex ratio
  printf("Running bluebottle seeder for hex-array particles...\n\n");
  fflush(stdout);
 // real xx, yy, zz;
 // int fits = 1;
 // int attempts = 1;
  int fail = 0;
  int redo = 1;
  int index = 0; 
  real ddz = 0.87; //ddz should be given in input
  
  // read domain input
  domain_read_input();
  domain_init();

  nparts = (int)((Nz/2)*(Nx*Ny + (Nx-1)*(Ny-1))) + (Nz%2)*Nx*Ny; //calculate the total number of particles
  printf("The total number of particle is %d \n\n",nparts);
  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);


  real dx = Dom.xl/Nx; //dx in the distance between centers of two nearby particles in x direction
  real dy = Dom.yl/Ny;
  redo = 1;
  
  
  
  if(Nz*ddz + 2*a > Dom.zl){
  printf("Too many layers in z direction");
  exit(EXIT_FAILURE);
  }
  if(Nx*2*a > Dom.xl){
  printf("Too many layers in z direction");
  exit(EXIT_FAILURE);
  }
  if(Ny*2*a > Dom.yl){
  printf("Too many layers in z direction");
  exit(EXIT_FAILURE);
  }  
  
  
  int nx, ny;
  int point;
    
  for(int k = 0; k < Nz; k++)
  {	
  	point = k%2;
  	if (point == 1){
  		nx = Nx-1;
  		ny = Ny-1;
  	}
  	else{
  		nx = Nx;
  		ny = Ny;
  	}  	
  	for (int j = 0; j < ny; j++)
  	{
  		for (int i = 0; i < nx; i++)
  		{		
  			parts[index].x = Dom.xs + dx/2 + i*dx + point*dx/2;
  			parts[index].y = Dom.ys + dy/2 + j*dy + point*dy/2;
  			parts[index].z = Dom.zs + a + k*ddz; 				
  			parts[index].r = a;
  			parts[index].u = 0;
     		parts[index].v = 0;
      		parts[index].w = 0;
      		parts[index].aFx = 0;
      		parts[index].aFy = 0;
      		parts[index].aFz = 0;
      		parts[index].aLx = 0;
      		parts[index].aLy = 0;
      		parts[index].aLz = 0;
      		parts[index].rho = rho;
      		parts[index].E = E;
      		parts[index].sigma = sigma;
      		parts[index].order = o;
      		parts[index].ncoeff = 0;
      		parts[index].translating = t;
      		parts[index].rotating = r;
		index = index + 1;
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
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
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
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
 

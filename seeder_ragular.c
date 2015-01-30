#include <time.h>
#include <math.h>

#include "bluebottle.h"
#include "domain.h"
#include "particle.h"

void seeder_regular(int Nx, int Ny, int Nz, real a, real rho, real E, real sigma, int o, int t, int r) {
  printf("Running bluebottle seeder for %d particles...\n\n", Nx*Ny*Nz);
  fflush(stdout);
  int fail = 0;
  // This function is used to give a random field. Each particle can move freely and it will check the distance between another particle
  // bias is the ratio of how  much it can move 
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
  
  // Set the initial regular domain
  for(int k = 0; k < Nz; k++)
  {
  	for (int j = 0; j < Ny; j++)
  	{
  		for (int i = 0; i < Nx; i++)
  		{
  			parts[i + j*Nx + k*(Nx*Ny)].x = Dom.xs + (2*i+1)*dx/2; 
  			parts[i + j*Nx + k*(Nx*Ny)].y = Dom.ys + (2*j+1)*dy/2;
  			parts[i + j*Nx + k*(Nx*Ny)].z = Dom.zs + (2*k+1)*dz/2;	
		}
	 }
   }
  
  real bias = 0.6;
  real x_new = 0.0;
  real y_new = 0.0;
  real z_new = 0.0;
  real d_min = 100*a;
  real d_pair= 0.0;
  
  int times = 10000000;
  for (int t = 0; t < times; t++)
  {
  for(int k = 0; k < Nz; k++)
  {
  	for (int j = 0; j < Ny; j++)
  	{
  		for (int i = 0; i < Nx; i++)
		{	
			d_min = 100*a;
			//printf("i,j,k = %d %d %d \n",i,j,k);	
			//tmp_x =  bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)));
			//printf("tmp_x = %f \n",tmp_x);
			x_new = parts[i + j*Nx + k*(Nx*Ny)].x +  bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)));
  			y_new = parts[i + j*Nx + k*(Nx*Ny)].y + bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)));
  			z_new = parts[i + j*Nx + k*(Nx*Ny)].z + bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)));
			//printf("x_new,y_new,z_new = %f %f %f \n",x_new,y_new,z_new);
			//printf("Dom.xs,Dom.ys,Dom.zs,Dom.xe,Dom.ye,Dom.ze = %f %f %f %f %f %f",Dom.xs,Dom.ys,Dom.zs,Dom.xe,Dom.ye,Dom.ze);
  			if (x_new > Dom.xs  && x_new < Dom.xe && y_new > Dom.ys && y_new < Dom.ye && z_new > Dom.zs && z_new < Dom.ze)
  			{	
  				for (int n = 0; n < Nz; n++)
  				{
  					for (int m = 0; m < Ny; m++)
  					{
  						for (int l = 0; l < Nx; l++)
  						{
  							//printf("l,m,n = %d %d %d \n",l,m,n);
							if(i == l && j == m && k==n)
								d_pair = 100*a;
							else
							{
							d_pair = (x_new - parts[l + m*Nx + n*(Nx*Ny)].x)*(x_new - parts[l + m*Nx + n*(Nx*Ny)].x) 
  							+ (y_new - parts[l + m*Nx + n*(Nx*Ny)].y)*(y_new - parts[l + m*Nx + n*(Nx*Ny)].y)
  							+ (z_new - parts[l + m*Nx + n*(Nx*Ny)].z)*(z_new - parts[l + m*Nx + n*(Nx*Ny)].z);
  							d_pair = sqrt(d_pair);
							}	
							//printf("d_pair = %f\n",d_pair);	
  							if (d_pair < d_min)
  							{	
								d_min = d_pair;
								//printf("d_min = %f\n",d_min);
							}	
						
  						}
  					}
  				}
  				if (d_min > 2*a)
  				{
  					parts[i + j*Nx + k*(Nx*Ny)].x = x_new;
  					parts[i + j*Nx + k*(Nx*Ny)].y = y_new;
  					parts[i + j*Nx + k*(Nx*Ny)].z = z_new;
					//printf("x of particle = %f\n",parts[i + j*Nx + k*(Nx*Ny)].x);
					//printf("y of particle = %f\n",parts[i + j*Nx + k*(Nx*Ny)].y);
					//printf("z of particle = %f\n",parts[i + j*Nx + k*(Nx*Ny)].z);
					
  				}	
  			}
  		}
  	 }
  }
 printf("time = %d\n",t);
  }
	
    for(int k = 0; k < Nz; k++)
  {
        for (int j = 0; j < Ny; j++)
        {
                for (int i = 0; i < Nx; i++)
                {
                        //srand(time(NULL));
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
 // fprintf(ofile, "n %d\n", nparts);
//	fprintf(ofile, "%f %f %f \n", Dom.xl,Dom.yl,Dom.zl);
  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
   // fprintf(ofile, "\n");
   // fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "%f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
   // fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
   // fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy, parts[i].aFz);
   // fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy, parts[i].aLz);
   // fprintf(ofile, "rho %f\n", parts[i].rho);
   // fprintf(ofile, "E %f\n", parts[i].E);
   // fprintf(ofile, "sigma %f\n", parts[i].sigma);
   // fprintf(ofile, "order %d\n", parts[i].order);
   // fprintf(ofile, "spring_k %f\n", 0.);//parts[i].spring_k);
   // fprintf(ofile, "spring (x, y, z) %f %f %f\n", 0., 0., 0.);
   // fprintf(ofile, "translating %d\n", parts[i].translating);
   // fprintf(ofile, "rotating %d\n", parts[i].rotating);
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
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
 

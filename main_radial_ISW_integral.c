/****************************************************************************************************
                       HEADERS
****************************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/****************************************************************************************************
                       INCLUDING SUPPORT FILES
****************************************************************************************************/

#include "variables.c"
#include "reading.c"
#include "functions.c"


/*******************************************************************
NAME: main
FUNCTION: main block
INPUT: parameters file
RETURN: none
******************************************************************/

int main(int argc, char *argv[])
{

  /****** Declaring variables ******/
  char *infile=NULL;


  /****** Reading parameters ******/
  if(argc < 2)
    {
      printf("Error: Incomplete number of parameters. Execute as follows:\n");
      printf("%s Parameters_file\n", argv[0]);
      exit(0);      
    }//if
  
  infile = argv[1];

  
  printf("Reading parameter's file\n");
  printf("--------------------------------------------------------------------------------------\n");
  read_parameters( infile );

  /*+++++ Other variables +++++*/
  GV.ZERO         = 1e-30;      
  
  /*+++++ Reading datafile +++++*/
  printf("Reading the binary file...\n");
  printf("-----------------------------------------\n");
  read_binary();

  GV.CellSize = GV.BoxSize/(1.0*GV.NCELLS);
  GV.c_SL     = 299792.458; // km/s
  GV.CMB_T0   = 2725480.0; // micro K
  GV.CellStep = 1.0*GV.CellSize / 2.0;

  printf("NCELLS=%d, CellSize=%lf, CellsStep=%lf\n", GV.NCELLS, GV.CellSize, GV.CellStep);
  printf("--------------------------------------------------------------------------------------\n");

  /*+++++ Generating rays +++++*/
  ray = (struct radial_rays *) calloc((size_t) GV.NRays, sizeof(struct radial_rays));
  
  rand_rays_coordinates();
   
  /*+++++ Finding intersections and integrating +++++*/
  intersect_integ();
  printf("--------------------------------------------------------------------------------------\n");

  /*+++++ Saving data file +++++*/
  printf("Saving data file");
  printf("--------------------------------------------------------------------------------------\n");
  write_binary();
  
  
  printf("Code finished successfully\n");
  printf("See ya, Mr. Barbarian\n");

  return 0;
}//main

/****************************************************************************************************
                       HEADERS
****************************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_float.h>


/*************************************************************************************
                           DEFINITION OF GLOBAL VARIABLES
*************************************************************************************/
extern double *z_depth=NULL, *PotDot=NULL;


/*************************************************************************************
                       INCLUDING SUPPORT FILES
*************************************************************************************/
#include "variables.c"
#include "reading.c"
#include "interp_integ_gsl.c"
//#include "interp_PotDot_of_Z.c"



/*******************************************************************
NAME: main
FUNCTION: executes the read_data(), interpolation() and SWintegral() functions
INPUT: data file
RETURN: none
******************************************************************/

int main(int argc, char *argv[])
{
  int i, j, k, n, m, count = 0;
  double z, *dT_dr=NULL; 
  char *infile=NULL;
  FILE *pf=NULL;
  char buff[1000];
  double down_lim, up_lim, aux_z, aux_dT;
  char buffer[50];
  char Proj_plane[2];
  double SW_Temp;



  if(argc < 2)
    {
      printf("Error: Incomplete number of parameters. Execute as follows:\n");
      printf("%s Parameters_file\n", argv[0]);
      exit(0);      
    }//if
    
  infile = argv[1];
  
  /*+++++ Reading parameters +++++*/
  printf("Reading parameters file\n");
  printf("-----------------------------------------\n");
  read_parameters( infile );

  /*+++++ Other variables +++++*/
  GV.ZERO         = 1e-30;
  GV.NTOTALCELLS  = GV.NCELLS*GV.NCELLS*GV.NCELLS;
  
  /*+++++ Memory allocation +++++*/
  gp     = (struct grid *) malloc((size_t) GV.NTOTALCELLS*sizeof(struct grid));
  printf("Memory allocated!\n");
  printf("--------------------------------------------------\n");
  

  /*+++++ Reading datafile +++++*/
  printf("Reading the binary file...\n");
  printf("-----------------------------------------\n");
  read_binary();


  GV.CellSize = GV.BoxSize/(1.0*GV.NCELLS);
  GV.c_SL     = 299792.458; // km/s
  GV.CMB_T0   = 2725480.0; // micro K
  GV.CellStep = 1.0*GV.CellSize / 2.0;

  printf("NCELLS=%d, CellSize=%lf, CellsStep=%lf\n", GV.NCELLS, GV.CellSize, GV.CellStep);
  
  /*----- Integration limits -----*/
  down_lim = 0.0;
  up_lim = GV.BoxSize;
  printf("Down_lim=%lf, Up_lim=%lf\n", down_lim, up_lim);
  
  printf("NCells=%d\n", GV.NCELLS);
  printf("--------------------------------------------------\n");
  

  printf("File read!\n");
  printf("--------------------------------------------------\n");
      

  //---------------------------------------------------------
  /*Interpolation of values from exact PotDot*/
  //---------------------------------------------------------  
  printf("Allocating memory\n");
  z_depth = (double *) calloc(GV.NCELLS, sizeof(double));
  PotDot  = (double *) calloc(GV.NCELLS, sizeof(double));
  printf("Memory allocated\n");

#ifdef XLOS  
  pf = fopen( "./../../Processed_data/XLOS/SWIntegral_Exact_XLOS.dat", "w" );  
  fprintf(pf, "#n\t x\t y\t SW_Integral\n");
#endif //XLOS

#ifdef YLOS  
  pf = fopen( "./../../Processed_data/YLOS/SWIntegral_Exact_YLOS.dat", "w" );  
  fprintf(pf, "#n\t x\t y\t SW_Integral\n");
#endif //YLOS

#ifdef ZLOS  
  pf = fopen( "./../../Processed_data/ZLOS/SWIntegral_Exact_ZLOS.dat", "w" );  
  fprintf(pf, "#n\t x\t y\t SW_Integral\n");
#endif //ZLOS


  printf("Beginning code\n");

  
  /*+++++ With GSL methods +++++*/  
  printf("Beginning interpolation and integration with gsl methods\n");
        
  /*----- Creating array with positions. As the box is cubic, 
    it does not matter to put x, y or z positions -----*/
  
   for(i=0; i<GV.NCELLS; i++)     
     z_depth[i] = i*GV.CellSize;
   
   /*::::: X-axis as LOS :::::*/
#ifdef XLOS     
   for(j=0; j<GV.NCELLS; j++)
     {                                                                                                          
       for(k=0; k<GV.NCELLS; k++)                                                                               
	 {                                                                                                      
	   n = INDEX_2D_XLOS(j,k); 
	   fill_potdot_yz(j, k); // this one builds pot_dot(z)	 
	   SW_Temp = GV.a_SF*interp_integ_potdot_dx(down_lim, up_lim);
	   
	   fprintf( pf, 
		    "%12d %16.8f %16.8f %16.8f\n",
		    n, j*GV.CellSize, k*GV.CellSize, SW_Temp);
	 }//for k 
     }//for j
 #endif //XLOS	


   /*::::: Y-axis as LOS :::::*/
 #ifdef YLOS      
   for(i=0; i<GV.NCELLS; i++)
     {                                                                                                          
       for(k=0; k<GV.NCELLS; k++)
	 {                                                                                                      
	   n = INDEX_2D_YLOS(i,k);                                                                                 
	   fill_potdot_xz(i, k); // this one builds pot_dot(z)	 
	   SW_Temp = GV.a_SF*interp_integ_potdot_dx(down_lim, up_lim);
	   
	   fprintf( pf, 
		    "%12d %16.8f %16.8f %16.8f\n",
		    n, i*GV.CellSize, k*GV.CellSize, SW_Temp);
	 }//for j 
     }//for i
 #endif //YLOS



   /*::::: Z-axis as LOS :::::*/
 #ifdef ZLOS      
   for(i=0; i<GV.NCELLS; i++)
     {
       for(j=0; j<GV.NCELLS; j++)
	 { 
	   n = INDEX_2D_ZLOS(i,j);
	   fill_potdot_xy(i, j); // this one builds pot_dot(z)
	   SW_Temp = GV.a_SF*interp_integ_potdot_dx(down_lim, up_lim);
	   
	   /*
	   if((n==0) || (n==256) || (n==1000) )
	     {
	       printf("***********************\n");
	       printf("n=%12d SW_temp=%16.8lf\n", n, SW_Temp);
	       
	       for(k=0; k<GV.NCELLS; k++)
		 {
		   if(k==0 || k==512 || k==1024)
		     printf("z=%lf, PotDot_interp=%lf\n", z_depth[k], PotDot[k]);		  
		 }//for k
	       printf("***********************\n");
	       
	     }//if

	   
	   if(n%100000000==0)
	     {
	       printf("Ready for n=%d with PotDot=%lf\n",
		      n, SW_Temp);
	     }//if 
	   */
	   
	   fprintf( pf, 
		    "%12d %16.8f %16.8f %16.8f\n",
		    n, i*GV.CellSize, j*GV.CellSize, SW_Temp);
	 }//for j 
     }//for i
 #endif //ZLOS			   
	     
  fclose(pf);
  
  printf("Interpolation finished\n");
  printf("-----------------------------------------\n");
  

  /*+++++ dT/dr +++++*/
#ifdef DTDR
  
  printf("Performing dT/dr\n");
  printf("-----------------------------------------\n");
  
  for(i=0; i<GV.NCELLS; i++)
    z_depth[i] = i*GV.CellSize;
  
  /*----- X-axis as LOS -----*/
 #ifdef XLOS

  printf("For XLOS\n");
  pf = fopen("./../../Processed_data/dTdr_XLOS.bin", "w" );
  
  Proj_plane[0] = 'Y';
  Proj_plane[1] = 'Z';

  /*
  snprintf(Proj_plane, sizeof(char)*2, "YZ");  
  printf("Saving dT_dr for plane %s", Proj_plane);
  */
  fwrite(&(GV.BoxSize),      sizeof(double),  1, pf);  // Box Size  
  fwrite((&GV.Omega_M0),     sizeof(double),  1, pf);  // Matter density parameter    
  fwrite((&GV.Omega_L0),     sizeof(double),  1, pf);  // Cosmological constant density parameter  
  fwrite((&GV.z_RS),         sizeof(double),  1, pf);  // Redshift                                   
  fwrite(&(GV.H0),           sizeof(double),  1, pf);  // Hubble parameter                                
  fwrite(&(GV.NCELLS),       sizeof(int),     1, pf);  // Number of cells in one axis
  fwrite(&(Proj_plane[0]),   sizeof(char),    1, pf);  //Projection plane
  fwrite(&(Proj_plane[1]),   sizeof(char),    1, pf);  //Projection plane

  for( j=0; j<GV.NCELLS; j++ )
    {
      for( k=0; k<GV.NCELLS; k++ )
	{	  
	  n = INDEX_2D_XLOS(j,k); 
	  
	  fill_potdot_yz(j, k);	
	  dT_dr = dT_dr_gsl_yz(j, k);
	  
	  fwrite(&(n), sizeof(int), 1, pf);
	  
	  for( i=0; i<GV.NCELLS; i++ )
	    {	
	      m = INDEX_C_ORDER(i,j,k);
	      //m = INDEX_C_ORDER(j,k,i);
	      //fprintf(pf, "%16.8lf %16.8lf\n", i*GV.CellSize, dT_dr[i]);	      
	      //fwrite(&(m),        sizeof(int),    1, pf);	      
	      fwrite(&(i),        sizeof(int),    1, pf); 
	      
	      /* The -1.0 factor is due to the fact that, although I'm integrating from 
		 0 to BoxSize, I'm reducing the distance from BoxSize to 0, then the dr must 
		 be multiplied by -1.0 */
	      dT_dr[i] *= -1.0;
	      fwrite(&(dT_dr[i]), sizeof(double), 1, pf); 
	    }//for k
	  	  
	  //fclose(pf);
	  
	  count += 1;
	  if(count%(4096*4096)==0)
	    {
	      printf("ready for count=%d of %d\n", count, GV.NCELLS*GV.NCELLS);
	    }
	  	  
	}//for j
    }//for i

  fclose(pf);
 #endif //XLOS

  
  /*----- Y-axis as LOS -----*/
 #ifdef YLOS
  
  pf = fopen("./../../Processed_data/dTdr_YLOS.bin", "w" );
  
  Proj_plane[0] = 'X';
  Proj_plane[1] = 'Z';

  /*
  snprintf(Proj_plane, sizeof(char)*2, "XZ");  
  printf("Saving dT_dr for plane %s", Proj_plane);
  */
  fwrite(&(GV.BoxSize),      sizeof(double),  1, pf);  // Box Size                                                     
  fwrite((&GV.Omega_M0),     sizeof(double),  1, pf);  // Matter density parameter                                     
  fwrite((&GV.Omega_L0),     sizeof(double),  1, pf);  // Cosmological constant density parameter                      
  fwrite((&GV.z_RS),         sizeof(double),  1, pf);  // Redshift                                                     
  fwrite(&(GV.H0),           sizeof(double),  1, pf);  // Hubble parameter                                             
  fwrite(&(GV.NCELLS),       sizeof(int),     1, pf);  // Number of cells in one axis
  fwrite(&(Proj_plane[0]),   sizeof(char),    1, pf);  //Projection plane
  fwrite(&(Proj_plane[1]),   sizeof(char),    1, pf);  //Projection plane
  
  
  for( i=0; i<GV.NCELLS; i++ )
    {
      for( k=0; k<GV.NCELLS; k++ )
	{
	  /*
	    if( (k==0 && i==0) 
	    || (k==(GV.NCELLS/8) && i==(GV.NCELLS/8)) 
	    || (k==(GV.NCELLS/4) && i==(GV.NCELLS/4)) 
	    || (k==(GV.NCELLS/2) && i==(GV.NCELLS/2)) 
	    || (k==(GV.NCELLS-1) && i==(GV.NCELLS-1)) )
	    {
	  */  
	  n = INDEX_2D_YLOS(i,k); 
	  
	  fill_potdot_xz(i, k);	
	  dT_dr = dT_dr_gsl_xz(i, k);
	  
	  fwrite(&(n),        sizeof(int),    1, pf);
	  
	  /*
	    if( (k==0 && i==0) 
	    || (k==(GV.NCELLS/8) && i==(GV.NCELLS/8)) 
	    || (k==(GV.NCELLS/4) && i==(GV.NCELLS/4)) 
	    || (k==(GV.NCELLS/2) && i==(GV.NCELLS/2)) 
	    || (k==(GV.NCELLS-1) && i==(GV.NCELLS-1)) )
	    {	      
	  */
	  //snprintf(buffer, sizeof(char)*50, "./../../Processed_data/YLOS/dT_dr_i%d_k%d.txt", i, k);
	  //pf = fopen(buffer, "w");
	      
	  for( j=0; j<GV.NCELLS; j++ )
	    {		  
	      //fprintf(pf, "%16.8lf %16.8lf\n", j*GV.CellSize, dT_dr[j]);
	      m = INDEX_C_ORDER(i,j,k);
	      //m = INDEX_C_ORDER(i,k,j);
	      //fprintf(pf, "%16.8lf %16.8lf\n", i*GV.CellSize, dT_dr[i]);
	      
	      fwrite(&(j),        sizeof(int),    1, pf);	      
	      /* The -1.0 factor is due to the fact that, although I'm integrating from 
		 0 to BoxSize, I'm reducing the distance from BoxSize to 0, then the dr must 
		 be multiplied by -1.0 */
	      dT_dr[j] *= -1.0;
	      fwrite(&(dT_dr[j]), sizeof(double), 1, pf);	      
	    }//for j
	  
	  count += 1;
	  if(count%(4096*4096)==0)
	    {
	      printf("ready for count=%d of %d\n", count, GV.NCELLS*GV.NCELLS);
	    }

    
	  //fclose(pf);
	  
	  //}//if
	  
	}//for k
    }//for i

  fclose(pf);

 #endif //YLOS
  

  /*----- Z-axis as LOS -----*/
  #ifdef ZLOS

  pf = fopen("./../../Processed_data/dTdr_ZLOS.bin", "w" );

  Proj_plane[0] = 'X';
  Proj_plane[1] = 'Y';
  
  /*
  snprintf(Proj_plane, sizeof(char)*2, "XY");  
  printf("Saving dT_dr for plane %s", Proj_plane);
  */

  fwrite(&(GV.BoxSize),      sizeof(double),  1, pf);  // Box Size                                                     
  fwrite((&GV.Omega_M0),     sizeof(double),  1, pf);  // Matter density parameter                                     
  fwrite((&GV.Omega_L0),     sizeof(double),  1, pf);  // Cosmological constant density parameter                      
  fwrite((&GV.z_RS),         sizeof(double),  1, pf);  // Redshift                                                     
  fwrite(&(GV.H0),           sizeof(double),  1, pf);  // Hubble parameter                                             
  fwrite(&(GV.NCELLS),       sizeof(int),     1, pf);  // Number of cells in one axis
  fwrite(&(Proj_plane[0]),   sizeof(char),    1, pf);  //Projection plane
  fwrite(&(Proj_plane[1]),   sizeof(char),    1, pf);  //Projection plane
  
  
  for( i=0; i<GV.NCELLS; i++ )
    {
      for( j=0; j<GV.NCELLS; j++ )
	{
	  /*
	  if( (j==0 && i==0) 
	      || (j==(GV.NCELLS/8) && i==(GV.NCELLS/8)) 
	      || (j==(GV.NCELLS/4) && i==(GV.NCELLS/4)) 
	      || (j==(GV.NCELLS/2) && i==(GV.NCELLS/2)) 
	      || (j==(GV.NCELLS-1) && i==(GV.NCELLS-1)) )
	    {	      
	  */  
	      n = INDEX_2D_ZLOS(i,j); 
	      
	      fill_potdot_xy(i, j);	
	      dT_dr = dT_dr_gsl_xy(i, j);

	      fwrite(&(n),        sizeof(int),    1, pf);
	      
	      /*
		if( (j==0 && i==0) 
		|| (j==(GV.NCELLS/8) && i==(GV.NCELLS/8)) 
		|| (j==(GV.NCELLS/4) && i==(GV.NCELLS/4)) 
		|| (j==(GV.NCELLS/2) && i==(GV.NCELLS/2)) 
		|| (j==(GV.NCELLS-1) && i==(GV.NCELLS-1)) )
		{	      
	      */
	      //snprintf(buffer, sizeof(char)*50, "./../../Processed_data/ZLOS/dT_dr_i%d_j%d.txt", i, j);
	      //pf = fopen(buffer, "w");
	      
	      for( k=0; k<GV.NCELLS; k++ )
		{		  
		  m = INDEX_C_ORDER(i,j,k);		  
		  //fprintf(pf, "%16.8lf %16.8lf\n", i*GV.CellSize, dT_dr[i]);
		  
		  fwrite(&(k),        sizeof(int),    1, pf);	      
		  dT_dr[k] *= -1.0;
		  fwrite(&(dT_dr[k]), sizeof(double), 1, pf);	      
		  //fprintf(pf, "%16.8lf %16.8lf\n", k*GV.CellSize, dT_dr[k]);
		}//for k
	     
	      count += 1;
	      if(count%(4096*4096)==0)
		{
		  printf("ready for count=%d of %d\n", count, GV.NCELLS*GV.NCELLS);
		}
	      
 
	      //fclose(pf);
	      
	      //}//if
	      
	}//for j
    }//for i
  fclose(pf);

 #endif //ZLOS

#endif //DTDR  

  free(z_depth);
  free(PotDot);
    
  printf("Code finished!\n");  
  printf("-----------------------------------------\n");

}//main

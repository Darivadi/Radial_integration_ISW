/****************************************************************************************************
NAME: rand_rays_coordinates
FUNCTION: Generates random spherical coordinates for each ray
INPUT: 
RETURN: 0
****************************************************************************************************/
int rand_rays_coordinates(void)
{
  int m;
  const gsl_rng_type * T; /*Define el tipo de generador de números aleatorios. No hay que liberarlo*/
  gsl_rng * r; /*Análogo al w. Puntero que contiene la info sobre cual generador se va a usar,cantidad de memoria a usar, etc.*/
  long seed;

  /*+++++ Initializing random generation of numbers +++++*/
  gsl_rng_env_setup();//Inicializa las rutinas de generación                                                    
  T = gsl_rng_default;/*Inicialización de T con esta variable de GSL que es la default*/
  r = gsl_rng_alloc(T);/*Alocación de memoria*/
  
  seed = time(NULL)*getpid();
  
  gsl_rng_set(r, seed);/*Recibe puntero de inicialización de generación y  un entero largo como semilla*/
  
  for(m=0; m<GV.NRays; m++)
    {     
      
      ray[m].rad   = 0.5*GV.BoxSize;
      ray[m].theta = acos(1.0 - 2.0 * gsl_rng_uniform (r)); //Polar
      ray[m].phi   = 2.0*M_PI * gsl_rng_uniform (r); //Azimuth
      
      
      /*----- Initial position's vector -----*/
      ray[m].vec_ini[X] = 0.0;
      ray[m].vec_ini[Y] = 0.0;
      ray[m].vec_ini[Z] = 0.0;
      
      /*----- Final position's vector -----*/
      ray[m].vec_end[X] = ray[m].rad * sin(ray[m].theta) * cos(ray[m].phi);
      ray[m].vec_end[Y] = ray[m].rad * sin(ray[m].theta) * sin(ray[m].phi);
      ray[m].vec_end[Z] = ray[m].rad * cos(ray[m].theta);
      
      if(m<2 || m==GV.NRays-1)
	{
	  printf("Ray m=%d generated at (theta, phi) = (%10.5lf, %10.5lf)\n", 
		 m, ray[m].theta, ray[m].phi);
	  printf("This corresponds to (x,y,z) = (%10.5lf, %10.5lf, %10.5lf)\n", 
		 ray[m].vec_end[X], ray[m].vec_end[Y], ray[m].vec_end[Z]);
	}//if
      
    }//for m
  
  gsl_rng_free (r);

  printf("Rays generated\n");
  printf("--------------------------------------------------------------------------------------\n");
 
  return 0;
}//rand_rays_coordinates


/****************************************************************************************************
NAME: intersect_integ(void)
FUNCTION: Finds the intersections of the rays with the grid and performs the integration as 
INPUT: 
RETURN: 0
****************************************************************************************************/
int intersect_integ(void)
{
  int i, j, k, n, m, p;
  double x0, y0, z0, xf, yf, zf, rad_max;
  double tMax_x, tMax_y, tMax_z, tMin_all; 
  double* dist_trav=NULL, *PotDot=NULL;

  printf("Let's find intersections");
  printf("--------------------------------------------------------------------------------------\n");

  dist_trav = (double *) calloc((size_t) (GV.HalfNCells), sizeof(double) );
  PotDot    = (double *) calloc((size_t) (GV.HalfNCells), sizeof(double) );

  printf("Memory allocated for dist_trav and PotDot arrays\n");
  printf("--------------------------------------------------------------------------------------\n");
  
  
  for(m=0; m<GV.NRays; m++)
    {      
      rad_max = 0.0;
      x0 = y0 = z0 = xf = yf = zf = 0.0;      
              
      for(p=0; p<GV.HalfNCells; p++)
	{
	  dist_trav[p] = 0.0;
	  PotDot[p]    = 0.0;
	}//for p
      
      p = 1;
      
      do 
	{
	  x0 = xf;
	  y0 = yf;
	  z0 = zf;	 

	  if(p>=GV.HalfNCells)
	    {
	      printf("Oops! p=%d >= HalfNCells=%d\n", p, GV.HalfNCells);
	      return 1;
	    }//if

	  /*+++++ Computing parameter t in equation \vec[u] + t * \vec[v] = vec_end +++++*/
	  if( (ray[m].vec_end[X] - ray[m].vec_ini[X]) > 0.0 || (ray[m].vec_end[X] - ray[m].vec_ini[X]) < 0.0 )
	    tMax_x = fabs( (p*GV.CellSize - ray[m].vec_ini[X]) / (ray[m].vec_end[X] - ray[m].vec_ini[X]) );
	  else
	    tMax_x = ray[m].rad;
	  
	  if( (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) > 0.0 || (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) < 0.0 )
	    tMax_y = fabs( (p*GV.CellSize - ray[m].vec_ini[Y]) / (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) );
	  else
	    tMax_y = ray[m].rad;
	  
	  if( (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) > 0.0 || (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) < 0.0 )
	    tMax_z = fabs( (p*GV.CellSize - ray[m].vec_ini[Z]) / (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) );
	  else
	    tMax_z = ray[m].rad;

	  /*
	  if( m<2 && p<2 )
	    {
	      printf("-------------------------------------------------------------------------------\n");
	      printf("tx, ty, tz = %10.5lf, %10.5lf, %10.5lf\n", tMax_x, tMax_y, tMax_z);
	    }//if
	  */
	  /*+++++ Finding the minimum of the 3 t parameters +++++*/
	  if (tMax_x < tMax_y)
	    {
	      tMin_all = tMax_x;
	      
	      if (tMax_z < tMin_all)
		{
		  tMin_all = tMax_z;	      
		}//if 1
	    }//if
	  else
	    {
	      tMin_all = tMax_y;
	      
	      if (tMax_z < tMin_all)
		{
		  tMin_all = tMax_z;
		}//if 1
	    }//else
	  /*
	  if(m<2 && p<2)
	    {	      
	      printf("For m=%d, p=%d, tMin_all=%lf\n", m, p, tMin_all);
	    }//if
	  */
	  xf = ray[m].vec_ini[X] + tMin_all * (ray[m].vec_end[X] - ray[m].vec_ini[X]);
	  yf = ray[m].vec_ini[Y] + tMin_all * (ray[m].vec_end[Y] - ray[m].vec_ini[Y]);
	  zf = ray[m].vec_ini[Z] + tMin_all * (ray[m].vec_end[Z] - ray[m].vec_ini[Z]);
	  	  
	  /*+++++ Computing index of the corresponding cell +++++*/
	  i = floor( ( (xf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  j = floor( ( (yf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  k = floor( ( (zf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  n = INDEX_C_ORDER(i,j,k);	  
	  
	  /*+++++ Assigning distance and PotDot to compute Integral  +++++*/
	  dist_trav[p-1] = sqrt( POW2(xf - x0) + POW2(yf - y0) + POW2(zf - z0) );
	  rad_max = sqrt( POW2(xf) + POW2(yf) + POW2(zf) );
	  PotDot[p-1] = gp[n].potDot_r;
	  /*
	  if(m<2 && p<2)
	    {	      
	      printf("Intersection at (x,y,z)=(%10.5lf,%10.5lf,%10.5lf) at cell n=%d\n", xf, yf, zf, n);
	      printf("Distance traveled: %10.5lf and PotDot = %10.5lf\n", dist_trav[p-1], PotDot[p-1]);
	    }//if	  	  
	  */	  	  

	  p++;	  
	} while( rad_max <= (ray[m].rad - 1e-2) );
      /*
      if(m<2)
	{
	  printf("For the last cell:\n");
	  printf("rad_max=%10.5lf\n", rad_max);
	  printf("Intersection at (x,y,z)=(%10.5lf,%10.5lf,%10.5lf) at cell n=%d\n", xf, yf, zf, n);
	  printf("Distance traveled: %10.5lf and PotDot = %10.5lf\n", dist_trav[p-1], PotDot[p-1]);
	}//if
      */
      /*+++++ Computing integral +++++*/
      /*
      if(m<2)
	printf("Computing integral for m=%d\n", m);
      */
      ray[m].ISW_temp = 0.0;
	
      for(p=0; p<(GV.NCELLS/2); p++)
	{
	  ray[m].ISW_temp += dist_trav[p] * PotDot[p] * GV.a_SF;
	}//for p

      if(m<2)
	printf("For ray m=%d, ISW temperature is %10.5lf\n", m, ( 2.0*GV.CMB_T0/(POW3(GV.c_SL)) ) * ray[m].ISW_temp);
      
      if(m%100==0)
	printf("Ready for m=%d\n", m);      
      
    }//for m

  free(dist_trav);
  free(PotDot);
  
  return 0;
}//intersect_integ

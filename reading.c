/****************************************************************************************************
NAME: conf2dump
FUNCTION: Reads the input file with parameters
INPUT: Parameters file
RETURN: 0
****************************************************************************************************/

int conf2dump( char filename[] )
{
  int nread;
  char cmd[1000];
  /*
  sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | gawk -F\"=\" '{print $2}' > %s.dump", 
	   filename, filename );
  */
  sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | awk -F\"=\" '{print $2}' > %s.dump", 
	   filename, filename );
  nread = system( cmd );
  return 0;
}


/****************************************************************************************************
NAME: read_parameters
FUNCTION: Reads the parameters
INPUT: Parameters file
RETURN: 0
****************************************************************************************************/
int read_parameters( char filename[] )
{
  int nread;
  char cmd[1000], filenamedump[1000];
  FILE *file;
  
  /*+++++ Loading the file +++++*/
  file = fopen( filename, "r" );
  if( file==NULL )
    {
      printf( "  * The file '%s' doesn't exist!\n", filename );
      return 1;
    }
  fclose(file);
  
  /*+++++ Converting to plain text +++++*/
  conf2dump( filename );
  sprintf( filenamedump, "%s.dump", filename );
  file = fopen( filenamedump, "r" );
  
  /*+++++ Parameters for binary data +++++*/
  nread = fscanf(file, "%d", &GV.NCELLS);
  nread = fscanf(file, "%s", GV.FILENAME);
  nread = fscanf(file, "%d", &GV.NRays);

  printf("Number of cells %10d\n", GV.NCELLS);
  printf("Number of rays %10d\n", GV.NRays);  
  printf("Data file at %s\n", GV.FILENAME);
  
  fclose( file );
  
  printf( "  * The file '%s' has been loaded!\n", filename );
  
  sprintf( cmd, "rm -rf %s.dump", filename );
  nread = system( cmd );
  
  return 0;
}



/**************************************************************************************************** 
NAME: read_binary
FUNCTION: Reads the binary data file
INPUT: None
RETURN: 0 
****************************************************************************************************/
int read_binary(void)
{
  int i, nread;
  double pos_aux[3], dummy;
  FILE *inFile=NULL;
  
  inFile = fopen(GV.FILENAME, "r");

  printf("Reading simulation parameters\n");

  
  /*+++++ Saving Simulation parameters +++++*/
  nread = fread(&GV.BoxSize,  sizeof(double), 1, inFile);  //Box Size
  nread = fread(&GV.Omega_M0, sizeof(double), 1, inFile);  //Matter density parameter
  nread = fread(&GV.Omega_L0, sizeof(double), 1, inFile);  //Cosmological constant density parameter
  nread = fread(&GV.z_RS,     sizeof(double), 1, inFile);  //Redshift
  nread = fread(&GV.H0,       sizeof(double), 1, inFile);  //Hubble parameter
  nread = fread(&GV.NCELLS,   sizeof(int),    1, inFile);  //Number or cells

  GV.a_SF = 1.0 / (1.0 + GV.z_RS);
  GV.NTOTALCELLS = GV.NCELLS * GV.NCELLS * GV.NCELLS;

  /*+++++ Memory allocation +++++*/
  gp     = (struct grid *) calloc((size_t) GV.NTOTALCELLS, sizeof(struct grid));
  printf("Memory allocated!\n");
  printf("--------------------------------------------------\n");
  

  
  printf("-----------------------------------------------\n");
  printf("Cosmological parameters:\n");
  printf("OmegaM0=%lf OmegaL0=%lf redshift=%lf HubbleParam=%lf\n",
	 GV.Omega_M0,
	 GV.Omega_L0,
	 GV.z_RS,
	 GV.H0);
  printf("-----------------------------------------------\n");
  
  printf("Simulation parameters:\n");
  printf("L=%lf\n",
	 GV.BoxSize);
  printf("-----------------------------------------------\n");

      
  for(i=0; i<GV.NTOTALCELLS; i++)
    {
      /*..... File app2 .....*/
        nread = fread(&gp[i].potDot_r, sizeof(double), 1, inFile);
	
	if(i%5000000==0)
	  {
	    printf("Ready for i=%d with PotDot=%lf\n", 
		   i, gp[i].potDot_r);
	  }//if
    }//for m	      

  printf("Data read!\n");
  fclose(inFile);
  
  return 0;
}//read_binary


/**************************************************************************************************** 
NAME: write_binary
FUNCTION: writes file in binary format
INPUT: None
RETURN: 0 
****************************************************************************************************/
int write_binary(void)
{
  FILE *outFile=NULL;
  int m, i, j, k;
  
  pf = fopen("./../../Processed_data/ISW_radial_app2.bin", "w");
  
  fwrite(&GV.BoxSize,  sizeof(double), 1, outFile);  // Box Size                                                     
  fwrite(&GV.Omega_M0, sizeof(double), 1, outFile);  // Matter density parameter                                     
  fwrite(&GV.Omega_L0, sizeof(double), 1, outFile);  // Cosmological constant density parameter                      
  fwrite(&GV.z_RS,     sizeof(double), 1, outFile);  // Redshift                                                     
  fwrite(&GV.H0,       sizeof(double), 1, outFile);  // Hubble parameter                                             
  fwrite(&GV.NCELLS,   sizeof(int),    1, outFile);  // Number of cells per axis
  fwrite(&GV.NRays,    sizeof(int),    1, outFile);  // Number of radial rays
  
  for(m=0; m<GV.NRays; m++)
    {
      fwrite(&(ray[m].rad),      sizeof(double), 1, outFile);
      fwrite(&(ray[m].theta),    sizeof(double), 1, outFile);
      fwrite(&(ray[m].phi),      sizeof(double), 1, outFile);
      fwrite(&(ray[m].ISW_temp), sizeof(double), 1, outFile);
    }//for m                                                                                                    

  fclose(outFile);


}//write_binary

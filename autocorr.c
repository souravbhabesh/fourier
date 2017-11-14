#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "fftw3.h"
#define M_PI 3.14159265358979323846

#define MAXLENGTH 10000000
#define PERIOD 10000
#define RUNMAX 20
#define NMAX 20000
#define MAXFRAMES 20000

double cnode[RUNMAX][MAXFRAMES];
double tseries[MAXFRAMES];
double autocorr[RUNMAX][MAXFRAMES];
double jk_blocks[NMAX][RUNMAX];
double jk_error[NMAX];
double error[NMAX]; //RMSE error in |hFT|^2

int NX,NY,RUNS,steps,frames,JK_BIN_COUNT;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin,*file;
   char data_file[1024],line[256];
   int i,j,n,signal_length;
   double *hd,*corr_h;
   fftw_complex *hFT;
   fftw_plan pdir,pinv;
   char const* fileName;

   switch (argc){
   case 6:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       steps = atoi(argv[4]);
       fileName = argv[5];
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa steps runlist.dat\n");
   }

   if(NULL==(file=fopen(fileName,"r")))  
	print_and_exit("I could not open file with simulation run numbers %s\n",fileName);
 
   frames = steps/PERIOD;
   n = frames/2; // Number of real data points for which FFT is taken
   signal_length = n;
   //JK_BIN_COUNT = RUNS; //Number of Jack Knife bins

   //for(int run=0;run<RUNS;run++)
 
   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/cnode.bin",NX,NY,KAPPA);
   //printf("%s\n",data_file);
   if(NULL==(Fin=fopen(data_file,"rb")))
   	print_and_exit("I could not open binary file %s\n",data_file);

   //printf("frames %d\n",frames);

   //We read the data file
   fread(cnode,sizeof(double),RUNMAX*MAXFRAMES,Fin);
   fclose(Fin);

   for(int i=0;i<frames;i++)
   {
	//printf("%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",i,cnode[0][i],cnode[1][i],cnode[2][i],cnode[3][i],cnode[4][i],cnode[5][i]);
	//tseries[i]=cnode[0][i];
   }

   int runnum,run_cnt=0;
   while (fscanf(file, "%d", &runnum) == 1)// 1 is returned if fscanf reads a number
   {
	   //printf("Run # %d \n",runnum);
	   hd = (double *) fftw_malloc(sizeof(double)*n*2);
	   corr_h = (double *) fftw_malloc(sizeof(double)*n*2);	
	   hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));

	   // We prepare the FFTW
	   pdir = fftw_plan_dft_r2c_1d(2*n,hd,hFT,FFTW_MEASURE);
	   pinv = fftw_plan_dft_c2r_1d(2*n,hFT,corr_h,FFTW_MEASURE);

	   // Now fill in the vector hd
	   for(i=0; i<n; i++)
	     hd[i] = cnode[runnum-1][i];
	   for(i=n; i<2*n; i++)
	     hd[i]=0;

	   // We execute the FFTW
	   fftw_execute(pdir); 

	   // hFT contains the FT of h, we need to compute | hFT | ^ 2
           for (i=0;i<n+1;i++){
      		hFT[i][0] = hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1];
      		hFT[i][1] = 0;
    		}
    	   fftw_execute(pinv);

	   // The inverse FT of | hFT | ^ 2 is autocorrelation, but we must
    	   //normalize
    	   for(i=0; i<n; i++)
	   {
      		corr_h[i] /= (n-i);
		autocorr[runnum-1][i] = corr_h[i]/corr_h[0];
	   }

	   //for (i=0;i<n;i++)
      		//printf("%d %.8g\n", i, (hFT[i][0])/(pow(n,2)));

	  //for (i=0;i<n;i++)
      		//printf("%d %.8g\n", i, corr_h[i]/corr_h[0]); 
	  run_cnt++;
   }
/*
   for (i=0;i<n;i++)
   	printf("%d %.8g %.8g %.8g %.8g %.8g\n", i, autocorr[0][i],autocorr[1][i],autocorr[2][i],autocorr[3][i],autocorr[4][i]);
*/
   fftw_destroy_plan(pdir);
   fftw_destroy_plan(pinv);
   fftw_free(hd);
   fftw_free(hFT);
   fftw_free(corr_h);

    JK_BIN_COUNT = run_cnt; //Number of Jack Knife bins
	   
    //		Jack Knife Error estimation	
    
    //		Total of Jack Knife blocks at each sampling interval	
    for(i=0;i<n;i++)
    {
	jk_blocks[i][JK_BIN_COUNT]=0;
	for(j=0;j<JK_BIN_COUNT;j++)
	{
		jk_blocks[i][JK_BIN_COUNT] += autocorr[j][i]; //summing IFT value at each i for different runs
	}
    }    

    //		Jack Knife Blocking	
    for(i=0;i<n;i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        { 
		jk_blocks[i][j]= (jk_blocks[i][JK_BIN_COUNT] - autocorr[j][i])/(JK_BIN_COUNT-1);	
	}
    }

    double jk_error_term1[NMAX],jk_error_term2[NMAX];
    //		Jack Knife Error	
    for(i=0;i<n;i++)
    {
	jk_error[i]=0;
	jk_error_term1[i]=0;
	jk_error_term2[i]=0;
    } 

    for(i=0;i<n;i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term1[i] += jk_blocks[i][j] * jk_blocks[i][j];
	}
	jk_error_term1[i] = (1.0/JK_BIN_COUNT) * jk_error_term1[i];
    }

    for(i=0;i<n;i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term2[i] += (1.0/JK_BIN_COUNT) * jk_blocks[i][j];
	}
	jk_error_term2[i] = jk_error_term2[i] * jk_error_term2[i];
        //	JK Error	
	jk_error[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1[i] - jk_error_term2[i]));
	printf ("%d\t%.8g\t%.8g\t%.8g\n",i,i*(2*M_PI/n),jk_blocks[i][JK_BIN_COUNT]/JK_BIN_COUNT,jk_error[i]);	
    }

    return 0;   
}

void print_and_exit(char *format, ...)
{
    va_list list;

    va_start(list,format);
    vprintf(format,list);
    va_end(list);
    exit(1);
}


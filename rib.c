#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "fftw3.h"

#define MAXLENGTH 10000000
#define PERIOD 10000
#define RUNMAX 20
#define NMAX 2000
#define MAXFRAMES 11000

double v[MAXLENGTH]; // The list of measures
double h_width[MAXFRAMES][NMAX];
double power_spectrum_frame[MAXFRAMES][NMAX];
double avg_hFT[RUNMAX][NMAX];
double avg_h[RUNMAX][NMAX];
double hFT_width_avg[NMAX];//Averaging fourier amplitude square across runs
double jk_blocks[NMAX][RUNMAX];
double jk_error[NMAX];
double error[NMAX]; //RMSE error in |hFT|^2

int NX,NY,RUNS,steps,frames,JK_BIN_COUNT;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin;
   char data_file[1024];
   int i,j,n,signal_length;
   double *hd;
   fftw_complex *hFT;
   fftw_plan pdir;

   switch (argc){
   case 6:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       RUNS = atoi(argv[4]);
       steps = atoi(argv[5]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa RUNS steps  \n");
   }

   frames = steps/PERIOD;
   n = 2*NX; // Number of real data points for which FFT is taken
   signal_length = n;
   JK_BIN_COUNT = RUNS; //Number of Jack Knife bins

   for(int run=0;run<RUNS;run++)
   {
	   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/width.bin",NX,NY,KAPPA,run+1);
	   //printf("%s\n",data_file);
	   if(NULL==(Fin=fopen(data_file,"rb")))
	      print_and_exit("I could not open binary file %s\n",data_file);

	    //We read the data file
	    fread(h_width,sizeof(double),frames*NX,Fin);
	    fclose(Fin); 

	    /*	Allocating memory for FFT	*/
	    hd = (double *) fftw_malloc(sizeof(double)*n);
	    hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*((n/2)+1));

	    // Plan for FFTW	
	    pdir = fftw_plan_dft_r2c_1d(n,hd,hFT,FFTW_PATIENT);

	    // Now fill in the vector hd; initializing of input array should be done after creating the plan
	    for(j=0;j<frames/2;j++)
	    {
		for(i=0; i<n; i++)
		{
			hd[i] = h_width[j][i];
			//if(run==0 && j==0)
				//printf("%d\t%.8f\n",i,h_width[j][i]);
			//printf("%.8g\n",hd[i]);
		}
	   	/*	Padding with 0's to make it periodic for Correlation evaluation	*/ 
	   	//for(i=n; i<2*n; i++)
			//hd[i]=0;	

	   	//Execute the FFTW
	   	fftw_execute(pdir);

	   	// hFT contains the FT of hd, we need to compute | hFT | ^ 2
	   	for (i=0;i<((n/2)+1);i++)
	   	{
			//hFT[i][0] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(signal_length/2,2));
			//hFT[i][1] = 0;
			power_spectrum_frame[j][i] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(n/2,2));
			//if(run==0 && j==0)
				//printf("%d\t%.8f\t%.8f\t%.8f\n",i,hFT[i][0]/(n),hFT[i][1]/(n),power_spectrum_frame[j][i]);
	   	}
	     }

	     //Adding fourier amplitudes at same x for different frames same run
	     for(i=0;i<((n/2)+1);i++)
           	{
                	avg_hFT[run][i]=0;
                	for(j=0;j<frames/2;j++)
                	{
                        	avg_hFT[run][i]+=power_spectrum_frame[j][i];
                	}
                	avg_hFT[run][i]/=(frames/2);
                	//printf("%d\t%.8f\n",run,avg_hFT[run][i]);
           	} 
    }
	   

    /* Average of power spectrum across runs	*/
    for(i=0;i<((n/2)+1);i++)
    {
	hFT_width_avg[i]=0;
	//printf("%d",i);
	for(int r=0;r<RUNS;r++)
	{
		hFT_width_avg[i]+=avg_hFT[r][i];
	}
	hFT_width_avg[i]/=RUNS;
	//printf("%.8f",hFT_width_avg[i]);
    }

   /*      RMSE error in power spectrum   */
    for(i=0;i<((n/2)+1);i++)
    {
	//printf("%d\t%.8g\t",i,hFT_width_avg[i]/hFT_width_avg[0]);
	//if (i==1)
		//printf("%.8f\t%.8f\n",sqrt(3)*(NY-1)/(2*NX-1),hFT_width_avg[i]/hFT_width_avg[i+1]);
	error[i]=0;
	for(int r=0;r<RUNS;r++)
	{
		error[i] += pow(((avg_hFT[r][i]/avg_hFT[r][0])-(hFT_width_avg[i]/hFT_width_avg[0])),2);
	}
	error[i] = sqrt(error[i]/RUNS);
	//printf("%.8g\n",error[i]);
    }  

    /*		Jack Knife Error estimation	*/
    
    /*		Total of Jack Knife blocks at each sampling interval	*/
    for(i=0;i<((n/2)+1);i++)
    {
	jk_blocks[i][JK_BIN_COUNT]=0;
	for(j=0;j<JK_BIN_COUNT;j++)
	{
		jk_blocks[i][JK_BIN_COUNT] += avg_hFT[j][i]/avg_hFT[j][0]; //summing avg_hFT for all runs at each i
	}
    }    

    /*		Jack Knife Blocking	*/
    for(i=0;i<((n/2)+1);i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        { 
		jk_blocks[i][j]= (jk_blocks[i][JK_BIN_COUNT]-(avg_hFT[j][i]/avg_hFT[j][0]))/(JK_BIN_COUNT-1);	
	}
    }

    double jk_error_term1[NMAX],jk_error_term2[NMAX];
    /*		Jack Knife Error	*/
    for(i=0;i<((n/2)+1);i++)
    {
	jk_error[i]=0;
	jk_error_term1[i]=0;
	jk_error_term2[i]=0;
    } 

    for(i=0;i<((n/2)+1);i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term1[i] += jk_blocks[i][j] * jk_blocks[i][j];
	}
	jk_error_term1[i] = (1.0/JK_BIN_COUNT) * jk_error_term1[i];
    }

    for(i=0;i<((n/2)+1);i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term2[i] += (1.0/JK_BIN_COUNT) * jk_blocks[i][j];
	}
	jk_error_term2[i] = jk_error_term2[i] * jk_error_term2[i];
        /*	JK Error	*/
	jk_error[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1[i] - jk_error_term2[i]));
	printf ("%d\t%.8g\t%.8g\n",i,jk_blocks[i][JK_BIN_COUNT]/JK_BIN_COUNT,jk_error[i]);	
    }

    //printf("Cleaning up\n");
    fftw_destroy_plan(pdir);
    fftw_free(hd);
    fftw_free(hFT);
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


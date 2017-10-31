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
double h_width[MAXFRAMES/2][2*NMAX];
double power_spectrum_frame[MAXFRAMES/2][2*NMAX];
double avg_hFT[RUNMAX][2*NMAX];
double avg_h[RUNMAX][2*NMAX];
double hFT_width_avg[2*NMAX];//Averaging fourier amplitude square across runs
double error[2*NMAX]; //RMSE error in |hFT|^2

int NX,NY,RUNS,steps,frames;
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
   n = 2*NX;
   signal_length = 2*n;//Taking into account padding with 0's

   for(int run=0;run<RUNS;run++)
   {
	   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/width.bin",NX,NY,KAPPA,run+1);
	   //printf("%s\n",data_file);
	   if(NULL==(Fin=fopen(data_file,"rb")))
	      print_and_exit("I could not open binary file %s\n",data_file);

	    //We read the data file
	    fread(h_width,sizeof(double),frames*NX,Fin);
	    fclose(Fin); 

	    hd = (double *) fftw_malloc(sizeof(double)*2*n);
	    hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));

	    // Plan for FFTW	
	    pdir = fftw_plan_dft_r2c_1d(2*n,hd,hFT,FFTW_PATIENT);

	    // Now fill in the vector hd; initializing of input array should be done after creating the plan
	    for(j=0;j<frames/2;j++)
	    {
		for(i=0; i<n; i++)
		{
			hd[i] =h_width[j][i];
			//if(run==0 && j==0)
				//printf("%d\t%.8f\n",i,h_width[j][i]);
			//printf("%.8g\n",hd[i]);
		}
		      
	   	for(i=n; i<2*n; i++)//Padding with 0's to make it periodic
			hd[i]=0;	

	   	//Execute the FFTW
	   	fftw_execute(pdir);

	   	// hFT contains the FT of h, we need to compute | hFT | ^ 2
	   	for (i=0;i<n+1;i++)
	   		{
				//hFT[i][0] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(signal_length/2,2));
				//hFT[i][1] = 0;
				power_spectrum_frame[j][i] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(n,2));
				//if(run==0 && j==0)
					//printf("%d\t%.8f\t%.8f\t%.8f\n",i,hFT[i][0]/(n),hFT[i][1]/(n),power_spectrum_frame[j][i]);
	   		}
	     }

	     //Adding fourier amplitudes at same x for different frames same run
	     for(i=0;i<n+1;i++)
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
    for(i=0;i<n+1;i++)
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
            for(i=0;i<n+1;i++)
            {
                printf("%d\t%.8g\t",i,hFT_width_avg[i]);
                //if (i==1)
                	//printf("%.8f\t%.8f\n",sqrt(3)*(NY-1)/(2*NX-1),hFT_width_avg[i]/hFT_width_avg[i+1]);
                error[i]=0;
                for(int r=0;r<RUNS;r++)
                {
                        error[i] += pow(avg_hFT[r][i]-hFT_width_avg[i],2);
                }
                error[i] = sqrt(error[i]/RUNS);
                printf("%.8g\n",error[i]);
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


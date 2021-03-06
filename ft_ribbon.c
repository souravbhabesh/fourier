#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "fftw3.h"

#define MAXLENGTH 10000000
//#define STEPS 100000000
#define PERIOD 10000
//#define FRAMES 	STEPS/PERIOD  //Total number of frames in the last half of the simulation
//#define NX 101
//#define NY 51
//#define KAPPA 5.0
#define RUNS 2
#define NMAX 2000
#define MAXFRAMES 11000

double v[MAXLENGTH]; // The list of measures
double h_width[MAXFRAMES/2][2*NMAX];
double width_hFT[MAXFRAMES/2][2*NMAX];
double avg_hFT[RUNS][2*NMAX];
double avg_h[RUNS][2*NMAX];
double hFT_width_avg[2*NMAX];//Averaging fourier amplitude square across runs
double error[2*NMAX]; //RMSE error in |hFT|^2

int NX,NY,steps,frames;
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
   case 5:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       steps = atoi(argv[4]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa steps  \n");
   }

   frames = steps/PERIOD;

   for(int run=0;run<RUNS;run++)
   {
	   sprintf(data_file,"width_L%d_W%d_k%.1f_r%d.bin",NX,NY,KAPPA,run+1);
	   //printf("%s\n",data_file);
	   if(NULL==(Fin=fopen(data_file,"rb")))
	      print_and_exit("I could not open binary file %s\n",data_file);

	    //We read the data file
	    fread(h_width,sizeof(double),frames*NX,Fin);
	    fclose(Fin); 

		/*    for(int i=0;i<2*NX;i++)
		    {
			printf("%d",i);
		for(int j=0;j<FRAMES;j++)
		{
			printf("\t%f",h_width[i][j]);
		}
		printf("\n");
	     }
	*/

	    n = 2*NX;
	    signal_length = 2*n;//Taking into account padding with 0's
	    //v_mean /= n; //divide by the total number of measure

	    hd = (double *) fftw_malloc(sizeof(double)*2*n);
	    hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));

	    // Now fill in the vector hd
	    //printf("%d\n",n);
	    for(j=0;j<frames/2;j++)
	    {
		for(i=0; i<n; i++)
		{
			hd[i] = h_width[j][i];//-v_mean;
			if(run==0 && j==0)
				printf("%d\t%.8f\n",i,h_width[j][i]);
		//printf("%.8g\n",hd[i]);
		}
		      
		for(i=n; i<2*n; i++)//Padding with 0's to make it periodic
			hd[i]=0;	

		// Plan for FFTW
		pdir = fftw_plan_dft_r2c_1d(2*n,hd,hFT,FFTW_ESTIMATE);

		//Execute the FFTW
		fftw_execute(pdir);

		// hFT contains the FT of h, we need to compute | hFT | ^ 2
		for (i=0;i<n+1;i++){
			//hFT[i][0] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(signal_length/2,2));
			//hFT[i][1] = 0;
			width_hFT[j][i] = hFT[i][0];
			if(run==0 && j==0)
				printf("%d\t%.8f\t%.8f\t%.8f\n",i,hFT[i][0]/(n),hFT[i][1]/(n),(hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(n,2)));
		}
	    }


	    for(i=0;i<n+1;i++)
	    {
		//printf("%d",i);
		avg_hFT[run][i]=0;
		for(j=0;j<frames/2;j++)
		{
			avg_hFT[run][i]+=width_hFT[j][i];//Adding fourier amplitudes at same x for different frames
		}
		avg_hFT[run][i]/=(frames/2);
		//printf("\t%f",avg_hFT[i]);
		//printf("\n");
	     } 

	    for(i=0;i<n+1;i++)
	    {
		//printf("%d",i);
		avg_h[run][i]=0;
		for(j=0;j<frames/2;j++)
		{
			avg_h[run][i]+=h_width[j][i];
		}
		avg_h[run][i]=avg_h[run][i]/(frames/2);
		//printf("\t%f",avg_h[i]);
		//printf("\n");
	     }
		 
		    //for (i=0;i<n;i++)
	      //printf("%d %.8g\n",i,hFT[i][0]);
		// (1000/signal_length)*double(i),2*hFT[i][0]); 
	}

	/*	Averaging |hFT|^2 across runs	*/
	    for(i=0;i<n+1;i++)
            {
		hFT_width_avg[i]=0;
		//printf("%d",i);
		for(int r=0;r<RUNS;r++)
		{
			hFT_width_avg[i]+=avg_hFT[r][i];
                }
                hFT_width_avg[i]/=RUNS;
                //printf("\t%f",hFT_width_avg[i]);
                //printf("\n");
             }

	/*	RMSE error in |hFT|^2	*/
	    for(i=0;i<n+1;i++)
            {
		//printf("%d\t%.8g\t",i,hFT_width_avg[i]);
		error[i]=0;
		for(int r=0;r<RUNS;r++)
                {
			error[i] += pow(avg_hFT[r][i]-hFT_width_avg[i],2);
		}
		error[i] = sqrt(error[i]/RUNS);
		//printf("%.8g\n",error[i]);
	    }
	   
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


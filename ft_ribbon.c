#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "fftw3.h"

#define MAXLENGTH 10000000
#define STEPS 100000000
#define PERIOD 10000
//#define FRAMES STEPS/(2*PERIOD)
#define FRAMES 100 //Total number of frames in the last half of the simulation
#define NX 101

double v[MAXLENGTH]; // The list of measures
double h_width[2*NX][FRAMES];
double width_hFT[2*NX][FRAMES];
double avg_hFT[NX];


void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin;
   char data_file[1024];
   char line[MAXLENGTH];
   int i,n,strip;
   double read_v,v_mean,read_v2;
   double *hd;
   fftw_complex *hFT;
   fftw_plan pdir;

   switch (argc){
   case 2:
       sscanf(argv[1],"%s",data_file); //Data file
       break;
   default:
       print_and_exit("Usage: %s data file \n",argv[0]);
   }

   if(NULL==(Fin=fopen(data_file,"rt")))
      print_and_exit("I could not open %s\n",data_file);

    //We read the data file
    n=0;
    i=0;
    char *token;
    int signal_length,col_cnt=0;
    while(fgets(line,MAXLENGTH, Fin) != NULL)//Read a line
    {
      i++;//Number of lines being read, Length of the signal
      token = strtok(line, "\t");
      
      //Iterating through the tokenized line items
      while( token != NULL )
      {
	if(col_cnt != 0 && col_cnt != FRAMES+1)
	{
		h_width[i-1][col_cnt-1] = atof(token);
		//printf( "%f %f %d %d\t",atof(token),h_width[i-1][col_cnt-1],i-1,col_cnt-1);
	}
      	//printf( "%s\t", token );
        token = strtok(NULL, "\t");
	col_cnt++;
      }
      //printf("\n");
      col_cnt=0;//reset colum counter in new line

      //sscanf(line,"%d\t%lf\t%lf",&strip,&read_v,&read_v2);
      //v[n]=read_v; //We record the measurements on the vector v
      //v_mean+= read_v; //We accumulate to evaluate the average value
      //printf("%.8g\n",v[n]);
      n++;
    }
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
    signal_length = i;
    //v_mean /= n; //divide by the total number of measure

    hd = (double *) fftw_malloc(sizeof(double)*2*n);
    hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));

    // Now fill in the vector hd
    //printf("%d\n",n);
    for(int j=0;j<FRAMES;j++)
    {
	for(i=0; i<n; i++)
	{
		hd[i] = h_width[i][j];//-v_mean;
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
		hFT[i][0] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(signal_length,2));
		hFT[i][1] = 0;
		width_hFT[i][j] = hFT[i][0];
		//printf("\t%f",width_hFT[i][j]);
    	}
    }

    for(int i=0;i<NX;i++)
    {
        printf("%d",i);
	avg_hFT[i]=0;
        for(int j=0;j<FRAMES;j++)
        {
		avg_hFT[i]+=width_hFT[i][j];//Adding fourier amplitudes at same x for different frames
                //printf("\t%f",width_hFT[i][j]);
        }
        avg_hFT[i]/=FRAMES;
        printf("\t%f",avg_hFT[i]);
        printf("\n");
     }   
 
    //for (i=0;i<n;i++)
      //printf("%d %.8g\n",i,hFT[i][0]);
	// (1000/signal_length)*double(i),2*hFT[i][0]); 
   
    fftw_destroy_plan(pdir);
    fftw_free(hd);
    fftw_free(hFT);
    fclose(Fin);

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


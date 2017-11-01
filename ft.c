#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "fftw3.h"
void print_and_exit(char *, ...); //Print out an error message and exits
#define MAXLENGTH 10000000
double v[MAXLENGTH]; // The list of measures

int main(int argc, char **argv)
{

   FILE *Fin;
   char data_file[1024];
   char line[1024];
   int i,n,strip;
   double read_v,v_mean;
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
    int signal_length;
    while(fgets(line,100, Fin) != NULL)//Read a line
    {
      i++;//Number of lined being read, Length of the signal
      sscanf(line,"%d\t%lf",&strip,&read_v);
      v[n]=read_v; //We record the measurements on the vector v
      v_mean+= read_v; //We accumulate to evaluate the average value
      //printf("%.8g\n",v[n]);
      n++;
    }

    signal_length = n;
    v_mean /= n; //divide by the total number of measure

    hd = (double *) fftw_malloc(sizeof(double)*n);
    hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));

   // Plan for FFTW
    pdir = fftw_plan_dft_r2c_1d(n,hd,hFT,FFTW_ESTIMATE);

    // Now fill in the vector hd
    //printf("%d\n",n);
    for(i=0; i<n; i++)
      {
	hd[i] = v[i];//-v_mean;
	//printf("%.8g\n",hd[i]);
      }
      
    //for(i=n; i<2*n; i++)//Padding with 0's to make it periodic
      //hd[i]=0;	

    //Execute the FFTW
    fftw_execute(pdir);

    // hFT contains the FT of h, we need to compute | hFT | ^ 2
    for (i=0;i<(n/2)+1;i++){
      hFT[i][0] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(signal_length/2,2));
      hFT[i][1] = 0;
    }

    for (i=0;i<n;i++)
      printf("%d %.8g\n",i,hFT[i][0]);
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


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
   char fichero[1024];
   char linea[100];
   int n, ndes,i;
   double lee_v,vmed;
   double *hd,*corr_h;
   fftw_complex *hFT;
   fftw_plan pdir,pinv;

   switch (argc){
     case 3:
       sscanf(argv[2],"%d",&ndes);    //Measures Discarded for Thermalization 
       sscanf(argv[1],"%s",fichero); //Data file
       break;
     default:
       print_and_exit("Usage: %s fichero ndes\n",
	   argv[0]);
   }

    if(NULL==(Fin=fopen(fichero,"rt")))
      print_and_exit("No he podido abrir %s\n",fichero);

    //We read the data file
    n=0;
    i=0;
    vmed=0;
    while(fgets(linea,100, Fin) != NULL)//Read a line
    {
      i++;
      if(i>ndes)//We discard the first measures
      {
	sscanf(linea,"%lf",&lee_v);
	v[n]=lee_v; //We record the measurements on the vector v
	vmed+= lee_v; //We accumulate the average value
	n++;
      }
    }
    vmed /= n; //We divide by the total number of measure

    //We will calculate the autocorrelation of v with an FFT, using
    // the Wiener-Kinchin theorem.
    // First we declare the necessary vectors (with double length,
    // must be filled with zeros to take into account the conditions
    // periodic contour of the FFT)

    hd     = (double *) fftw_malloc(sizeof(double)*n*2);
    corr_h = (double *) fftw_malloc(sizeof(double)*n*2);
    hFT    = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));


    // Now fill in the vector hd
    for(i=0; i<n; i++)
      hd[i] = v[i]-vmed;
    for(i=n; i<2*n; i++)
      hd[i]=0;

    
    // We prepare the FFTW
    pdir = fftw_plan_dft_r2c_1d(2*n,hd,hFT,FFTW_ESTIMATE);
    pinv = fftw_plan_dft_c2r_1d(2*n,hFT,corr_h,FFTW_ESTIMATE);

    
    // We execute the FFTW
    fftw_execute(pdir);

    // hFT contains the FT of h, we need to compute | hFT | ^ 2
    for (i=0;i<n+1;i++){
      hFT[i][0] = hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1];
      hFT[i][1] = 0;
    }
    fftw_execute(pinv);


    // The inverse TF of | hFT | ^ 2 is autocorrelation, but we must
    //normalize
    for(i=0; i<n; i++)
      corr_h[i] /= (n-i);

    //for (i=0;i<n;i++)
      //printf("%d %.8g\n", i, corr_h[i]/corr_h[0]);

    for (i=0;i<n;i++)
      printf("%d %.8g\n", i, (hFT[i][0])/(pow(n,2)));//hFT[i][0] is already the modulus

    fftw_destroy_plan(pdir);
    fftw_destroy_plan(pinv);
    fftw_free(hd);
    fftw_free(hFT);
    fftw_free(corr_h);
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

#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include <cstdlib>
#include <cmath>
#include <limits>


#define MAXFRAMES 11000
#define NMAX 200
#define PERIOD 10000
#define RUN 10


double h_width[MAXFRAMES/2][2*NMAX];

void print_and_exit(char *format, ...)
{
    va_list list;

    va_start(list,format);
    vprintf(format,list);
    va_end(list);
    exit(1);
}

double generateGaussianNoise(double mu, double sigma)
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;

	thread_local double z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

int main(int argc, char **argv)
{
	FILE *Fout;
	char outfile[1024];
	int i,j,r,NX,NY,frames,steps;
	double KAPPA;

	switch (argc){
	   case 5:
	       NX = atoi(argv[1]);
	       NY = atoi(argv[2]);
	       KAPPA = atof(argv[3]);
	       steps = atoi(argv[4]);
	       break;
	default:
	       print_and_exit("Usage Pass command line arguments: NX NY KAPPA steps\n"); //total steps
	   }
	frames = steps/PERIOD;
	//printf("frames: %d\n",frames);       

	for(r=0;r<RUN;r++)
	{

		sprintf(outfile,"width_L%d_W%d_k%.1f_r%d.bin",NX,NY,KAPPA,r+1);
		if(NULL==(Fout=fopen(outfile,"wb")))
		      print_and_exit("I could not open binary file %s\n",outfile);

		for(i=0;i<frames/2;i++)
		{
			for(j=0;j<2*NX;j++)
			{
				h_width[i][j] = 0.7*sin(2*M_PI*10*j/(2*NX))+ sin(2*M_PI*120*j/(2*NX))+generateGaussianNoise(0,.1);
				//generateGaussianNoise(0,1.0);// 0.7*sin(2*M_PI*50*j/(2*NX))+ sin(2*M_PI*120*j/(2*NX));
				if(i==0 && r==0)
					printf("%d\t%.8f\n",j,h_width[i][j]);
				//printf("%d\t%.8f\n",i,generateGaussianNoise(0,1.0));
				//printf("%.8f\n",generateGaussianNoise(0,1.0));
			}
		}

		fwrite(h_width,sizeof(double),frames*NX,Fout);
		fclose(Fout);
	}
	return 0;
}

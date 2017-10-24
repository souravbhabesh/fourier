#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include <cstdlib>
#include <cmath>
#include <limits>

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
	int i,num;
	switch (argc){
	   case 2:
	       num = atoi(argv[1]);
	       break;
	default:
	       print_and_exit("Usage Pass command line arguments:Number of Gaussian RN to be generated  \n");
	   }       
	
	for(i=0;i<num;i++)
	{
		//printf("%d\t%.8f\n",i,generateGaussianNoise(0,0.1));
		printf("%.8f\n",generateGaussianNoise(0,0.1));
	}
	return 0;
}

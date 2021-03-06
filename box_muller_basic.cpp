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
	int i,num,toggle;
	switch (argc){
	   case 3:
	       num = atoi(argv[1]);
	       toggle = atoi(argv[2]);
	       break;
	default:
	       print_and_exit("Usage Pass command line arguments:Number of Gaussian RN to be generated  and toggle value 1 2 3 4 5\n");
	   }       
	
	for(i=0;i<num;i++)
	{
		if(toggle == 1)
			printf("%d\t%.8f\n",i,generateGaussianNoise(0,0.1));
		if(toggle == 2)
			printf("%d\t%.8f\n",i,sin(2*M_PI*10*i/(2*num)));
		if(toggle == 3)
			printf("%d\t%.8f\n",i,0.7*sin(2*M_PI*120*i/(2*num)));
		if(toggle == 4)
			printf("%d\t%.8f\n",i,0.7*sin(2*M_PI*60*i/(2*num))+ sin(2*M_PI*10*i/(2*num)));
		if(toggle == 5)
			printf("%d\t%.8f\n",i,0.7*sin(2*M_PI*60*i/(2*num))+ sin(2*M_PI*10*i/(2*num))+ generateGaussianNoise(0,0.1));

	}
	return 0;
}

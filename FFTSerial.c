//**************************************************************
// Names: Anthony Enem and Ali Khalid
//***************************************************************
// This is a serial implementation of the radix-2 FFT algorithm.
// This implementation first calculates all the euler values (e^(-itheta))
// and stores those in a complex array. Then the euler values are re-used in
// computing the FFT values from k = 0 to N-1. The main loop only goes from
// k = 1 to N/2 - 1 because the even and odd components for the first half
// can be reused to calculate the values for the last half.
//*****************************************************************

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

#define N 16384

typedef struct {
	double real;
	double imag;
} Complex;

//Arrays for Input and Results
Complex Input[N];
Complex Result[N];

//Array for Euler values
Complex Euler[N/2];

//***************************************************************************
//computeEulers()
//Parameters: None
//The function precomputes the euler values from x = 0 to N/2 of e^(-4*PI*x/N) 
//to avoid recomputation of sines and cosines in the FFT
//***************************************************************************
void computeEulers()
{
	int x = 0;
	float theta;
	float n = (4.0*M_PI)/N;

	for(x = 0; x < (N>>1); x++){
		theta = x*n;
		Euler[x].real = cos(theta);
		Euler[x].imag = -sin(theta);
	}
}

//***************************************************************************
//multiply()
//Parameters: a and b. Two pointers to complex structs
//The function multiplies the complexes referenced by a and b and returns the result
//***************************************************************************
Complex multiply(Complex* a, Complex* b)
{
	Complex c;
	c.real = (a->real * b->real) - (a->imag * b->imag);
	c.imag = (a->real * b->imag) + (b->real * a->imag);
	return c;
}

//***************************************************************************
//add()
//Parameters: a and b. Two pointers to complex structs
//The function adds the complexes referenced by a and b and returns the result
//***************************************************************************
Complex add(Complex* a, Complex* b)
{
	Complex c;
	c.real = a->real+b->real;
	c.imag = a->imag+b->imag;
	return c;
}

//Main
int main(int argc, char** argv)
{
	//Redirect output to file
	freopen("EnemKhalidSerial.txt", "w", stdout);

	//Input values
	Input[0].real = 3.6; Input[0].imag = 2.6;
	Input[1].real = 2.9; Input[1].imag = 6.3;
	Input[2].real = 5.6; Input[2].imag = 4.0;
	Input[3].real = 4.8; Input[3].imag = 9.1;
	Input[4].real = 3.3; Input[4].imag = 0.4;
	Input[5].real = 5.9; Input[5].imag = 4.8;
	Input[6].real = 5.0; Input[6].imag = 2.6;
	Input[7].real = 4.3; Input[7].imag = 4.1;

	int n, k;
	for(n = 8; n < N; n++){
		Input[n].real = Input[n].imag = 0.0;
	}

	Complex even, odd;
	struct timespec now, tmstart;
	clock_gettime(CLOCK_REALTIME, &tmstart);

	//compute N/2 Euler values (cos(theta) - isin(theta))
	computeEulers();

	Complex twiddle, temp, euler;

	double theta, PI2_by_N = 2*M_PI/N;

	int diff, idx;

	for(k = 0; k < (N>>1); k++)
	{
		even.real = even.imag = odd.real = odd.imag = 0.0;

		//Get difference between indices of Euler values for current k
		diff = (k - 1 +(N>>1)) % (N>>1);
		idx = 0; //start index is 0

		for(n = 0; n < (N>>1); n++){
			//get current euler component
			euler = Euler[idx];

			//multiply even input with euler component
			temp = multiply(&Input[n<<1], &euler);    
			//add result to even         
			even = add(&even, &temp);

			//multiply odd component with euler input
			temp = multiply(&Input[(n<<1) +1], &euler);  
			//add result to odd
			odd = add(&odd, &temp);

			//compute index for next euler component
			idx = (idx + diff + 1) % (N>>1);
		}

		//Compute twiddle
		theta = k*PI2_by_N;
		twiddle.real = cos(theta);
		twiddle.imag = -sin(theta);
		//multiply twiddle with odd component
		temp = multiply(&odd, &twiddle);

		//Add even and odd to result
		Result[k] = add(&even, &temp);

		//Subtract odd from even to get k+N/2 component
		temp.real = -temp.real;
		temp.imag = -temp.imag;
		Result[k+(N>>1)] = add(&even, &temp);
	}

	clock_gettime(CLOCK_REALTIME, &now);
	double seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
	printf("C time: %f seconds\n", seconds);

	printf("TOTAL PROCESSED SAMPLES: %d\n", N);
	printf("============================================\n");
	for(k = 0; k <= 10; k++){
		printf("XR[%d]: %.5f  \t XI[%d]: %.5f\n", k, Result[k].real, k, Result[k].imag);
		printf("============================================\n");
	}


}

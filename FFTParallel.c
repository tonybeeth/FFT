//**************************************************************
// Name: Anthony Enem and Ali Khalid
//***************************************************************
// This is a parallel implementation of the radix-2 FFT algorithm.
// All processes first calculate (N/2)/comm_size of the euler values
// (e^(-theta)) and then distribute to every other process using 
// MPI_Allgather. Then each process computes (N/2)/comm_size of 
// the result and gathers their results to process 0 which ouputs it.
//
// Compilation: mpicc EnemKhalidParallel.c -lm -o parallel.out
// Execution: qsub script
// **Modification to the script: we added a change directory command 
// so the output is created in the same directory.
//
//*****************************************************************

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include<stdlib.h>

#define N 16384

typedef struct {
	double real;
	double imag;
} Complex;


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
	int size, rank;
	Complex Input[N];

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

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int localN = (N/2)/size;
	double EulerR[N/2], EulerI[N/2];
	double *ResultR = NULL, *ResultI = NULL;
	double tempEulerI[localN], tempEulerR[localN];
	double tempResultR[localN*2], tempResultI[localN*2];	
	
	struct timespec now, tmstart;
	double mpist, mpiend;

	if(rank==0)
	{
		ResultR = malloc(sizeof(double)*N);
		ResultI = malloc(sizeof(double)*N);

		//start timers
		clock_gettime(CLOCK_REALTIME, &tmstart);
		mpist = MPI_Wtime();
	}

	//Compute Euler values
	double theta, ang = 4.0*M_PI/N;
	int x;
	for(x = 0; x < localN; x++){
		theta = (localN*rank + x)*ang;
		tempEulerR[x] = cos(theta);
		tempEulerI[x] = -sin(theta);
	}

	//Distribute euler values to all processes
	MPI_Allgather(tempEulerR, localN, MPI_DOUBLE, EulerR, localN, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(tempEulerI, localN, MPI_DOUBLE, EulerI, localN, MPI_DOUBLE, MPI_COMM_WORLD);

	Complex even, odd;

	Complex twiddle, temp, euler, result;

	double PI2_by_N = 2*M_PI/N;

	int diff, idx;

	for(x = 0; x < localN; x++)
	{
		k = rank*localN + x;
		even.real = even.imag = odd.real = odd.imag = 0.0;

		//Get difference between indices of Euler values for current k
		diff = (k - 1 +(N>>1)) % (N>>1);
		idx = 0; //start index is 0

		for(n = 0; n < (N>>1); n++){
			//get current euler component
			euler.real = EulerR[idx];
			euler.imag = EulerI[idx];

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
		result = add(&even, &temp);
		tempResultR[x] = result.real;
		tempResultI[x] = result.imag;

		//Subtract odd from even to get k+N/2 component
		temp.real = -temp.real;
		temp.imag = -temp.imag;
		result = add(&even, &temp);
		tempResultR[x+localN] = result.real;
		tempResultI[x+localN] = result.imag;
	}

	//Gather results to process 0
	
	//First N/2 results
	MPI_Gather(tempResultR, localN, MPI_DOUBLE, ResultR, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(tempResultI, localN, MPI_DOUBLE, ResultI, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//Second N/2 results
	MPI_Gather(tempResultR+localN, localN, MPI_DOUBLE, ResultR+(N/2), localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(tempResultI+localN, localN, MPI_DOUBLE, ResultI+(N/2), localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank == 0)
	{
		//Redirect output to file
		freopen("EnemKhalidParallel.txt", "w", stdout);
	
		//End timers
		mpiend = MPI_Wtime();
		clock_gettime(CLOCK_REALTIME, &now);
		double seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));

		printf("C time: %f secs\n", seconds);
		printf("MPI time: %f secs\n", mpiend-mpist);

		printf("TOTAL PROCESSED SAMPLES: %d\n", N);
		printf("============================================\n");
		for(k = 0; k <= 10; k++){
			printf("XR[%d]: %.5f  \t XI[%d]: %.5f\n", k, ResultR[k], k, ResultI[k]);
			printf("============================================\n");
		}
	}
	MPI_Finalize();
}

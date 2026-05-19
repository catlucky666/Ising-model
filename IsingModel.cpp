#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define SIZE 10
#define K_B  1
#define J    1 
//-----------------------------------------------------------------------------
// calculation  energy
//-----------------------------------------------------------------------------
int CalculateEnergy(int iSpin[SIZE][SIZE])
{
	int iE_Total, iNeighbors;

	iE_Total = 0;
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			iNeighbors = iSpin[(i - 1 + SIZE) % SIZE][j]
			     	   + iSpin[(i + 1) % SIZE][j]
				       + iSpin[i][(j - 1 + SIZE) % SIZE]
				       + iSpin[i][(j + 1) % SIZE];
			iE_Total = iE_Total - J * iSpin[i][j] * iNeighbors;
		}
	}
	return (iE_Total / 2);
}
//-----------------------------------------------------------------------------
// Print spin configuration
//-----------------------------------------------------------------------------
void PrintSpin(int iSpin[SIZE][SIZE])
{
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			printf("%2d ", iSpin[i][j]);
		}
		printf("\n");
	}
}
//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	int iSpin[SIZE][SIZE];
	int iE_Total, iE_LocalOld, iE_LocalNew, iE_Delta;
	int iStep, iStep_Max, iNeighbors;
	int iX, iY;
	double dTemperature, dProb;
	double dProb_Target[3];
//-----------------------------------------------------------------------------
// Initialization
//-----------------------------------------------------------------------------
	srand((unsigned int)(time(NULL)));
	iStep_Max = 5000;
	dTemperature = 1.0;
//-------------------------------------------------------------------
// Generate spin configuration
//-------------------------------------------------------------------
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			iSpin[i][j] = (rand() % 2) * 2 - 1;
		}
	}
//-------------------------------------------------------------------
// Pre-calculate target probability
//-------------------------------------------------------------------
	dProb_Target[0] = 1.0;
	dProb_Target[1] = exp(-4.0 * J / (K_B * dTemperature));
	dProb_Target[2] = exp(-8.0 * J / (K_B * dTemperature));
//-----------------------------------------------------------------------------
// Calculate total energy
//-----------------------------------------------------------------------------
	PrintSpin(iSpin);
	iE_Total = CalculateEnergy(iSpin);
	printf("Initial energy = %d\n", iE_Total);
//-----------------------------------------------------------------------------
// Do Monte Carlo steps
//-----------------------------------------------------------------------------
	for (iStep = 1; iStep <= iStep_Max; iStep++)
	{
//-------------------------------------------------------------------
// Do move
//-------------------------------------------------------------------
		iX = (int)(SIZE * rand() / (RAND_MAX + 1.0));
		iY = (int)(SIZE * rand() / (RAND_MAX + 1.0));

		iNeighbors = iSpin[(iX - 1 + SIZE) % SIZE][iY]
		           + iSpin[(iX + 1) % SIZE][iY]
			       + iSpin[iX][(iY - 1 + SIZE) % SIZE]
			       + iSpin[iX][(iY + 1) % SIZE];

		iE_Delta = 2 * J * iSpin[iX][iY] * iNeighbors;
//-------------------------------------------------------------------
// Accept or not
//-------------------------------------------------------------------
		if (iE_Delta <= 0)
		{
//---------------------------------------------------------
// Accept the move and update the spin and energy
//---------------------------------------------------------
			iE_Total += iE_Delta;
			iSpin[iX][iY] = -iSpin[iX][iY];
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		}
		else
		{
//---------------------------------------------------------
// Reject or not
//---------------------------------------------------------
			dProb = (double)(rand() / (RAND_MAX + 1.0));

			if (dProb < dProb_Target[iE_Delta / (4 * J)])
			{
//-----------------------------------------------
// Accept the move and update the spin and energy
//-----------------------------------------------
				iE_Total += iE_Delta;
				iSpin[iX][iY] = -iSpin[iX][iY];
//- - - - - - - - - - - - - - - - - - - - - - - -
			}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		}
	}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	printf("\n");
	printf("Final energy = %d\n", iE_Total);
	printf("\n");
//-----------------------------------------------------------------------------
// Output final results
//-----------------------------------------------------------------------------
	PrintSpin(iSpin);
	//iE_Total = CalculateEnergy(iSpin);
	//printf("%d\n", iE_Total);
//-----------------------------------------------------------------------------
// The operation completed successfully
//-----------------------------------------------------------------------------
	return 0;
}

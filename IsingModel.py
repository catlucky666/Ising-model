import random as rd
import numpy as np

class Ising_Model:
    def __init__(self, size=10, J=1, dTemperature=1.1, iStep_Max=5000, seed=None):
# -----------------------------------------------------------------------------
#        Constructor initializes the physical system and immediately
#        performs a full Monte Carlo simulation (as in the original main()).

#        Parameters:
#        - size        : Lattice dimension (size x size)
#        - J           : Coupling constant
#        - dTemperature: Temperature (k_B assumed = 1)
#        - iStep_Max   : Total Monte Carlo steps
#        - seed        : Random seed for reproducibility (optional)
# -----------------------------------------------------------------------------
#******************************************************************************#
# -----------------------------------------------------------------------------
# Store basic simulation parameters
# -----------------------------------------------------------------------------
        self.SIZE = size
        self.J = J
        self.dTemperature = dTemperature
        self.iStep_Max = iStep_Max
        if seed is not None:
            rd.seed(seed)

#-------------------------------------------------------------------
# Initialize spin lattice
#-------------------------------------------------------------------
        # 2D list representing spin configuration (+1 or -1)
        self.iSpin2D = [[0] * self.SIZE for _ in range(self.SIZE)]
#---------------------------------------------------------
# Random initialization: each spin is either +1 or -1
#---------------------------------------------------------
        for i in range(self.SIZE):
            for j in range(self.SIZE):
                r = rd.randrange(2)      # randomly choose 0 or 1
                self.iSpin2D[i][j] = r * 2 - 1  # map 0->-1, 1->+1

#-------------------------------------------------------------------
        # Precompute Boltzmann probabilities
#-------------------------------------------------------------------
        # In 2D Ising with nearest neighbors,
        # possible ΔE values are {0, 4J, 8J} (positive cases only).
        # We precompute exp(-ΔE / T) to avoid recomputing repeatedly.
        self.dProb_Target = [0.0] * 3
        self.dProb_Target[0] = 1.0  # ΔE = 0 → always accept
        self.dProb_Target[1] = np.exp(-4.0 * self.J / self.dTemperature)
        self.dProb_Target[2] = np.exp(-8.0 * self.J / self.dTemperature)

# --------------------------------------------------
# Compute initial energy
# --------------------------------------------------
        self.iE_Total = self.CalculateEnergy()

        # Record energy evolution (time series)
        self.iEnergy_list = []
        # Record magnetization evolution (normalized -1..+1)
        self.dMagnetization_list = []
# --------------------------------------------------
 # Run Metropolis Monte Carlo
# --------------------------------------------------
    def MetropolisStep(self):
        # Randomly select a lattice site
        iX = rd.randrange(self.SIZE)
        iY = rd.randrange(self.SIZE)
        # Compute sum of 4 nearest neighbors (periodic boundary condition)
        iNeighbors = (
            self.iSpin2D[(iX - 1 + self.SIZE) % self.SIZE][iY]
          + self.iSpin2D[(iX + 1) % self.SIZE][iY]
          + self.iSpin2D[iX][(iY - 1 + self.SIZE) % self.SIZE]
          + self.iSpin2D[iX][(iY + 1) % self.SIZE]
         )
         # Energy change if this spin is flipped
         # ΔE = 2 * J * s_i * (sum of neighbors)
        iE_Delta = 2 * self.J * self.iSpin2D[iX][iY] * iNeighbors

             # Metropolis acceptance rule
        if iE_Delta <= 0:
            # If energy decreases (or unchanged), always accept flip
            self.iE_Total += iE_Delta
            self.iSpin2D[iX][iY] = -self.iSpin2D[iX][iY]
        else:
            # Otherwise accept with probability exp(-ΔE / T)
            dProb = rd.random()

            # Since ΔE can only be 4J or 8J here,
            # index = ΔE / (4J) gives 1 or 2
            dIndex = int(iE_Delta / (4 * self.J))

            if dProb < self.dProb_Target[dIndex]:
                self.iE_Total += iE_Delta
                self.iSpin2D[iX][iY] = -self.iSpin2D[iX][iY]
    def Sweep(self):
        for _ in range(self.SIZE * self.SIZE):
            self.MetropolisStep()

#            for iStep in range(1, self.iStep_Max):
#                # Record current energy before attempting flip
#                self.iEnergy_list.append(self.iE_Total)
                # Record current magnetization (normalized)
#                self.dMagnetization_list.append(self.CalculateMagnetization())


#    def Run(self):

        #print("Initial energy = ", self.iE_Total)
        #print("Initial magnetization = ", self.CalculateMagnetization())
#        self.PrintSpin()

#        for iStep in range(1, self.iStep_Max):
    
#            self.iEnergy_list.append(self.iE_Total)
#            self.dMagnetization_list.append(self.CalculateMagnetization())
    
#            self.Sweep()
        # Append final energy and magnetization
#        self.iEnergy_list.append(self.iE_Total)
#        self.dMagnetization_list.append(self.CalculateMagnetization())
    
       # print("Final energy = ", self.iE_Total)
        #print("Final magnetization = ", self.CalculateMagnetization())
#        self.PrintSpin()

    def Run(self, record=True):
        for iStep in range(1, self.iStep_Max):
            self.Sweep()
    
            if record:
                self.iEnergy_list.append(self.iE_Total)
                self.dMagnetization_list.append(self.CalculateMagnetization())
    
        if record:
            self.iEnergy_list.append(self.iE_Total)
            self.dMagnetization_list.append(self.CalculateMagnetization())
# ------------------------------------------------------------
# Energy calculation (Hamiltonian)
# ------------------------------------------------------------
    def CalculateEnergy(self):
        """
        Compute total energy of current lattice configuration.
        Hamiltonian:
            E = -J Σ s_i s_j
        Each bond is counted twice in naive summation,
        so we divide final result by 2.
        """
        iE_Total = 0

        for i in range(self.SIZE):
            for j in range(self.SIZE):
                # Sum of 4 nearest neighbors with periodic boundary
                iNeighbors = (
                    self.iSpin2D[(i - 1 + self.SIZE) % self.SIZE][j]
                    + self.iSpin2D[(i + 1) % self.SIZE][j]
                    + self.iSpin2D[i][(j - 1 + self.SIZE) % self.SIZE]
                    + self.iSpin2D[i][(j + 1) % self.SIZE]
                )

                iE_Total -= self.J * self.iSpin2D[i][j] * iNeighbors

        # Divide by 2 to correct double counting of bonds
        return iE_Total / 2
# ------------------------------------------------------------
# Magnetization calculation
# ------------------------------------------------------------
    def CalculateMagnetization(self):
        """
        Compute normalized magnetization (average spin per site).
        Returns a float in the range [-1, +1].
        """
        iTotal = 0
        for i in range(self.SIZE):
            for j in range(self.SIZE):
                iTotal += self.iSpin2D[i][j]

        return iTotal / (self.SIZE * self.SIZE)
# ------------------------------------------------------------
# Utility: Print spin lattice
# ------------------------------------------------------------
    def PrintSpin(self):
        """
        Print the current spin configuration in matrix form.
        """
        for i in range(self.SIZE):
            for j in range(self.SIZE):
                print(f"{self.iSpin2D[i][j]:2d}", end=" ")
            print()
# ------------------------------------------------------------
# Execution entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    # Create model instance (simulation runs automatically)
    model = Ising_Model()
    model.Run()
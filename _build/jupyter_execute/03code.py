#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')


# # The code
# 
# ## Introduction and basic structure
# 
# The main objective of this code is to implement the strategy presented in the previous sections to generate a library for SINDy. In order to do so the final output needs to be a list of callable functions.
# 
# In order to keep the code as tidy as possible the computation of the coefficients and that of the alpha matrix have been separated in two different modules. 
# 
# The main functionality of the codes are nothing more than a translation from matlab to python of the arbitrary Polynomial Chaos code (aPC) created by Oladyshkin {cite}`aPCMatlab`. The part that generates the list of functions has been, instead, written from scratch.

# ## The alpha Matrix
# 
# This first part of the code outputs the alpha matrix given the maximum degree of the expansion and the number of variables.

# In[2]:


import numpy as np
from math import factorial

def MultivariatePolynomialIndex(
        N, # number of variables
        d  # degree of the polynomial expansion 
        ):
    M = int(factorial(N+d) / (factorial(N)*factorial(d))) # total number of inputs
    
    PossibleDegree = [list(range(d+1)) for i in range(N)]
    UniqueDegreeCombinations = np.array(np.meshgrid(*PossibleDegree))
    UniqueDegreeCombinations = UniqueDegreeCombinations.reshape((N,-1))
    UniqueDegreeCombinations = UniqueDegreeCombinations.T
    
    DegreeWeight = np.zeros(len(UniqueDegreeCombinations))
    for i in range(len(UniqueDegreeCombinations)):
        DegreeWeight[i] = sum(UniqueDegreeCombinations[i,:])
            
    ID = np.argsort(DegreeWeight)
    SortDegreeCombination = UniqueDegreeCombinations[ID,:]
    
    return SortDegreeCombination[0:M,:]


# As an example the matrix for the 2-2 case, presented before, is calculated.

# In[3]:


d = 2 # maximum degree of the expansion
N = 2 # number of variables

alpha = MultivariatePolynomialIndex(N,d)

print(alpha)


# It has to be noted that a diffierence in the construction of the `UniqueDegreeCombinations` matrix has as a result a sligthly different ordering of the indices with respect to the original code. However, it doesn't effect the role of the matrix in the final expansion.

# ## The coefficients and the Library
# 
# The main part of the code computes the coefficients of the expansion and output the list of functions to be transformed in a SINDy library.
# 
# The main class of the code deals with the calculation of the coefficients and uses the previous code to compute the alpha matrix. This informations are enough to fully characterize the expansion. The last function can take this data and output the desired function list.

# In[4]:


import sympy as sym

class PolynomialChaos():
    '''
    Basic class for the Polynomial Chaos Expansion
    '''
    def __init__(self, 
                 distribution,
                 expansionDegree,
                 numberOfInputs):
        self.expansionDegree = expansionDegree
        self.distribution = distribution
        self.numberOfInputs = numberOfInputs
        
    def ComputeMoments(self, distribution_1D):
        '''
        Compute the statistical moments of a distribution (1D) up to the 
        2*expansionDegree (2*expansionDegree-1 would be sufficient, the last 
        could be useful for a furteher versions that implement normalization in 
        order to have orthonormal base)
        '''
        numberOfSamples = len(distribution_1D)
        DistributionMoments = np.array([np.sum(distribution_1D**i)/numberOfSamples for i in range((self.expansionDegree+1)*2)])
        return DistributionMoments
        
    def MomentMatrix(self, distribution_1D, polynomialDegree):
        '''
        Generate the moment matrix to compute the coefficients for a polynomial 
        of degree polynomialDegree, as explained in the reference paper [1]
        '''
        d = polynomialDegree + 1
        Hankel = np.zeros((d,d)) # moments matrix initialization
        moments = self.ComputeMoments(distribution_1D)
        for i in range(polynomialDegree+1):
            for j in range(polynomialDegree+1):
                if i < polynomialDegree:
                    Hankel[i,j] = moments[i+j]
                else:
                    Hankel[i,-1] = 1
        return Hankel
    
    def aPC_OneDimensional(self, distribution_1D):
        '''
        Computes and returns the coefficient matrix for a 1D distribution from the 0-degree
        polynomial up to the one of degree expansionDegree.
        '''
        d = self.expansionDegree + 1
        coefficients = np.zeros((d,d))
        for i in range(d):
            H = self.MomentMatrix(distribution_1D,i)
            v = np.zeros(i+1)
            v[-1] = 1
            coefficients[0:i+1,i] = np.linalg.solve(H,v)
        # coefficients = np.reshape(coefficients,(d,d,1))
        return coefficients

    def ComputeCoefficients(self):
        '''
        Computes the coefficient for the PC expansion (in general multidimensional).
        Makes use of the MultivariatePolynomialsIndex function to generate the 
        Alpha matrix useful for the construction of the base.
        The coefficient tensor and the Alpha matrix are enough to fully characterize
        the PC expansion
        
        The coefficient tensor has three dimensions, as:
            - the first represents the order of the sigle term (from 0 to expansionDegree)
            - the second the total degree of the polynomial (from 0 to expansionDegree)
            - the third the variable (from 1 to numberOfInputs)
        '''
        d = self.expansionDegree + 1
        if self.numberOfInputs == 1:
            self.coefficients = np.reshape(self.aPC_OneDimensional(self.distribution), (d,d,1))
            self.AlphaMatrix = np.array([range(d)]).T
        else:
            self.coefficients = np.zeros((d,d, numberOfInputs))
            for i in range(numberOfInputs):
                self.coefficients[:,:,i] = self.aPC_OneDimensional(self.distribution[:,i])
            self.AlphaMatrix = MultivariatePolynomialIndex(numberOfInputs, d-1)

            
def GenerateLibraryList(
        expansionDegree,
        coefficients,
        AlphaMatrix
        ):
    '''
    Given the Alpha matrix and the coefficient tensor conputes a list of functions
    ready to be transformed into a SINDy library.
    '''
    M , numberOfInputs = AlphaMatrix.shape # M = total number of terms in the expansion
    x = []
    for i in range(numberOfInputs): x.append(sym.symbols(f'x{i}')) # list of symbolic variables
    LibraryList = []
    for i in range(M): # order
        index = AlphaMatrix[i,:]
        MultivariatePolynomial = 1
        for j in range(numberOfInputs): # variable
            coeff = coefficients[:, index[j], j] 
            coeff = np.flip(coeff) # The MultivariatePolynomials function gives the coefficients from 0 to max_deg, while Poly starts from max_deg and goes to 0
            Polynomial1D = sym.Poly(coeff, x[j])
            MultivariatePolynomial = MultivariatePolynomial * Polynomial1D # multivaried polynomial object
            MultivariatePolynomial = MultivariatePolynomial.as_expr() 
            
        LibraryList.append(sym.lambdify(x, MultivariatePolynomial, 'numpy'))
        
    return LibraryList
        


# As an example, the coefficients for the Legendre polynomials are calculated from an evely distributed $x$ over $\mathcal X \in [-1,1]$.

# In[5]:


import matplotlib.pyplot as plt
from matplotlib import cm

n = 1000 # number of samples
data_uniform = np.linspace(-1,1,n)
data_uniform = np.array([data_uniform]) # to have the correct shape

expansionDegree = 5
numberOfInputs = 1

aPC = PolynomialChaos(data_uniform, expansionDegree, numberOfInputs)
aPC.ComputeCoefficients()
coefficients = aPC.coefficients
A = aPC.AlphaMatrix

LibraryList = GenerateLibraryList(5, coefficients, A)

print(coefficients[:,:,0].T)
fig, ax = plt.subplots()
for i in range(1,6):
    ax.plot(data_uniform[0], LibraryList[i](data_uniform).T, '-')
ax.legend(['deg = 1','deg = 2','deg = 3','deg = 4','deg = 5'])
ax.set_xlabel('x')
ax.set_title('Legendre Polynomials')
plt.show()


# Note that the coefficient are the same as the classical ones up to a multiplicative constant, since the original Legendre polynomials are meant to be equal to $\pm1$ for $x=\pm1$
# 
# As a proof of concept a second variable following a gaussian distribution will be added. This one should generate Hermite Polynomials. The overall expansion will have 2 variables polynomials.

# In[6]:


np.random.seed(43)
data_gaussian = np.random.randn(n)
x_min = np.min(data_gaussian)
data = np.zeros((n,2))
data[:,0] = data_uniform
data[:,1] = data_gaussian

numberOfInputs = 2

aPC = PolynomialChaos(data, expansionDegree, numberOfInputs)
aPC.ComputeCoefficients()
coefficients = aPC.coefficients
A = aPC.AlphaMatrix

LibraryList = GenerateLibraryList(5, coefficients, A)

# FIGURES

x_min = np.min(data_gaussian)
x_max = np.max(data_gaussian)
X, Y = np.meshgrid(data_uniform, np.linspace(x_min,x_max,n)) 
Z = LibraryList[-1](X,Y)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X,Y,Z, cmap=cm.coolwarm)
ax.set_xlabel('x uniform')
ax.set_ylabel('x gaussian')
ax.set_title(f'Polynomial {A[-1]}')

plt.show()
print('All terms in the expansion:')
print(A)


# In[ ]:





# Polynomial Chaos Expansion: better and faster

## What changes in this version?

In the previous jupyter notebook the necessity of using an orthogonal set of basis functions has been outlined and a first code for computing them was  added in. The code worked for different propability distributions and was able to retrieve the known Legendre polynomials where provided the right data (a uniformly distributed vector in [-1,1]). However in was very inefficient computationally and, more importantly, couldn't handle multivariate polynomials, which will be necessary when dealing with the Theodorsen funtion identification (since there will be at least two inputs, namely $\ddot{\alpha}$ and $\ddot{h}$). 

The paper by Oladyshkin {cite}`oladyshkin2012data` provides a much better framework for calculating the expansion coefficients and hints to a way of keeping track of the indices for multivariate polynomials. Thus the necessity to rewrite the previous code. 

After a brief introduction to the new way of computing the coefficients and of the index system the code will be displayed in the last section and some example are going to be made to illustrate how the code works and the precision in identifying known polynomials.

```{tableofcontents}
```

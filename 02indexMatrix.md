# The alpha matrix: dealing with multivariate polynomials

The previous section only deals with single valued random variables. But what happens when the input is multidimensional?

Let the input be $X = [x_0, x_1, \dots, x_N]\in\mathcal X ^N$. Under the assumption that all $x_i$ are independent, the resulting multivariate polynomial expansion of maximum degree $d$ is just the multiplication of the single variable ones:

$$
\begin{gather}
	Y(x) = \sum_i \beta_i\Phi_i(X)\\
	\Phi_i(X) = \prod_{j= 0}^N \pi_{\alpha^i_j}^j(X)\\
	\mbox{with: }\sum_{j=0}^N \alpha^i_j \leq d
\end{gather}
$$

where the index matrix $\alpha^i_j$ contains all the combinatorial information needed to accurately multiply the right coefficients. It must be noted that now the number of terms in the expansion is increased to $M = (N+d)!/(N!d!)$ to accomodate for all the possible combinations of order less than $d$.

The matrix can be constructed in the following way. We let every column represent a variable and every row a term in the expansion. In this way the entire matrix will be of dimension $d\times N$. Every entry $\alpha^i_j$ will represent the degree of the polynomial of the variable $j$ in the term $i$. An example will hopefully clarify this idea.

Let $X = [x_0, x_1]$ and the maximum degree of the expansion be $M = 2$. Let the polynomial of degree $i$ of the variable $x_j$ be called $\pi_i^j$. This will generate an expansion with $M = 6$ terms, in the form:

$$
	Y(X) = \beta_0\pi_0^0\pi_0^1 + \beta_1\pi_1^0\pi_0^1 + \beta_2\pi_0^0\pi_1^1 + \beta_3\pi_2^0\pi_0^1 + \beta_4\pi_1^0\pi_1^1 + \beta_5\pi_0^0\pi_2^1
$$

The alpha matrix that indicates such a combination would thus be:

$$
\alpha = \begin{bmatrix}
	0&0\\1&0\\0&1\\2&0\\1&1\\0&2
\end{bmatrix}
$$
/* MIPoD 1.0 basic matrix operations header file */
/* copyright 2007 Paul Hohenlohe */
/* permission is granted to use and modify this software _provided that_ any
publication resulting from use of this software or any part or modification thereof
contains at least one of the following citations:

Hohenlohe, P.A.  2007.  MIPoD: Microevolutionary Inference from Patterns of Divergence.
Software available at http://oregonstate.edu/~hohenlop/software

Hohenlohe, P.A. & S.J. Arnold.  2008.  MIPoD: A Hypothesis-Testing Framework for Microevolutionary
Inference from Patterns of Divergence.  American Naturalist 171(3): 366-385.

questions/comments: hohenlop@science.oregonstate.edu */

/* Functions defined herein:
minimum(x,y)
maximum(x,y)
kroneckerdelta(x,y)
bn logdeterminant(double **a, int n)
int SquareKroneckerProduct(double **a, double **b, double **c, int n, int m)
int LUinverse(double **a, double **b, int n)
int linearsolve(double **mat, double *x, double *y, int n)
int rowechelon(double **mat, int r, int c)
double chisquare(double x, int d)
int eigen(double **mat, double **out, int n)
int inveigen(double **evecs, double **out, int n, int y) 
int ordereigen(double **evecs, int n) */

#define minimum(x,y) ( ((x) < (y)) ? (x) : (y) )
#define maximum(x,y) ( ((x) > (y)) ? (x) : (y) )
#define kroneckerdelta(x,y) ( ((x) == (y)) ? (1) : (0) )

typedef struct {
	int i;
	double d;
} bn;  /* store big numbers as i*e^d, so d = ln(|bn|) and i is the sign */

bn logdeterminant(double **a, int n) {
	/* uses PLU decomposition algorithm from Demmel 1997; returns determinant of a as a bn */
	int i, j, k, sign, zero, row, col;
	double m, tmp, **b;
	bn det;
	if ((b = malloc((n+1)*sizeof(double *))) == NULL) {
		printf ("memory error"); exit(0); }
	for (i=0; i<(n+1); i++)
		if ((b[i] = malloc((n+1)*sizeof(double))) == NULL) {
			printf ("memory error"); exit(0); }
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			b[i][j] = a[i][j];
	zero = 0; sign = 1;	
	for(i=0; i<n; i++) {
		m=0; row=i; col=i;
		for(k=i; k<n; k++) {
			for(j=i; j<n; j++) {
				tmp = fabs(b[k][j]);
				if(tmp > m) {
					m = tmp;
					row = k; col = j;
				}
			}
		}
		if(m == 0) {
			zero = 1; continue; }
		if(row != i) {
			for(j=0; j<n; j++) {
				m = b[row][j]; b[row][j] = b[i][j]; b[i][j] = m; }
			sign *= -1;
		}
		if(col != i) {
			for(j=0; j<n; j++) {
				m = b[j][col]; b[j][col] = b[j][i]; b[j][i] = m; }
			sign *= -1;
		}
		for(j=i+1; j<n; j++)
			b[j][i] /= b[i][i];
		for(j=i+1; j<n; j++) {
			for(k=i+1; k<n; k++)
				b[j][k] -= b[j][i]*b[i][k];
		}
	}
	if(zero == 1) {
		det.i = 0; det.d = 0;
	}
	else {
		det.d = 0;
		for (i=0; i<n; i++) {
			if (b[i][i] == 0) {
				det.i = 0; det.d = 0; continue; }
			else {
				det.d += log(fabs(b[i][i]));
				if (b[i][i] < 0) sign *= -1;
			}
		}
		if(sign > 0) det.i = 1;
		else det.i = -1;
	}
	for (i=0; i<(n+1); i++) free(b[i]);
	free (b);
	return (det);
}

int SquareKroneckerProduct(double **a, double **b, double **c, int n, int m) {
	/* Populates square matrix c (size (n*m)x(n*m)) with the Kronecker product 
	of square matrices a (size nxn) and b (size mxm) */
	int i,j,k,l;
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			for(k=0; k<m; k++)
				for(l=0; l<m; l++)
					c[i*m + k][j*m + l] = a[i][j]*b[k][l];
	return 0;				
}

int LUinverse(double **a, double **b, int n) {  
	/* Inverts the matrix a into the matrix b */
	int i, j, k, q, row, col, *rows, *cols, zero;
	double **UL, *temprow, m, tmp;
	if ((rows = malloc(n*sizeof(int))) == NULL) {
		printf("memory error"); exit(0); }
	if ((cols = malloc(n*sizeof(int))) == NULL) {
		printf("memory error"); exit(0); }
	if ((temprow = malloc(n*sizeof(double))) == NULL) {
		printf("memory error"); exit(0); }
	if ((UL = malloc(n*sizeof(double *))) == NULL) {
		printf("memory error"); exit(0); }
	for (i=0; i<n; i++) if ((UL[i] = malloc(n*sizeof(double))) == NULL) {
		printf("memory error"); exit(0); }
	for(i=0; i<n; i++) {
		rows[i] = i; cols[i] = i;
		for(j=0; j<n; j++) b[i][j] = a[i][j];
	}
	zero = 0;
	for(i=0; i<n; i++) {
		m=0; row=i; col=i;
		for(k=i; k<n; k++) {  /* Find largest element */
			for(j=i; j<n; j++) {
				tmp = fabs(b[k][j]);
				if(tmp > m) {m = tmp; row = k; col = j; }
			}
		}
		if(m == 0) {zero = 1; continue; }
		if(row != i) {  /*switch rows */
			for(k=0; k<n; k++) {temprow[k] = b[row][k]; b[row][k] = b[i][k]; b[i][k] = temprow[k]; }
			q = rows[i]; rows[i] = rows[row]; rows[row] = q;
		}
		if(col != i) {  /*switch columns */
			for(j=0; j<n; j++) {
				m = b[j][col]; b[j][col] = b[j][i]; b[j][i] = m;
			}
			q = cols[i]; cols[i] = cols[col]; cols[col] = q;
		}
		for(j=i+1; j<n; j++)  /* LU decomposition */
			b[j][i] = b[j][i]/b[i][i];
		for(j=i+1; j<n; j++) {
			for(k=i+1; k<n; k++)
				b[j][k] = b[j][k] - b[j][i]*b[i][k];
		}
	}
	if(zero == 1) return zero;
	for(i=1; i<n; i++) { /* invert lower factor */
		for(j=0; j<i; j++) {
			tmp = 0; for(k=(j+1); k<i; k++) tmp += b[i][k]*UL[k][j];
			UL[i][j] = - tmp - b[i][j];
		}
	}
	for(i=(n-1); i>=0; i--) { /* invert upper factor */
		UL[i][i] = 1/b[i][i];
		for(j=(i+1); j<n; j++) {
			tmp = 0; for(k=(i+1); k<(j+1); k++) tmp += b[i][k]*UL[k][j];
			UL[i][j] = -1 * UL[i][i] * tmp;
		}
	}
	for(i=0; i<n; i++) { /*multiply inverses of U and L */
		for(j=0; j<n; j++) {
			tmp = 0; for(k=maximum((j+1),i); k<n; k++) tmp += UL[i][k]*UL[k][j];
			if(i<=j) b[i][j] = UL[i][j] + tmp;
			else b[i][j] = tmp;
		}
	}
	for(i=0; i<n; i++) for(j=0; j<n; j++) UL[cols[i]][j] = b[i][j];
	for(j=0; j<n; j++) for(i=0; i<n; i++) b[i][rows[j]] = UL[i][j];
	for (i=0; i<n; i++) free(UL[i]);
	free(rows); free (cols); free(temprow); free (UL);
	return zero;
}

int rowechelon(double **mat, int r, int c) {
	/* returns reduced row echelon form of matrix by Gaussian elimination */
	int i, j, k, l;  
	double tmp, *temprow;
	if ((temprow = malloc(c*sizeof(double))) == NULL) {
		printf("memory error"); exit(0); }
	i=0; j=0;
	while (i < minimum(r,(c-1))) {
		while (mat[i][j] == 0 && j < c) {
			k = i+1;
			while (mat[k][i] == 0 && k < r)
				k++;
			if (k < r) {
				temprow = mat[i];
				mat[i] = mat[k];
				mat[k] = temprow;
			}
			else j++;
		}
		tmp = mat[i][j];
		for(l=j; l<c; l++)
			mat[i][l] /= tmp;
		for(k=0; k<r; k++)
			if(k != i && mat[k][j] != 0) {
				tmp = mat[k][j];
				for(l=j; l<c; l++)
					mat[k][l] -= tmp*mat[i][l];
			}
		i++; if(j < i) j=i;
	}
	free(temprow);
	return 0;
}

int linearsolve(double **mat, double *x, double *y, int n) {
	/* Solves the linear equations mat*x=y for x */
	int i, j, zero;   
	double **tempmat;
	if ((tempmat = malloc(n*sizeof(double *))) == NULL) {
		printf("memory error"); exit(0); }
	for (i=0; i<n; i++)
		if ((tempmat[i] = malloc((1+n)*sizeof(double))) == NULL) {
			printf ("memory error"); exit(0); }
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++)
			tempmat[i][j] = mat[i][j];
		tempmat[i][n] = y[i];
	}
	zero = rowechelon(tempmat, n, (n+1));
	for(i=0; i<n; i++) x[i] = tempmat[i][n];
	for (i=0; i<n; i++) free(tempmat[i]); free (tempmat);
	return zero;
}

double chisquare(double x, int d) { 
	/* Returns p-value for chi-square=x and df=d */
	if (x <= 0 || d < 1) return 1.0;
	int i, j;
	double a, b, q, sum, prod, out;
	a = (double) d/2; b = x/2;
	sum = 0;
	if ((int) a == a) { /* If df is even, p can be calculated exactly: */
		for(i=0; i<a; i++) {
			prod = 1;
			for(j=1; j<(i+1); j++)
				prod *= j;
			sum += pow(b,i)/prod;
		}
		out = sum/exp(b);
	}
	else {  /* If df is odd, use converging approximation: */
		if(x < 38) {  /* small x -- the cut-off at 38 is arbitrary but works well */
			for(i=0; i<100; i++) {
				prod = 1;
				for(j=1; j<(i+1); j++)
					prod *= b/j;
				sum += pow(-1,i)*prod/(2*i + 1);
			}
			sum *= 2*sqrt(b);
			if (d == 1) out = 1 - sum/sqrt(M_PI);
		}
		else {  /* large x */
			for(i=0; i<100; i++) {
				prod = 1;
				for(j=1; j<(i+1); j++)
					prod *= (2*j - 1)/(2*x);
				sum += pow(-1,i)*prod;
			}
			out = (1 + sum)/(exp(b)*sqrt(M_PI*b));
			sum = sqrt(M_PI)*(1 - out);
		}
		if (d > 1) { /* If df is odd and >1, sum forms the kernel of larger exact calc. */
			for(q=0.5; q<a; q++)
				sum = q*sum - pow(b,q)/exp(b);
			prod = sqrt(M_PI);
			for(q=0.5; q<a; q++)
				prod *= q;
			out = 1 - sum/prod;
		}
	}
	return out;
}

int eigen(double **mat, double **out, int n) {
	/* returns eigenvectors and eigenvalues of mat as columns of out -- symmetric only!!! */
	int i, j, k, l;   
	double c, s, tmp, tmp1, rot, **tmpmat;
	if ((tmpmat = malloc((n+1)*sizeof(double *))) == NULL) {
		printf("memory error"); exit(0); }
	for (i=0; i<(n+1); i++)
		if ((tmpmat[i] = malloc((n+1)*sizeof(double))) == NULL) {
			printf ("memory error"); exit(0); }
	for(i=0; i<n; i++)
		for(j=(i+1); j<n; j++)
			if(mat[i][j] != mat[j][i]) return 1;  /* return 1 if not symmetric */
	if (n == 1) {out[0][0] = 1; out[0][1] = mat[0][0]; }
	else if (n == 2) {
		tmp = sqrt(4*mat[0][1]*mat[0][1] + (mat[0][0] - mat[1][1])*(mat[0][0] - mat[1][1]));
		out[2][0] = 0.5*(mat[0][0] + mat[1][1] + tmp);
		out[2][1] = 0.5*(mat[0][0] + mat[1][1] - tmp);
		for(i=0; i<2; i++) {
			if(mat[0][1] != 0) {
				tmp = (out[2][i] - mat[0][0])/mat[0][1];
				tmp1 = sqrt(1 + tmp*tmp);
				out[0][i] = 1/tmp1; out[1][i] = tmp/tmp1;
			}
			else { out[i][i] = 1; out[1-i][1-i] = 0; }
		}
	}
	else {  /* Jacobi algorithm for finding eigenvalues */
		for(i=0; i<n; i++)
			for(j=0; j<n; j++)
				tmpmat[i][j] = mat[i][j];
		k = 0; l = 1; tmp = tmpmat[k][l];
		for (i=0; i<n; i++) {
			for (j=(i+1); j<n; j++) {
				if (fabs(tmpmat[i][j]) > fabs(tmp)) {
					k = i; l = j; tmp = tmpmat[i][j]; }
			}
		}
		while (fabs(tmp) > 0.00001*sqrt(fabs(tmpmat[k][k]*tmpmat[l][l]))) {
			if (tmpmat[k][k] == tmpmat[l][l]) rot = 0.25*M_PI;
			else rot = 0.5*atan(2*tmpmat[k][l]/(tmpmat[k][k]-tmpmat[l][l]));
			c = cos(rot); s = sin(rot);  
			/* lower triangle is rotated, then transfered to upper; k<l */
			for (i=0; i<k; i++) {
				tmpmat[k][i] = c*tmpmat[i][k] + s*tmpmat[i][l];
				tmpmat[l][i] = -1*s*tmpmat[i][k] + c*tmpmat[i][l];
			}
			for (i=(k+1); i<l; i++) {
				tmpmat[i][k] = c*tmpmat[k][i] + s*tmpmat[i][l];
				tmpmat[l][i] = -1*s*tmpmat[k][i] + c*tmpmat[i][l];
			}
			for (i=(l+1); i<n; i++) {
				tmpmat[i][k] = c*tmpmat[k][i] + s*tmpmat[l][i];
				tmpmat[i][l] = -1*s*tmpmat[k][i] + c*tmpmat[l][i];
			}
			tmp = tmpmat[k][k]; tmp1 = tmpmat[l][l];
			tmpmat[l][k] = (c*c - s*s)*tmpmat[k][l] + c*s*(tmp1 - tmp);
			tmpmat[k][k] = c*c*tmp + 2*c*s*tmpmat[k][l] + s*s*tmp1;
			tmpmat[l][l] = s*s*tmp - 2*s*c*tmpmat[k][l] + c*c*tmp1;
			for (i=0; i<(n-1); i++)
				for (j=(i+1); j<n; j++)
					tmpmat[i][j] = tmpmat[j][i];
			k = 0; l = 1; tmp = tmpmat[k][l];
			for (i=0; i<n; i++) {
				for (j=(i+1); j<n; j++) {
					if (fabs(tmpmat[i][j]) > fabs(tmp)) {
						k = i; l = j; tmp = tmpmat[i][j]; }
				}
			}
		}
		for(i=0; i<n; i++) {  /* order eigenvalues by absolute value */
			k=0; tmp = fabs(tmpmat[k][k]);
			for(j=1; j<n; j++) {
				if(fabs(tmpmat[j][j]) > tmp) {
					tmp = fabs(tmpmat[j][j]); k = j;
				}
			}	
			out[n][i] = tmpmat[k][k];
			tmpmat[k][k] = 0;
		}
		for(i=0; i<n; i++)  {  /* Calculate eigenvectors */
			for(j=0; j<n; j++) {
				for(k=0; k<n; k++)
					tmpmat[j+1][k] = mat[j][k];
				tmpmat[j+1][j] -= out[n][i];
				tmpmat[0][j] = 1;
				tmpmat[j+1][n] = 0;
			}
			tmpmat[0][n] = 1;
			j = rowechelon(tmpmat, (n+1), (n+1));
			tmp = 0;
			for(j=0; j<n; j++) tmp += tmpmat[j][n]*tmpmat[j][n];
			for(j=0; j<n; j++)
				out[j][i] = tmpmat[j][n]/sqrt(tmp);
		}
	}
	for(i=0; i<(n+1); i++) free(tmpmat[i]); free(tmpmat);
	return 0;
}

int inveigen(double **evecs, double **out, int n, int y) {
	/* Takes eigenvectors and eigenvalues from eigen and produces the matrix out */
	/* If y=0, evecs is an orthonormal basis, i.e. the matrix is real symmetric, and calculation is quicker */
	int i, j, k; 
	double **inv, **tmp, sum;
	if ((inv = malloc(n*sizeof(double *))) == NULL) {
		printf("memory error"); exit(0); }
	for (i=0; i<n; i++)
		if ((inv[i] = malloc(n*sizeof(double))) == NULL) {
			printf ("memory error"); exit(0); }
	if ((tmp = malloc(n*sizeof(double *))) == NULL) {
		printf("memory error"); exit(0); }
	for (i=0; i<n; i++)
		if ((tmp[i] = malloc(n*sizeof(double))) == NULL) {
			printf ("memory error"); exit(0); }
	for(i=0; i<n; i++) 
		for(j=0; j<n; j++)
			tmp[i][j] = evecs[i][j]*evecs[n][j];
	if (y != 0) {
		if(LUinverse(evecs, inv, n) == 1) return 1; }
	else {
		for(i=0; i<n; i++)
			for(j=0; j<n; j++)
				inv[i][j] = evecs[j][i];
	}
	for(i=0; i<n; i++)
		for(j=0; j<n; j++) {
			sum = 0;
			for(k=0; k<n; k++)
				sum += tmp[i][k]*inv[k][j];
			out[i][j] = sum;
		}
	for(i=0; i<n; i++) {free(inv[i]); free(tmp[i]); }
	free(inv); free(tmp);
	return 0;
}

int ordereigen(double **evecs, int n) {
	/* Order eigenvectors in evecs from largest to smallest abs value of e-value
	evecs stored as e-vectors in columns, corresponding e-value in evecs[i][n] */
	int i,j,k,l;
	double tmp, **tmpevecs;
	if ((tmpevecs = malloc((n+1)*sizeof(double *))) == NULL) {
		printf("memory error"); exit(0); }
	for (i=0; i<(n+1); i++)
		if ((tmpevecs[i] = malloc((n+1)*sizeof(double))) == NULL) {
			printf ("memory error"); exit(0); }
	for(i=0; i<(n+1); i++)
		for(j=0; j<n; j++)
			tmpevecs[i][j] = evecs[i][j];
	l = n; /* l is number of non-zero evals */
	for(i=0; i<n; i++) { /* scan first for evals == 0 */
		if(tmpevecs[n][i] == 0) {
			for(j=0; j<(n+1); j++)
				evecs[j][l-1] = tmpevecs[j][i];
			l--;
		}
	}
	for(i=0; i<l; i++) {
		tmp = 0; k = 0;
		for(j=0; j<n; j++) {
			if(fabs(tmpevecs[n][j]) > tmp) {
				tmp = fabs(tmpevecs[n][j]);
				k = j;
			}
		}
		for(j=0; j<(n+1); j++)
			evecs[j][i] = tmpevecs[j][k];
		tmpevecs[n][k] = 0;
	}
	for(i=0; i<(n+1); i++) free(tmpevecs[i]);
	free(tmpevecs);
	return 0;
}


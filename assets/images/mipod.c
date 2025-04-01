/* MIPoD 1.0 Neutral Module */
/* copyright 2008 Paul Hohenlohe */
/* permission is granted to use and modify this software _provided that_ any
* publication resulting from use of this software or any part or modification thereof
* contains at least one of the following citations:
*
* Hohenlohe, P.A.  2007.  MIPoD: Microevolutionary Inference from Patterns of Divergence.
* Software available at http://oregonstate.edu/~hohenlop/software
*
* Hohenlohe, P.A. & S.J. Arnold.  2008.  MIPoD: A Hypothesis-Testing Framework for Microevolutionary
* Inference from Patterns of Divergence.  American Naturalist 171(3):366-385.
*
* questions/comments: hohenlop@science.oregonstate.edu */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MatrixFunx.h"

/* structures */
typedef struct {
	double branch;
	int parent;
} node;
typedef struct {
	char name[10];
	double *trait;
	double branch;
	int parent;
} taxon;

/*functions: */
void getdata();
void output();
double lnlikelihood(double **Inv, double N, bn deter);
double maxNe(double **Inv);
double confN(double **Inv, double N, bn deter, int d);
double est_sigma(double **ellipse, double **evecs, double N, int y);
double est_epsilon(double **ellipse, double **evecs, double N, int which, int y);
double est_phi(double **ellipse, double **evecs, double N, int which);
int rotate(double **evecs, double phi, int k, int l, int n);
double vecmatvec(double **Inv);
int populate(double **G, double **Inv, bn *deter);
int mattopar(double **mat, double **ellipse, double **evecs, int n);
int partomat(double **mat, double **ellipse, double **evecs, int n);
void memerror();

/*global variables -- these are the DATA from the user and don't change: */
int spp, traits, sizeA;
double gens, Ne, acc, *testphi, *mean;
taxon *taxa;
node *nodes;
double **Tree, **invTree;
bn Treedeter;
double **Gmat;

int main() {
	acc = 3;
	printf("\n\n** MIPoD 1.0 Neutral Module **\n\n");
	getdata();		/* Read files */
	output();		/* Analyze and print to Outfile.txt */	
	return 0;
}

void getdata() {  /*This function reads in the data and populates the global variables */
	int i, j, k, l, q;
	double d, e, f;
	i=j=l=0;
	char *datafilename;
	FILE *fd;
	datafilename = (char *)malloc(20*sizeof(char));
	printf ("\nWhat is the name of the data file?");
	scanf ("%s", datafilename);
	fd = fopen (datafilename,"r");
	if (!fd) {
		printf ("Can't find the file\n");
		exit (0);
	}
	printf ("Reading data from the file %s ...\n", datafilename);
	fscanf (fd, "%d", &spp);
	fscanf (fd, "%d", &traits);
	fscanf (fd, "%lf", &Ne);
	fscanf (fd, "%lf", &gens);
	sizeA = spp*traits;	
		/* Assign memory for global variables: */
	node *tmpnode;
	if ((tmpnode = realloc(nodes, (spp - 1)*sizeof(node))) == NULL) memerror();
		nodes = tmpnode;
	taxon *tmptaxa;
	if ((tmptaxa = realloc(taxa, spp*sizeof(taxon))) == NULL) memerror();
	for(i=0; i<spp; i++) if ((tmptaxa[i].trait = malloc(traits*sizeof(double))) == NULL) memerror();
	taxa = tmptaxa;
	double **tmp;	
	if ((tmp = realloc(Tree, spp*sizeof(double *))) == NULL) memerror();
	for (i=0; i<spp; i++) if ((tmp[i] = malloc(spp*sizeof(double))) == NULL) memerror();
	Tree = tmp;
	if ((tmp = realloc(invTree, spp*sizeof(double *))) == NULL) memerror();
	for (i=0; i<spp; i++) if ((tmp[i] = malloc(spp*sizeof(double))) == NULL) memerror();
	invTree = tmp;
	if ((tmp = realloc(Gmat, traits*sizeof(double *))) == NULL) memerror();
	for (i=0; i<traits; i++) if ((tmp[i] = malloc(traits*sizeof(double))) == NULL) memerror();
	Gmat = tmp;
	if ((testphi = realloc(testphi, traits*sizeof(double))) == NULL) memerror();
	if ((mean = realloc(mean, traits*sizeof(double))) == NULL) memerror();
		/* Read testphi vector */
	for (i=0; i<traits; i++)	
		fscanf (fd, "%lf", &testphi[i]);
		/* Read the G-matrix */
	fscanf (fd, "%lf", &Gmat[0][0]);
	if(Gmat[0][0] != 0) {
		for (j=1; j<traits; j++)
			fscanf (fd, "%lf", &Gmat[0][j]);
		for (i=1; i<traits; i++) 
			for (j=0; j<traits; j++)
				fscanf (fd,"%lf",&Gmat[i][j]);
	}
		/* Read the tree */
	char inc; 
	int whichnode[spp - 1];
	for (i=0; i<(spp-1); i++)
		whichnode[i] = 0;
	inc = fgetc(fd);
	i = q = 0;
	k = -1;
	while (inc != ';') {
		if (inc == '(') {
			k++;
			whichnode[k] = q;
			l = 1; q++;  /*expect a species name next */
		}
		if (inc == ':') { /* info for internal node */
			fscanf (fd, "%lf", &nodes[whichnode[k+1]].branch);
			nodes[whichnode[k+1]].parent = whichnode[k];
		}	
		if (l == 1 && inc != '(') {
			j = 0;
			while (inc != ':') {  /*read the species name */
				taxa[i].name[j] = inc;
				j++;
				inc = fgetc (fd);
			}
			taxa[i].parent = whichnode[k]; /* info for taxon i */
			if (inc == ':')
				fscanf (fd, "%lf", &taxa[i].branch);
			l = 0; i++;
		}
		if (inc == ')') k--;
		if (inc == ',') l = 1;     /*expect a species name next */
		inc = fgetc(fd);
	}
		/* Build the tree matrix; q = number of nodes */
	int anc[spp][q]; 
	for (i=0; i<spp; i++) { 
		for (j=0; j<spp; j++) Tree[i][j] = 0;
		for (j=0; j<q; j++) anc[i][j] = 0;
	}
	for (i=0; i<spp; i++) Tree[i][i] += taxa[i].branch;
	for (i=0; i<spp; i++) {
		if ((j = taxa[i].parent) == 0) continue;
		while (j != 0) {
			anc[i][j] = 1;
			j = nodes[j].parent;
		}		
	}
	for (i=0; i<spp; i++) {
		for (k=0; k<q; k++) {
			if (anc[i][k] == 1) Tree[i][i] += nodes[k].branch;
			for (j=(i+1); j<spp; j++) {
				if (anc[i][k] == 1 && anc[j][k] == 1) {
					Tree[i][j] += nodes[k].branch;
					Tree[j][i] += nodes[k].branch;
				}	
			}
		}
	}
	for(i=0; i<spp; i++)
		for(j=0; j<spp; j++)
			Tree[i][j] *= gens;
		/* Read trait data */
	char checkname[10];
	for(i=0; i<spp; i++) {
		fscanf(fd, "%s", checkname);
		if (strcmp(checkname, taxa[i].name) != 0) {
			j=0;
			while (strcmp(checkname, taxa[j].name) != 0) {
				if(j > (spp-1)) {
					printf ("Can't find the taxon %s in the tree\n\n", checkname);
					exit(0);
				}
				j++;
			}
			for(k=0; k<traits; k++)
				fscanf(fd, "%lf", &taxa[j].trait[k]);
		}
		else for(k=0; k<traits; k++)
			fscanf(fd, "%lf", &taxa[i].trait[k]);
	}
	fclose (fd);
	i = LUinverse(Tree, invTree, spp);		
	if (i == 1) {printf ("Can't invert Tree matrix\n\n"); exit (0);}
	Treedeter = logdeterminant(Tree, spp);
		/* Calculate phylogenetically-corrected means */
	d = 0;
	for(i=0; i<spp; i++)
		for(j=0; j<spp; j++)
			d += invTree[i][j];
	for(i=0; i<traits; i++) {
		f = 0;
		for(j=0; j<spp; j++) {
			e = 0;
			for(k=0; k<spp; k++)
				e += invTree[j][k];
			f += e*taxa[j].trait[i];
		}
		mean[i] = f/d;
	}
	free(datafilename);
}

void output() {
	int i,j,k;
	FILE *f;
	double *like, popsize, tmpLR;
	double ***ell, ***evecs;
	/* vector with like[0] = user-supplied, like[1] = size, like[2]-like[traits] = shape, 
	* like[traits+1]-like[2*traits-1] = phi,  like[2*traits] = test phi */
	if ((like = malloc((2*traits+1)*sizeof(double))) == NULL) memerror();
	/* matrices with column 1 = (sigma, eps1, eps2,..), cols 0 and 2 = 95% conf limits
	*** NOTE: conf limits for epsilons await further statistical work ***
	*** done for sigma and for epsilon when traits == 2 ***
	* ell[0] = user-supplied, ell[1] = size, ell[2]-ell[traits] = shape, ell[traits+1]-ell[2*traits-1] = phi, ell[2*traits] = test phi */
	if ((ell = malloc((2*traits+1)*sizeof(double **))) == NULL) memerror();
	for (i=0; i<(2*traits+1); i++) { if ((ell[i] = malloc(traits*sizeof(double *))) == NULL) memerror();
		for (j=0; j<traits; j++) if ((ell[i][j] = malloc(3*sizeof(double))) == NULL) memerror(); }
	/* eigenvector/eigenvalue matrix evecs[traits][traits+1][traits]:
	* evecs[0] = user-input G, evecs[1] through evecs[traits-1] = ML estimates of phi, evecs[traits] = testphi */
	if ((evecs = malloc((traits+1)*sizeof(double **))) == NULL) memerror();
	for(i=0; i<(traits+1); i++) { if ((evecs[i] = malloc((traits+1)*sizeof(double *))) == NULL) memerror();
		for (j=0; j<(traits+1); j++) if ((evecs[i][j] = malloc(traits*sizeof(double))) == NULL) memerror(); }

	if (Ne == 0) popsize = 1;
	else popsize = Ne;
	if (Gmat[0][0] == 0) {
		for(i=0; i<traits; i++)
			for(j=0; j<traits; j++)
				evecs[traits - 1][i][j] = kroneckerdelta(i,j);
		printf ("Estimating G-matrix parameters...\n");
		like[2*traits - 1] = est_phi(ell[2*traits - 1], evecs[traits-1], popsize, 0);
	}
	else if (Gmat[0][0] > 0) {
		if (mattopar(Gmat, ell[0], evecs[0], traits) == 1) {
			printf ("\nUnable to calculate eigenvectors of G\n\n\n"); exit(0); }
		if (Ne > 0) like[0] = est_sigma(ell[0], evecs[0], popsize, 0);
		printf ("Estimating size...\n");
		for (i=0; i<traits; i++) ell[1][i][1] = ell[0][i][1];
		like[1] = est_sigma(ell[1], evecs[0], popsize, 2);
		for(i=1; i<traits; i++) {
			printf ("Estimating shape(%d)...\n",i);
			for(j=0; j<traits; j++) ell[i+1][j][1] = ell[i][j][1];
			like[i+1] = est_epsilon(ell[i+1], evecs[0], popsize, (traits - i), 0);
		}
		for (i=1; i<traits; i++) {
			printf ("Estimating orientation(%d)...\n",i);		
			for(j=0; j<(traits+1); j++)
				for(k=0; k<traits; k++)
					evecs[i][j][k] = evecs[0][j][k];
			like[traits+i] = est_phi(ell[traits+i], evecs[i], popsize, (traits - 1 - i));		
		}
	}
	printf ("Comparing test vector...\n");
	for (i=0; i<(traits+1); i++)
		for (j=0; j<traits; j++)
			evecs[traits][i][j] = kroneckerdelta(i,j);
	tmpLR = 0;
	for (i=0; i<traits; i++) tmpLR += testphi[i]*testphi[i];
	tmpLR = sqrt(tmpLR);
	for (i=0; i<traits; i++) testphi[i] /= tmpLR;
	if (testphi[0] < 0) evecs[traits][0][0] = -1;
	for (i=1; i<traits; i++) {
		tmpLR = asin(testphi[i]);
		rotate(evecs[traits], tmpLR, 0, i, traits);
	}
	if (traits == 2) like[2*traits] = est_epsilon(ell[2*traits], evecs[traits], popsize, 1, 0);
	else like[2*traits] = est_phi(ell[2*traits], evecs[traits], popsize, 1);
		/* Display the results: */
	f = fopen ("Outfile.txt","w");
	fprintf (f,"%d taxa\n%d traits\n\n", spp, traits);
	if (Ne != 0)
		fprintf (f,"Ne = %4.2lf\n\n", Ne);
	if (Gmat[0][0] != 0) {
		fprintf (f,"Direct estimate of G: \n");
			for (i=0; i<traits; i++) {
				for (j=0; j<traits; j++)
					fprintf(f,"%4.2lf\t",Gmat[i][j]);
				fprintf(f,"\n");
			}
		fprintf(f,"\n\n");
	}
		/* Summary statistics table: */
	if (traits == 2) {
		if (Ne == 0) fprintf(f,"Step\t\tsig/Ne\teps\tphi\t\tlnL\tLR\tdf\tp\n");
		else fprintf(f,"Step\t\tsig\teps\tphi\t\tlnL\tLR\tdf\tp\n");
		fprintf(f,"------------------------------------------------------------------------------\n");
		if (Gmat[0][0] != 0) {
			if (Ne > 0) {
				fprintf(f,"Direct\t\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[0][0][1],ell[0][1][1],acos(evecs[0][0][0]));
				fprintf(f,"%4.2lf\t-\t-\t-\n",like[0]);
				tmpLR = 2*(like[1] - like[0]);
				fprintf(f,"1.Size\t\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[1][0][1],ell[1][1][1],acos(evecs[0][0][0]));
				fprintf(f,"%4.2lf\t%4.2lf\t1\t%4.4lf\n",like[1],tmpLR,chisquare(tmpLR,1));
			}
			else {
				fprintf(f,"1.Size\t\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[1][0][1],ell[1][1][1],acos(evecs[0][0][0]));
				fprintf(f,"%4.2lf\t-\t-\t-\n",like[1]);
			}
			tmpLR = 2*(like[2] - like[1]);
			fprintf(f,"2.Shape\t\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[2][0][1],ell[2][1][1],acos(evecs[0][0][0]));
			fprintf(f,"%4.2lf\t%4.2lf\t1\t%4.4lf\n",like[2],tmpLR,chisquare(tmpLR,1));
			tmpLR = 2*(like[3] - like[2]);
			fprintf(f,"3.Orientation\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[3][0][1],ell[3][1][1],acos(evecs[1][0][0]));
			fprintf(f,"%4.2lf\t%4.2lf\t1\t%4.4lf\n",like[3],tmpLR,chisquare(tmpLR,1));
		}
		else if (Gmat[0][0] == 0) {
			fprintf(f,"ML best fit\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[3][0][1],ell[3][1][1],acos(evecs[1][0][0]));
			fprintf(f,"%4.2lf\n",like[3]);
		}
		tmpLR = 2*(like[3] - like[4]);
		fprintf(f,"Test vector\t%4.4lf\t%4.4lf\t%4.4lf\t\t",ell[4][0][1],ell[4][1][1],acos(evecs[2][0][0]));
		fprintf(f,"%4.2lf\t%4.2lf\t1\t%4.4lf\n\n\n",like[4],tmpLR,chisquare(tmpLR,1));
	}
	else if (traits > 2) {
		fprintf(f,"Step\t\tlnL\tLR\tdf\tp\n");
		fprintf(f,"--------------------------------------------------------\n");
		if (Gmat[0][0] != 0) {
			if (Ne > 0) {
				fprintf(f,"Direct\t\t%4.2lf\t-\t-\t-\n",like[0]);
				tmpLR = 2*(like[1] - like[0]);
				fprintf(f,"1.Size\t\t%4.2lf\t%4.2lf\t1\t%4.4lf\n",like[1],tmpLR,chisquare(tmpLR,1));
			}
			else fprintf(f,"1.Size\t\t%4.2lf\t-\t-\t-\n",like[1]);
			for(i=1; i<traits; i++) {
				tmpLR = 2*(like[i+1] - like[i]);
				fprintf(f,"%d.Shape(%d)\t%4.2lf\t%4.2lf\t1\t%4.4lf\n",(i+1),i,like[i+1],tmpLR,chisquare(tmpLR,1));
			}
			for(i=1; i<traits; i++) {
				tmpLR = 2*(like[traits+i] - like[traits+i-1]);
				fprintf(f,"%d.Orient.(%d)\t%4.2lf\t%4.2lf\t%d\t%4.4lf\n",(traits+i),i,like[traits+i],tmpLR,i,chisquare(tmpLR,i));
			}
		}
		else {
			fprintf(f,"ML best fit\t%4.2lf\n",like[2*traits-1]);
		}
		tmpLR = 2*(like[2*traits - 1] - like[2*traits]);
		fprintf(f,"Test vector\t%4.2lf\t%4.2lf\t%d\t%4.4lf\n\n\n",like[2*traits],tmpLR,(traits-1),chisquare(tmpLR,(traits-1)));
	}
	fprintf(f,"Eigenvalues (eigenvectors) for G-matrix estimates:\n");
	fprintf(f,"--------------------------------------------------\n");
	if (Gmat[0][0] != 0) {
		if (Ne > 0) {
			fprintf(f,"Direct:\nSigma = %4.4lf\n", ell[0][0][1]);
			tmpLR = 0;
			for(i=0; i<(traits-1); i++) {
				fprintf(f,"%4.4lf\t(",ell[0][i+1][1]);
				tmpLR += ell[0][i+1][1];
				for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[0][j][i]);
				fprintf(f,"%4.4lf)\n",evecs[0][traits-1][i]);
			}
			fprintf(f,"%4.4lf\t(", (1 - tmpLR));
			for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[0][j][traits-1]);
			fprintf(f,"%4.4lf)\n\n",evecs[0][traits-1][traits-1]);
		}
		fprintf(f,"1.Size (with Ne = %4.2lf):\nSigma = %4.4lf [%4.4lf to %4.4lf]\n", popsize, ell[1][0][1],ell[1][0][0],ell[1][0][2]);
		tmpLR = 0;
		for(i=0; i<(traits-1); i++) {
			fprintf(f,"%4.4lf\t(",ell[1][i+1][1]);
			tmpLR += ell[1][i+1][1];
			for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[0][j][i]);
			fprintf(f,"%4.4lf)\n",evecs[0][traits-1][i]);
		}
		fprintf(f,"%4.4lf\t(", (1 - tmpLR));
		for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[0][j][traits-1]);
		fprintf(f,"%4.4lf)\n\n",evecs[0][traits-1][traits-1]);
		for(k=2; k<(traits+1); k++) {
			fprintf(f,"%d.Shape(%d):\nSigma = %4.4lf [%4.4lf to %4.4lf]\n",k,(k-1), ell[k][0][1],ell[k][0][0],ell[k][0][2]);
			tmpLR = 0;
			for(i=0; i<(traits-1); i++) {
				fprintf(f,"%4.4lf\t(",ell[k][i+1][1]);
				tmpLR += ell[k][i+1][1];
				for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[0][j][i]);
				fprintf(f,"%4.4lf)\n",evecs[0][traits-1][i]);
			}
			if (traits == 2) fprintf(f,"[%4.4lf to %4.4lf]\n", ell[k][1][0], ell[k][1][2]);
			fprintf(f,"%4.4lf\t(", (1 - tmpLR));
			for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[0][j][traits-1]);
			fprintf(f,"%4.4lf)\n\n",evecs[0][traits-1][traits-1]);
		}
		for(k=1; k<traits; k++) {
			fprintf(f,"%d.Orientation(%d):\nSigma = %4.4lf [%4.4lf to %4.4lf]\n",(traits+k),k,ell[traits+k][0][1],ell[traits+k][0][0],ell[traits+k][0][2]);
			tmpLR = 0;
			for(i=0; i<(traits-1); i++) {
				fprintf(f,"%4.4lf\t(",ell[traits+k][i+1][1]);
				tmpLR += ell[traits+k][i+1][1];
				for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[k][j][i]);
				fprintf(f,"%4.4lf)\n",evecs[k][traits-1][i]);
			}
			if (traits == 2) fprintf(f,"[%4.4lf to %4.4lf]\n", ell[traits+k][1][0], ell[traits+k][1][2]);
			fprintf(f,"%4.4lf\t(", (1 - tmpLR));
			for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[k][j][traits-1]);
			fprintf(f,"%4.4lf)\n\n",evecs[k][traits-1][traits-1]);
		}
	}
	else {
		fprintf(f,"ML best fit:\nSigma = %4.4lf [%4.4lf to %4.4lf]\n",ell[2*traits-1][0][1],ell[2*traits-1][0][0],ell[2*traits-1][0][2]);
		tmpLR = 0;
		for(i=0; i<(traits-1); i++) {
			fprintf(f,"%4.4lf\t(",ell[2*traits-1][i+1][1]);
			tmpLR += ell[2*traits-1][i+1][1];
			for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[traits-1][j][i]);
			fprintf(f,"%4.4lf)\n",evecs[traits-1][traits-1][i]);
		}
		if (traits == 2) fprintf(f,"[%4.4lf to %4.4lf]\n", ell[2*traits-1][1][0], ell[2*traits-1][1][2]);
		fprintf(f,"%4.4lf\t(", (1 - tmpLR));
		for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[traits-1][j][traits-1]);
		fprintf(f,"%4.4lf)\n\n",evecs[traits-1][traits-1][traits-1]);
	
	}
	fprintf(f,"Test vector:\nSigma = %4.4lf [%4.4lf to %4.4lf]\n",ell[2*traits][0][1],ell[2*traits][0][0],ell[2*traits][0][2]);
	tmpLR = 0;
	for(i=0; i<(traits-1); i++) {
		fprintf(f,"%4.4lf\t(",ell[2*traits][i+1][1]);
		tmpLR += ell[2*traits][i+1][1];
		for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[traits][j][i]);
		fprintf(f,"%4.4lf)\n",evecs[traits][traits-1][i]);
	}
	if (traits == 2) fprintf(f,"[%4.4lf to %4.4lf]\n", ell[2*traits][1][0], ell[2*traits][1][2]);
	fprintf(f,"%4.4lf\t(", (1 - tmpLR));
	for(j=0; j<(traits-1); j++) fprintf(f,"%4.4lf , ",evecs[traits][j][traits-1]);
	fprintf(f,"%4.4lf)\n\n",evecs[traits][traits-1][traits-1]);
	fclose (f);
	printf ("\nDone -- check Outfile.txt\n\n");
	for (i=0; i<(2*traits+1); i++) {
		for (j=0; j<traits; j++) free(ell[i][j]);
		free(ell[i]);
	}
	for (i=0; i<(traits+1); i++) {
		for (j=0; j<(traits+1); j++) free(evecs[i][j]);
		free(evecs[i]);
	}
	free(like); free(ell); free(evecs);
}

double lnlikelihood(double **Inv, double N, bn deter) {
	/* Returns log likelihood for input parameters and data */
	double out, numer, denom;
	numer = vecmatvec(Inv);
	numer *= -0.5*N;	
	denom = 0.5*(deter.d + sizeA*log(2*M_PI/N));
	out = numer - denom;
	return out;
}

double maxNe(double **Inv) {
	/* Returns MLE for Ne */
	double out;
	out = sizeA/vecmatvec(Inv);
	return out;
}

double confN(double **Inv, double N, bn deter, int d) {
	int i;	/* Newton-Raphson method to find Ne at 1.92 likelihood units from N
		in direction d, using explicit derivatives of likelihood function */
	double vmv, lr, adjmax, tryN, epsilon;
	vmv = vecmatvec(Inv);
	adjmax = 0.5*(sizeA*log(N) - N*vmv) - 1.92;
	tryN = N + d*(N/2);
	i=0; lr = 5;
	while (i< 100 && (lr < (-1*pow(0.1,acc)) || lr > (pow(0.1,acc)))) {
		lr = adjmax - 0.5*(sizeA*log(tryN) - tryN*vmv);
		epsilon = 2*lr/(sizeA/tryN - vmv);
		if((tryN + epsilon) < 0) tryN *= 0.1;
		else if(d*epsilon < d*(N - tryN)) tryN = N + d*(N - tryN);
		else tryN += epsilon;
		if (tryN < 0.001) return 0;
		i++;
	}
	return tryN;
}

double est_sigma(double **ellipse, double **evecs, double N, int y) {
	int i, j;  /* if y == 0, returns log likelihood for G-matrix parameter set in ellipse and evecs
		if y==1, re-estimates sigma and populates ellipse[0][1] 
		if y==2, also estimates confidence limits for sigma into ellipse[0] */
	bn *deter;
	double **tmpG, **tmpinvAmat, maxN, temp, like;
	if ((deter = malloc(sizeof(bn))) == NULL) memerror();
	if ((tmpG = malloc(traits*sizeof(double *))) == NULL) memerror();
	if ((tmpinvAmat = malloc(sizeA*sizeof(double *))) == NULL) memerror();
	for (i=0; i<traits; i++) if ((tmpG[i] = malloc(traits*sizeof(double))) == NULL) memerror();
	for (i=0; i<sizeA; i++) if ((tmpinvAmat[i] = malloc(sizeA*sizeof(double))) == NULL) memerror();
	j = partomat(tmpG, ellipse, evecs, traits);
	j += populate(tmpG, tmpinvAmat, deter);
	like = lnlikelihood(tmpinvAmat, N, deter[0]);
	if (y > 0) {  /* Re-estimates sigma based on sigma/Ne functioning as a single parameter */
		maxN = maxNe(tmpinvAmat);
		like = lnlikelihood(tmpinvAmat, maxN, deter[0]);
		ellipse[0][1] *= N/maxN;
		if (y == 2) {
			temp = confN(tmpinvAmat, maxN, deter[0], 1);
			ellipse[0][0] = ellipse[0][1]*maxN/temp;
			temp = confN(tmpinvAmat, maxN, deter[0], -1);
			ellipse[0][2] = ellipse[0][1]*maxN/temp;
		}
	}
	for(i=0; i< traits; i++) free(tmpG[i]);
	for(i=0; i<sizeA; i++) free(tmpinvAmat[i]);
	free(tmpG); free(tmpinvAmat); free(deter);
	if (j == 0) return like;
	else return 1;
}

double est_epsilon(double **ellipse, double **evecs, double N, int which, int y) {
	/* starts with ellipse, estimates sigma, epsilon[which] and above
	* populates ellipse[0][1], ellipse[which and above][1]
	* doesn't touch evecs or ellipse[1][1] through ellipse[which-1][1]
	* y == 0 produce most accurate, y == 1 quicker and less accurate
	* if y == 0 and traits == 2 calculates marginal confidence limits, populates ellipse[1][0] and [1][2]
	* returns likelihood of final model */
	/* NOTE: if any epsilon == 0, the Amat matrix is not invertible */
	int i, j, k;  
	double *like, *dx, err, temp, maxtotal, maxhere;
	if ((like = malloc(3*sizeof(double))) == NULL) memerror();
	if ((dx = malloc(traits*sizeof(double))) == NULL) memerror();
	temp = 0;
	for (i=1; i<which; i++) temp += ellipse[i][1];
	maxtotal = (1 - temp);
	ellipse[0][1] = (double) traits;
	for (i=which; i<traits; i++) {
		ellipse[i][1] = maxtotal/(traits - which + 1); /* arbitrary seed */
		dx[i] = 0.5*ellipse[i][1];
	}
	like[1] = est_sigma(ellipse, evecs, N, 1);
	k = 0; err = 1;
	while (err > pow(0.1,(acc - y)) && k < pow(5,(acc - y))) { /* estimation */
		for (i=which; i<traits; i++) {
			temp = 0;
			for (j=which; j<traits; j++) if(i != j) temp += ellipse[j][1];
			maxhere = maxtotal - temp;
			dx[i] = minimum(0.5*ellipse[i][1], dx[i]);
			dx[i] = minimum(0.5*(maxhere - ellipse[i][1]), dx[i]);
			ellipse[i][1] -= dx[i]; 
			like[0] = est_sigma(ellipse, evecs, N, 0);
			ellipse[i][1] += 2*dx[i]; 
			like[2] = est_sigma(ellipse, evecs, N, 0);
			ellipse[i][1] -= dx[i]; 
			like[1] = est_sigma(ellipse, evecs, N, 0);
			err = 0;
			if ((temp = like[0] + like[2] - 2*like[1]) < 0) {
				ellipse[i][1] -= 0.5*dx[i]*(like[2] - like[0])/temp;
				if (ellipse[i][1] <= 0) ellipse[i][1] = 0.6*dx[i];
				if (ellipse[i][1] >= maxhere) ellipse[i][1] = maxhere - 0.6*dx[i];
				err = maximum(-1*temp, err);
				dx[i] *= 0.5;  
			}
			else if (like[2] > like[0]) {ellipse[i][1] += dx[i]; err = 1; }
			else if (like[2] < like[0]) {ellipse[i][1] -= dx[i]; err = 1; }
		}
		like[1] = est_sigma(ellipse, evecs, N, 1);
		k++;
	}
	if (traits == 2 && y == 0) {
		temp = ellipse[1][1];
		for (j=-1; j<2; j += 2) {
			if (j == -1) dx[0] = 0.5*temp; else dx[0] = 0.5*(1 - temp);
			err = -1; k = 0;
			while (fabs(err) > pow(0.1,acc) && k < pow(5,acc)) {
				if (err > 0) ellipse[1][1] -= j*dx[0];
				else ellipse[1][1] += j*dx[0];
				like[0] = est_sigma(ellipse, evecs, N, 0);
				err = like[1] - like[0] - 1.92;
				dx[0] *= 0.5; 
				k++;
			}
			ellipse[1][j + 1] = ellipse[1][1];
			ellipse[1][1] = temp;
		}
	}
	temp = est_sigma(ellipse, evecs, N, 2);
	free(like); free(dx);
	return temp;
}

double est_phi(double **ellipse, double **evecs, double N, int which) {
	/* Estimates best fit for phi[which] and above, epsilons, sigma, with conf limits for epsilon; 
	populates ellipse and evecs; returns likelihood model */
	int i, j, k, l, m, yet;
	double dx, phi, temp, *like, **tmpevecs;
	if ((like = malloc(3*sizeof(double))) == NULL) memerror();
	if ((tmpevecs = malloc((traits+1)*sizeof(double *))) == NULL) memerror();
	for (i=0; i<(traits+1); i++) if ((tmpevecs[i] = malloc(traits*sizeof(double))) == NULL) memerror();
	phi = 0.1*M_PI;
	for (k=which; k<traits; k++) {
		if (k == (traits-1)) l = which;
		else l = k + 1;
		for (i=0; i<traits; i++) 
			for (j=0; j<traits; j++)
				tmpevecs[i][j] = evecs[i][j];
		like[0] = est_epsilon(ellipse, tmpevecs, N, 1, 1);
		for (i=0; i<4; i++) { /* rotate through 0 to 0.5pi at 0.1pi for each pair; seed process with best */
			rotate(tmpevecs, phi, k, l, traits);
			like[1] = est_epsilon(ellipse, tmpevecs, N, 1, 1);
			if (like[1] > like[0]) {
				for (j=0; j<traits; j++) {
					evecs[j][k] = tmpevecs[j][k];
					evecs[j][l] = tmpevecs[j][l];
				}
				like[0] = like[1];
			}	
		}
	}	
	m = 0; yet = 0;
	while (yet < (traits - which) && m < pow(5,acc)) {  /* loop until all evecs are within accuracy */
		yet = 0; dx = pow(0.5,m)*0.025*M_PI;
		for (k=which; k<traits; k++) {  /* loop through all evec pairs */
			if (k == (traits-1)) l = which;
			else l = k + 1;
			temp = 1; like[0] = 1; like[1] = 0; like[2] = 1;
			while (temp > 0 && like[0] > like[1] && like[2] > like[1]) {
				like[1] = est_epsilon(ellipse, evecs, N, 1, 0);  /* calculate second derivatives wrt phi */
				for (i=0; i<traits; i++) 
					for (j=0; j<traits; j++)
						tmpevecs[i][j] = evecs[i][j];
				rotate (tmpevecs, -1*dx, k, l, traits);
				like[0] = est_epsilon(ellipse, tmpevecs, N, 1, 1);
				rotate (tmpevecs, 2*dx, k, l, traits);
				like[2] = est_epsilon(ellipse, tmpevecs, N, 1, 1);
				if((temp = like[0] + like[2] - 2*like[1]) < 0)  /* rotate acc. to second derivative */
					phi = dx - 0.5*dx*(like[2] - like[0])/temp;
				else if (like[2] > like[0]) phi = dx;
				else if (like[2] < like[0]) phi = -dx;
				else phi = 0;
				rotate (evecs, phi, k, l, traits);
			}
			if ((like[1] - like[0]) < pow(0.1,acc) && (like[1] - like[2]) < pow(0.1,acc)) yet += 1;
		}
		m++;
	}
	if (which == 0) {  /* If which==0, re-order eigenvectors according to eigenvalues */
		temp = 0;
		for(i=0; i<(traits-1); i++) {
			evecs[traits][i] = ellipse[i+1][1];
			temp += ellipse[i+1][1];
		}
		evecs[traits][traits-1] = (1 - temp);
		j = ordereigen(evecs, traits);
		for(i=0; i<(traits-1); i++) ellipse[i+1][1] = evecs[traits][i];
	}
	temp = est_epsilon(ellipse, evecs, N, 1, 0);
	for (i=0; i<(traits+1); i++) free(tmpevecs[i]);
	free (tmpevecs); free (like);
	return temp;
}

int rotate(double **evecs, double phi, int k, int l, int n) {
	int i;   /* Rotates eigenvectors k and l by angle phi, doesn't change eigenvalues */
	double c, s, **tmpevecs;
	if ((tmpevecs = malloc(n*sizeof(double *))) == NULL) memerror();
	for (i=0; i<n; i++) if ((tmpevecs[i] = malloc(2*sizeof(double))) == NULL) memerror();
	c = cos(phi); s = sin(phi);
	for (i=0; i<n; i++) {
		tmpevecs[i][0] = c*evecs[i][k] + s*evecs[i][l];
		tmpevecs[i][1] = c*evecs[i][l] - s*evecs[i][k];
	}
	for (i=0; i<n; i++) {
		evecs[i][k] = tmpevecs[i][0];
		evecs[i][l] = tmpevecs[i][1];
	}
	for (i=0; i<n; i++) free(tmpevecs[i]);
	free(tmpevecs);
	return 0;
}

double vecmatvec(double **Inv) {
	/* Returns the vector*matrix*vector product:
	 (obs - exp)*(inverse covariance)*(obs - exp) */
	int i, j;
	double out, temp, *vector;
	if ((vector = (double *)malloc(sizeA*sizeof(double))) == NULL) memerror();
	for(i=0; i<spp; i++) {
		for(j=0; j<traits; j++)
			vector[i*traits + j] = taxa[i].trait[j] - mean[j];
	}
	out = 0;
	for(i=0; i<sizeA; i++) {
		temp = 0;
		for(j=0; j<sizeA; j++)
			temp += vector[j]*Inv[j][i];
		out += temp*vector[i];
	}
	free (vector);
	return out;
}

int populate(double **G, double **Inv, bn *deter) {
	int i,j;
	double **invG;
	bn tmp;
	if ((invG = malloc(traits*sizeof(double *))) == NULL) memerror();
	for (i=0; i<traits; i++) if ((invG[i] = malloc(traits*sizeof(double))) == NULL) memerror();
	i = 0;
	i += LUinverse(G, invG, traits);
	tmp = logdeterminant(G, traits);
	deter[0].i = pow(Treedeter.i,traits)*pow(tmp.i,spp);
	deter[0].d = traits*Treedeter.d + spp*tmp.d;
	i += SquareKroneckerProduct(invTree, invG, Inv, spp, traits);
	for(j=0; j<traits; j++) free(invG[j]);
	free(invG);
	if (i > 0) return 1;
	else return 0;
}

int mattopar(double **mat, double **ellipse, double **evecs, int n) {
	int i, zero; double tmp;
	zero = eigen(mat, evecs, n);
	tmp = 0;
	for(i=0; i<n; i++) tmp += evecs[n][i];
	ellipse[0][1] = tmp;
	for(i=0; i<(n-1); i++)
		ellipse[i+1][1] = evecs[n][i]/tmp;
	return zero;
}

int partomat(double **mat, double **ellipse, double **evecs, int n) {
	int i, j, zero; double tmp, **tmpevecs;
	if ((tmpevecs = malloc((n+1)*sizeof(double *))) == NULL) memerror();
	for (i=0; i<(n+1); i++) if ((tmpevecs[i] = malloc(n*sizeof(double))) == NULL) memerror();
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			tmpevecs[i][j] = evecs[i][j];
	for(i=0; i<(n-1); i++)
		tmpevecs[n][i] = ellipse[i+1][1]*ellipse[0][1];
	tmp = 0;
	for(i=0; i<(n-1); i++) tmp += tmpevecs[n][i];
	tmpevecs[n][n-1] = ellipse[0][1] - tmp;
	zero = inveigen(tmpevecs, mat, n, 0);
	for(i=0; i<(n+1); i++) free(tmpevecs[i]);
	free(tmpevecs);
	return zero;
}

void memerror() {
	printf("\n\n** Memory Error **\n\n\n\n");
	exit(1);
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 0.00001

#define TS 0.00005

void print(const double *x, const int n)
{
	int index;
	for (index = 0; index < n; ++index)
	{
		printf("%5.4f ", x[index]);
	}
	printf("\n\n");
}

void print_matrix(const double **x, const int n, const int m)
{
	int index1, index2;
	for (index1 = 0; index1 < n; ++index1)
	{
		for (index2 = 0; index2 < m; ++index2)
		{
			printf("%5.4f ", x[index1][index2]);
		}
		printf("\n");
	}
	printf("\n\n");
}

double* diff(const double *x, const int n)
{
	int index;
	double *temp = (double*)malloc((n - 1)*sizeof(double));
	for (index = 0; index < n - 1; ++index)
	{
		temp[index] = x[index + 1] - x[index];
	}
	return temp;
}

double* allocate(int n)
{
	double *temp = (double*)calloc(n, sizeof(double));
	return temp;
}

double** matrix(int n, int m)
{
	double **temp = (double**)malloc(n*sizeof(double*));
	int index;
	for (index = 0; index < n; ++index)
	{
		temp[index] = (double*)calloc(m, sizeof(double));
	}
	return temp;
}

double* solve_tridiagonal_equation(double **x, double *y, int n)
{
	int index;
	double *eredmeny = allocate(n);
	double temp;
	// elso sorra vetitetese a masodikra
	temp = x[1][0] / x[0][0];
	x[1][0] = 0;
	x[1][1] -= temp*x[0][1];
	x[1][2] -= temp*x[0][2];
	y[1] -= temp*y[0];

	// kozepso elemek
	for (index = 2; index < n - 1; ++index)
	{
		temp = x[index][0] / x[index - 1][1];
		x[index][0] = 0;
		x[index][1] -= temp*x[index - 1][2];
		y[index] -= temp*y[index - 1];
	}

	// utolso elem
	temp = x[n - 1][0] / x[n - 3][1];
	x[n - 1][0] = 0;
	x[n - 1][1] -= temp*x[n - 2][2];
	y[n - 1] -= temp*y[n - 3];

	temp = x[n - 1][1] / x[n - 2][1];
	x[n - 1][1] = 0;
	x[n - 1][2] -= temp*x[n - 2][2];
	y[n - 1] -= temp*y[n - 2];

	// c ertekeinek kiszamitasa
	eredmeny[n - 1] = y[n - 1] / x[n - 1][2];
	for (index = n - 2; index > 0; --index)
	{
		eredmeny[index] = (y[index] - (x[index][2] * eredmeny[index + 1])) / x[index][1];
	}
	eredmeny[0] = (y[0] - (x[0][1] * eredmeny[index + 1]) - (x[0][2] * eredmeny[index + 2])) / x[index][0];
	//print(eredmeny, n);
	return eredmeny;
}

int floatcmp(const void *elem1, const void *elem2)
{
	if (*(const double*)elem1 < *(const double*)elem2)
		return -1;
	else
		return *(const double*)elem1 > *(const double*)elem2;
}

void evaluate_spline(const double *x, const double *a, const double *b, const double *c, const double *d, const int n, const double *time, double *yout, int nt)
{
	int index;
	int klo, khi, compute_index;

	// rendezzuk a time tombben levo ertekeket, biztos ami biztos
	qsort(time, nt, sizeof(double), floatcmp);
	klo = 0; khi = n - 1;
	while (khi - klo > 1)
	{
		compute_index = (khi + klo) >> 1;
		if (x[compute_index] > time[0]) khi = compute_index;
		else klo = compute_index;
	}
	compute_index = klo;
	double xi_x;
	for (index = 0; index < nt; ++index)
	{
		if ((compute_index + 1 < n) && (time[index] + EPS >= x[compute_index + 1]))
		{
			while ((compute_index + 1 < n) && (time[index] + EPS >= x[compute_index + 1]))
				++compute_index;
		}
		// kiszamolni, behejettesiteni: S(index) = a(index) + b(index(x(c_i) - time(index)) + ...
		xi_x = time[index] - x[compute_index];
		yout[index] = a[compute_index] + b[compute_index] * xi_x + c[compute_index] * xi_x*xi_x +
			d[compute_index] * xi_x*xi_x*xi_x;
	}
}

void spline(const double *x, double *y, int n, double *time, double *yout, int nt)
{
	int index;

	if (yout == NULL)
		yout = allocate(nt);

	if (n == 2) // ha 2 elemu a bemenet, egy egyszeru egyenes a kimenet
	{
		for (index = 0; index < nt; ++index)
		{
			yout[index] = ((y[1] - y[0]) / (x[1] - x[0]))*(time[index] - x[0]) + y[0];
		}
	}
	else
	{
		double *a = y;
		double *b = allocate(n);
		double *c = allocate(n);
		double *d = allocate(n);
		if (n == 3)
		{
			// natural cubic spline -> kezzel megoldhato egyenletek
			// printf("Erre meg nincs megoldas!\n\n");
			double h1 = x[1] - x[0];
			double h2 = x[2] - x[1];
			double dy1 = y[1] - y[0];
			double dy2 = y[2] - y[1];

			c[1] = (3 * dy2 / h2 - 3 * dy1 - h1) / (2 * (x[0] + x[2]));
			d[0] = c[1] / (3 * h1); d[1] = -c[1] / (3 * h2);
			b[0] = dy1 / h1 - d[0] * h1*h1;
			b[1] = dy2 / h2 - c[1] * h2 - d[1] * h2*h2;
			evaluate_spline(x, a, b, c, d, n - 1, time, yout, nt);
		}
		else
		{
			// S(x) = aj + bj(x-xj) + cj(x-xj)^2 + dj(x-xj)^3

			// aj = yj
			//memcpy(a, y, n*sizeof(float));
			double *h = diff(x, n);
			double *dy = diff(y, n);
			double *my = allocate(n);

			double **mx = matrix(n, 3);

			mx[0][0] = h[1]; mx[0][1] = -(h[0] + h[1]); mx[0][2] = h[0];
			mx[n - 1][0] = h[n - 2];
			mx[n - 1][1] = -(h[n - 3] + h[n - 2]);
			mx[n - 1][2] = h[n - 3];
			for (index = 1; index < n - 1; ++index)
			{
				mx[index][0] = h[index - 1];
				mx[index][1] = 2 * (h[index - 1] + h[index]);
				mx[index][2] = h[index];
				my[index] = 3 * dy[index] / h[index] - 3 * dy[index - 1] / h[index - 1];
			}
			// print_matrix(c, n, 3);
			// print(my, n);
			c = solve_tridiagonal_equation(mx, my, n);
			free(my); //free(mx) -> it's a matrix
			for (index = 0; index < n; ++index)
			{
				d[index] = (c[index + 1] - c[index]) / (3 * h[index]);
				b[index] = dy[index] / h[index] - c[index] * h[index] - d[index] * h[index] * h[index];
			}

			evaluate_spline(x, a, b, c, d, n - 1, time, yout, nt);
			free(b); free(c); free(d); free(h); free(dy);
			/*print(a, n);
			print(b,n);
			print(c, n);
			print(d, n);*/
		}
	}
}

int find_max_peaks(const double *x, const double *y, int n, double **max_x, double **max_y)
{
	int index;
	int counter = 0;
	
	*max_x = allocate(n);
	*max_y = allocate(n);

	for (index = 1; index < n - 1; ++index)
	{
		if (y[index] > y[index - 1] && y[index] > y[index + 1]){
			(*max_x)[counter] = index*TS;
			(*max_y)[counter] = y[index];
			++counter;
		}
	}
	//max_x = (double*)realloc(max_x, counter*sizeof(double));
	//max_y = (double*)realloc(max_y, counter*sizeof(double));
	return counter;
}

int find_min_peaks(const double *x, const double *y, int n, double *min_x, double *min_y)
{
	int index;
	int counter = 0;

	min_x = allocate(n);
	min_y = allocate(n);

	for (index = 1; index < n - 1; ++index)
	{
		if (y[index] > y[index - 1] && y[index] > y[index - 1]){
			min_x[counter] = index*TS;
			min_y[counter] = y[index];
		}
	}
	min_x = (double*)realloc(min_x, counter*sizeof(double));
	min_y = (double*)realloc(min_y, counter*sizeof(double));
	return counter;
}

void emd(double *x, double *y, int n)
{

}

int main()
{
	int index;
	/*	float *x = (float*)malloc(N*sizeof(float));
	float *y = (float*)malloc(N*sizeof(float));
	for(index=0; index < N; ++index)
	{
	x[index] = index * (PI/10);
	y[index] = sin(x[index]);
	}*/

	/*float *x = allocate(2);
	float *y = allocate(2);
	float *time = allocate(2);
	float *yout = allocate(2);
	x[0] = 1; x[1] = 2;
	y[0] = 2; y[1] = 3;
	time[0] = 0; time[1] = 4;
	spline(x, y, 2, time, yout, 2);*/

	//double *x = allocate(5);
	//double *y = allocate(5);
	// x[0] = 1; x[1] = 2; x[2] = 3; x[3] = 4; x[4] = 5;
	// y[0] = 2; y[1] = 3; y[2] = 4; y[3] = 5; y[4] = 6;

	/*x[0] = 0.1; x[1] = 0.4; x[2] = 1.2; x[3] = 1.8; x[4] = 2.0;
	y[0] = 0.1; y[1] = 0.7; y[2] = 0.6; y[3] = 1.1; y[4] = 0.9;

	double *time = allocate(61);
	double *yout = allocate(61);
	double val;
	for (index = 0, val = -0.5; val <= 2.5; val += TS, ++index)
	{
		time[index] = val;
	}
	//time[0] = 0.1; time[1] = 0.4; time[2] = 1.2; time[3] = 1.8; time[4] = 2;
	spline(x, y, 5, time, yout, 61);

	/*FILE *f = fopen("data.txt", "w");
	for (index = 0; index < 61; ++index)
		fprintf(f, "%5.5f ", yout[index]);*/
	
	// free(x); free(y); free(time); free(yout);

	double *x = allocate(99946);
	double *y = allocate(99946);
	double *max_x = allocate(99946);
	double *max_y = allocate(99946);
	double *time = allocate(99946);
	double *yout = allocate(99946);

	int n = 99946;
	double val;
	int i;
	FILE *fin = fopen("jel.txt", "r");
	fscanf(fin, "%lf", &val);
	fscanf(fin, "%lf", &val);

	for (i = 0; i < 99946; ++i){
		fscanf(fin, "%lf", &y[i]);
		time[i] = i*TS;
	}
	fclose(fin);
	
	int counter = find_max_peaks(x, y, n, &max_x, &max_y);
	spline(max_x, max_y, i, time, yout, 99946);

	//printf("%5.6lf %5.6lf\n", time[0], yout[0]);

	fin = fopen("data.txt", "w");

	for (i = 0; i < 99946; ++i)
	{
		fprintf(fin, "%5.6lf ", yout[i]);
	}
	fclose(fin);

	printf("Nr. of maxima: %d\n", counter);
	getch(); 
	return 0;
}
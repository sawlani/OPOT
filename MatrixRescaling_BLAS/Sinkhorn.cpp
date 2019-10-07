#include <chrono>
#include "cblas.h"
#include<bits/stdc++.h>
using namespace std;
typedef double D;
typedef vector<D> vD;
mt19937 gene(233);
D EPS = 1;
bool _W = false; //warning
int MAX_ITER = 300000;
D DIV = 1000000;
//D DIV = 200000000;
D * onec, * oner;
	
vD mm(int n, int m, int p, vD A, vD B) {
	//printf("%d %d %d %d %d\n", n, m, p, A.size(), B.size());
	assert(A.size() == n * m);
	assert(B.size() == m * p);
	vD C(n * p);
	for(int i(0); i < n; i++) {
		for(int j(0); j < p; j++) {
			for(int k(0); k < m; k++) {
				C[i * m + j] += A[i * m + k] * B[k * p + j];
			}
		}
	}
	return C;
}
vD operator / (vD c, vD d) {
	for(int i(0); i < (int)c.size(); i++)
		c[i] /= d[i];
	return c;
}
vD operator * (D x, vD C) {
	for(D & y : C) y *= x;
	return C;
}
void exp(int len, D * c, D eta, D * A) {
	for(int i(0); i < len; i++) A[i] = exp(eta * c[i]);
}
vD abs(vD a) {
	for(D & x : a) x = abs(x);
	return a;
}
D sum(vD a) {
	D res = 0;
	for(D x : a) res += x;
	return res;
}
vD diag(vD r) {
	int n = r.size();
	vD res(n * n);
	for(int i(0); i < n; i++) res[i * n + i] = r[i];
	return res;
}

void dist(int n, D * a, D * b, D * c) {
	for(int i(0); i < n; i++) {
		c[i] = b[i] - a[i] + a[i] * log(a[i] / b[i]);
	}
}
vD operator - (const vD & a, const vD & b) {
	vD res(a.size());
	for(int i(0); i < (int)a.size(); i++) res[i] = a[i] - b[i];
	return res;
}
vD operator + (const vD & a, const vD & b) {
	vD res(a.size());
	for(int i(0); i < (int)a.size(); i++) res[i] = a[i] + b[i];
	return res;
}
void print(const vD & x) {
	for(D y : x) printf("%lf ", y);
	printf("\n");
}
void print(int n, D * x) {
	for(int i(0); i < n; i++) printf("%f ", x[i]);
	printf("\n");
}
void eval_dist(int r_dim, int c_dim, D* A, D* r, D* c) {
	static D * row_sum = new D[r_dim];
	static D * col_sum = new D[c_dim];
	cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);
	cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);
	static D * row_dif = new D[r_dim];
	memcpy(row_dif, row_sum, sizeof(D) * r_dim);
	cblas_daxpy(r_dim, -1, r, 1, row_dif, 1);
	static D * col_dif = new D[c_dim];
	memcpy(col_dif, col_sum, sizeof(D) * c_dim);
	cblas_daxpy(c_dim, -1, c, 1, col_dif, 1);
	
	
	print(r_dim, row_dif);
	print(c_dim, col_dif);
}
int argmax(vD c) {
	int res = 0;
	D mx = c[0];
	for(int i(1); i < (int)c.size(); i++) {
		if(c[i] > mx) {
			mx = c[i];
			res = i;
		}
	}
	return res;
}
double max_(D* c, int n) {
	D mx = c[0];
	for(int i(1); i < n; i++) {
		if(c[i] > mx) {
			mx = c[i];			
		}
	}
	return mx;
}

vD getDiag(int n, vD x) {
	vD res(n);
	for(int i(0); i < n; i++) 
		res[i] = x[i * n + i];
	return res;
}
vD mx_trunc(vD x, D mx) {
	for(D & y : x) y = min(y, mx);
	return x;
}
D dot(vD x, vD y) {
	int n = x.size();
	D res = 0;
	for(int i(0); i < n; i++) res += x[i] * y[i];
	return res;
}
struct Rec {
	D eta, eps;
};
D capsum = 0;
D _answer, _eta;
chrono::duration<double> _time;
int _iter;
double sinkhorn(int r_dim, int c_dim, D * C, D * r, D * c, double eta, D * A) {

  	D * C0 = new D[r_dim * c_dim];
	memcpy(C0, C, sizeof(D) * r_dim * c_dim);
	//int time0 = clock();
	
	auto start = chrono::high_resolution_clock::now();
	D C_max = max_(C, r_dim*c_dim);
	cblas_daxpy (r_dim * c_dim, -1 + (1. / max(C_max, (D)0.01)), C, 1, C, 1);	
	exp(r_dim * c_dim, C, -eta, A);

	/*D * r_diag = new D[r_dim * r_dim];
	for(int i=0; i<r_dim*r_dim; i++){		 
	  r_diag[i] = 0;
	}
	D * c_diag = new D[c_dim * c_dim];
	for(int i=0; i<c_dim*c_dim; i++){		 
	  c_diag[i] = 0;
	  }
	D *row0 = new D[c_dim];
	for (int i=0; i<c_dim; i++){
	  row0[i] = 0;
	}
	D *col0 = new D[r_dim];
	for (int i=0; i<r_dim; i++){
	  col0[i] = 0;
	}*/
	  
	/*D * A0 = new D[r_dim * c_dim];
	memset(A0, 0, r_dim*c_dim*sizeof(D));
	D * A_temp;*/
	
	D * row_sum = new D[r_dim];
	D * col_sum = new D[c_dim];
	cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);
	
 	cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);
	D * row_dif = new D[r_dim];
	memcpy(row_dif, row_sum, sizeof(D) * r_dim);
	cblas_daxpy(r_dim, -1, r, 1, row_dif, 1);
	D * col_dif = new D[c_dim];
	memcpy(col_dif, col_sum, sizeof(D) * c_dim);
	cblas_daxpy(c_dim, -1, c, 1, col_dif, 1);
		
	D dist_AU = cblas_dasum(r_dim, row_dif, 1) + cblas_dasum(c_dim, col_dif, 1);
	
	int n_iter = 0;
	vD cost_l;
	
	double eps = 10;
	D * r_ratio = new D[r_dim];
	D * c_ratio = new D[c_dim];
	int row = 1;
	
	while(dist_AU > eps) {
	  
	  if(n_iter > MAX_ITER) {
	    //printf("MAX ITER reached with dist %f.\n", dist_AU);
	    break;
	  }	  
		n_iter += 1;
		if(row) {
		  ///cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);
		  		  
		  	for(int i=0; i<r_dim; i++){
			  //r_diag[i*r_dim + i] = r[i] / row_sum[i];
			  //cblas_daxpy(c_dim, r[i] / row_sum[i], &A[i*c_dim], 1, row0, 1);
			  cblas_daxpy(c_dim, r[i] / row_sum[i] - 1, &A[i*c_dim], 1, &A[i*c_dim], 1);			  
			  /*for(int k=0; k<c_dim; k++){
			    A[i * c_dim + k] = row0[k];
			    row0[k] = 0;
			    }*/
			}
			
			//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, r_dim, c_dim, r_dim, 1, r_diag, r_dim, A, c_dim, 0, A0, c_dim);
			/*for(int i=0; i<r_dim; i++){			 
			  r_ratio[i] = r[i] / row_sum[i];  
			  for(int j(0); j < c_dim; j++){
			    A[i * c_dim + j] = A[i * c_dim + j] * r_ratio[i];
			  }
			  }*/
		}else {
		  ///cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);
		  for(int j=0; j<c_dim; j++){
		    ///c_diag[j*c_dim + j] = c[j] / col_sum[j];			  
		    cblas_daxpy(r_dim, c[j] / col_sum[j] - 1, &A[j], c_dim, &A[j], c_dim);
		    /*for(int k=0; k<r_dim; k++){
		      A[ k * c_dim + j] = col0[k];
		      col0[k] = 0;
		      }*/		
		  }
		  //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, r_dim, c_dim, c_dim, 1, A, c_dim, c_diag, c_dim, 0, A0, c_dim);
		  
			/*
			for(int j=0; j<c_dim; j++){
			  c_ratio[j] = c[j] / col_sum[j];
			  for(int i(0); i < r_dim; i++){
			    A[i * c_dim + j] = A[i * c_dim + j] * c_ratio[j];
			  }
			  }*/
		}
		
		/*A_temp = A;
		A = A0;
		A0 = A_temp;
		memset(A0, 0, r_dim*c_dim*sizeof(D));*/
		cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);		  
		cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);
		//row_dif = new D[r_dim];
			memset(row_dif, 0, r_dim);
			memcpy(row_dif, row_sum, sizeof(D) * r_dim);
			//this subtracts r from row_diff
			cblas_daxpy(r_dim, -1, r, 1, row_dif, 1);
			//col_dif = new D[c_dim];
			memset(col_dif, 0, c_dim);
			memcpy(col_dif, col_sum, sizeof(D) * c_dim);
			cblas_daxpy(c_dim, -1, c, 1, col_dif, 1);
			//take the l1 norm
			dist_AU = cblas_dasum(r_dim, row_dif, 1) + cblas_dasum(c_dim, col_dif, 1);
		
			row ^= 1;
	}
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	_time = elapsed;

	D cur_cost = cblas_ddot(c_dim * r_dim, C0, 1, A, 1);// * C_max;
	//printf("final dist_AU: %f\n", dist_AU);
	//printf("%f %d %f\n", 
	
	_answer = cur_cost; _iter = n_iter;
	//_time = (clock() - time0) / (double)CLOCKS_PER_SEC;
	
	printf("%f %d %f\n", _answer, _iter, _time);
	//printf("dist to marginals after iteration: ");
	memcpy(C, C0, sizeof(double) * r_dim * c_dim);
	return cur_cost;
}

void round(int r_dim, int c_dim, D* A, D* r, D* c, D* res) {
  //assert(res == A);
	D * row_sum = new D[r_dim];
	D * col_sum = new D[c_dim];
	cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);
	cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);
	
	D * sr = new D[r_dim];
	D * sc = new D[c_dim];
	for(int i(0); i < r_dim; i++) sr[i] = min((D)1., r[i] / row_sum[i]);
	for(int i(0); i < c_dim; i++) sc[i] = min((D)1., c[i] / col_sum[i]);
	
	for(int i(0); i < r_dim; i++)
		for(int j(0); j < c_dim; j++)
			A[i * c_dim + j] = A[i * c_dim + j] * sr[i] * sc[j];
/*	vD x = r / sumrow(r_dim, c_dim, A);
	x = mx_trunc(x, 1);
	A = mm(r_dim, r_dim, c_dim, diag(x), A);
	vD y = c / sumcol(r_dim, c_dim, A);
	y = mx_trunc(y, 1);
	A = mm(r_dim, c_dim, c_dim, A, diag(y));*/

	cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);
	cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);

	D * row_dif = new D[r_dim];
	memcpy(row_dif, row_sum, sizeof(D) * r_dim);
	cblas_daxpy(r_dim, -1, r, 1, row_dif, 1);
	D * col_dif = new D[c_dim];
	memcpy(col_dif, col_sum, sizeof(D) * c_dim);
	cblas_daxpy(c_dim, -1, c, 1, col_dif, 1);	
	
//	vD err_r = (r - sumrow(r_dim, c_dim, A));
//	vD err_c = (c - sumcol(r_dim, c_dim, A));
	D sa = cblas_dasum(r_dim, row_dif, 1);
	for(int i(0); i < r_dim; i++)
		for(int j(0); j < c_dim; j++)
			A[i * c_dim + j] = A[i * c_dim + j] + row_dif[i] * col_dif[j] / sa;
	//A = A + (1./ sum(abs(err_r))) * mm(r_dim, 1, c_dim, err_r, err_c);
}
int main(int argc, char ** argv) {
	int n, m;
	double std = -1;
	if(argc >= 2) {
		sscanf(argv[1], "%lf", &std);
	}
	while(scanf("%d%d", &n, &m) == 2) {
	  //scanf("%d%d", &n, &m);
	D * r = new D[n], * c = new D[m];
	capsum = 0;
	for(int i(0); i < n; i++) {
		double _r;
		scanf("%lf", &_r);
		r[i] = _r;
		capsum += r[i];
	}
	/*
	for(int i(0); i < n; i++) {
		r[i] /= capsum;
		}*/
	capsum = 0;
	for(int i(0); i < m; i++) {
		double _r;
		scanf("%lf", &_r);
		c[i] = _r;
		capsum += c[i];
	}
	printf("IO done\n");
	fflush(stdout);
	//for(int i(0); i < m ; i++) c[i] /= capsum;

	capsum = 0;
	int r_dim = n, c_dim = m;

	oner = new D[r_dim];
	onec = new D[c_dim];
	for(int i(0); i < r_dim; i++) oner[i] = 1;
	for(int i(0); i < c_dim; i++) onec[i] = 1;
	
	D * C = new D[n * m], * A = new D[n * m];

	for(int i(0); i < r_dim * c_dim; i++) {
		int a;
		scanf("%d", &a);
		C[i] = a;// / DIV;
		//C[i] = i % n == i / n ? 1 : 2;
	}
	
	D eps = 1;
	//Rec opt{4 * log(n) / eps, eps};
	double eta = 50; //3000 for geo
	//printf("eta = %f\n", eta);
	/*for(int i(0); i < r_dim; i++) {
		for(int j(0); j < c_dim; j++) {
			//printf("%lf%c", C[i * c_dim + j], j == c_dim - 1 ? '\n' : ' ');
		}
		}*/
	clock_t time1 = clock();
	/*above: initialization of r_dim, c_dim, r, c, C*/
	//double le = 5368, ri = 5377;
	double le = 0, ri = eta * 2;
	double ans = sinkhorn(r_dim, c_dim, C, r, c, 100, A);
	printf("%f %d %f %f\n", _answer, _iter, _time, _eta);
	exit(0);
	do {
		double mid = (le + ri) / 2;
		//printf("eta = %f\n", mid);
		_eta = mid;
		double ans = sinkhorn(r_dim, c_dim, C, r, c, mid, A);
		if(abs(ans - std) < 0.1 * std) ri = mid;
		else le = mid;
		if(std == -1) break;
	} while((ri - le) > min(1., le * 0.01));
	printf("%f %d %f %f\n", _answer, _iter, _time, _eta);
		//printf("answer after rounding: %f\n", cblas_ddot (n * m, C, 1, A, 1) * DIV);
	//printf("clock() = %f seconds\n", (clock()-time1)/(double)CLOCKS_PER_SEC );
	//std::cout<<"total time: "<< (clock()-time1)/(double)CLOCKS_PER_SEC <<'\n';
	
	//printf("dist to marginals:\n");
	//eval_dist(r_dim, c_dim, A, r, c);

	/*D ans1 = 0;
	static bool vst[5000];
	memset(vst, 0, sizeof(vst));
	for(int i(0); i < n; i++) {
		int ij = -1;
		D mx = -1;
		for(int j(0); j < n; j++) {
			if(A[i * n + j] > mx) {
				mx = A[i * n + j];
				ij = j;
			}
		}
		ans1 += C[i * n + ij];
		//assert(vst[ij] == 0); 
		vst[ij] = 1;
	}
	cout << ans1 << endl;*/
	
	}
}


#include "cblas.h"
#include<bits/stdc++.h>
using namespace std;
typedef double D;
typedef vector<D> vD;
mt19937 gene(233);
D EPS = 1;
bool _W = false; //warning
int MAX_ITER = 10000000; //3000;
D DIV = 1000000;
//D DIV = 200000000;
D * onec, * oner;

D _eps = 0.5, _eta = 50;

void round(int r_dim, int c_dim, D* A, D* r, D* c);



vD mm(int n, int m, int p, vD A, vD B) {
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
	cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, r_dim, onec, 1, 0, row_sum, 1);
	cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, r_dim, oner, 1, 0, col_sum, 1);
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
double min_(D* c, int n) {
	D mx = c[0];
	for(int i(1); i < n; i++) {
		if(c[i] < mx) {
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
int _iter;
D _res, _clock;
D capsum = 0;
double greenkhorn(int r_dim, int c_dim, D * C, D * r, D * c, double eta, D * A) {
	D * C0 = new D[r_dim * c_dim];
	memcpy(C0, C, sizeof(D) * r_dim * c_dim);
	int time0 = clock();
	D C_max = max_(C, r_dim*c_dim);

	D C_min = min_(C, r_dim*c_dim);
	//printf("max C %f   %f\n", C_max, C_min);

	cblas_daxpy (r_dim * c_dim, -1 + (1. / max(C_max, (D)0.01)), C, 1, C, 1);
	exp(r_dim * c_dim, C, -eta, A);

	D maxa0 = -1, mina0 = 1e50;
	for(int i(0); i < r_dim * c_dim; i++) {
		maxa0 = max(maxa0, A[i]);
		mina0 =min(mina0, A[i]);
	}
	//printf("max/min of A0: %f %f\n", maxa0, mina0);
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
	//for(int i(0); i < r_dim; i++) printf("%f ", row_dif[i]); printf("\n");

	D dist_AU = cblas_dasum(r_dim, row_dif, 1) + cblas_dasum(c_dim, col_dif, 1);
	//printf("initial dist_AU: %f\n", dist_AU);
	int n_iter = 0;
	//opt.eps /= 8 * C[cblas_idamax(r_dim * c_dim, C, 1)];

	double eps = _eps;
	D * dist_r = new D[r_dim];
	D * dist_c = new D[c_dim];
	D * r_ratio = new D[r_dim];
	D * c_ratio = new D[c_dim];
	//printf("eps = %f\n", eps);

	while(dist_AU > eps) {
		if(n_iter > MAX_ITER) {
			printf("MAX ITER reached with dist %f.\n", dist_AU);
			break;
		}		
		n_iter += 1;
		dist(r_dim, r, row_sum, dist_r);
		dist(c_dim, c, col_sum, dist_c);
		int mxri = cblas_idamax(r_dim, dist_r, 1);
		int mxci = cblas_idamax(c_dim, dist_c, 1);
		D mxr = dist_r[mxri];
		D mxc = dist_c[mxci];


		if(mxr > mxc) {
			r_ratio[mxri] = r[mxri] / row_sum[mxri];
			for(int j(0); j < c_dim; j++){
				double & aij = A[mxri * c_dim + j];
				row_sum[mxri] -= aij;
				col_sum[j] -= aij;
				aij = aij * r_ratio[mxri];
				row_sum[mxri] += aij;
				col_sum[j] += aij;
			}
		}else {
			c_ratio[mxci] = c[mxci] / col_sum[mxci];
			for(int i(0); i < r_dim; i++){
				double & aij = A[mxci + c_dim * i];
			       	row_sum[i] -= aij;
				col_sum[mxci] -= aij;	
				aij *= c_ratio[mxci];
				row_sum[i] += aij;
				col_sum[mxci] += aij;
			}
		}


		row_dif = new D[r_dim];
		memcpy(row_dif, row_sum, sizeof(D) * r_dim);
		cblas_daxpy(r_dim, -1, r, 1, row_dif, 1);
		col_dif = new D[c_dim];
		memcpy(col_dif, col_sum, sizeof(D) * c_dim);
		cblas_daxpy(c_dim, -1, c, 1, col_dif, 1);
		dist_AU = cblas_dasum(r_dim, row_dif, 1) + cblas_dasum(c_dim, col_dif, 1);

		if (n_iter % 10000 == 0) {
			//D cur_cost = cblas_ddot (c_dim * r_dim, C0, 1, A, 1);		  
			//printf("distAU after iter %d:  %f, mincost = %f\n", n_iter, dist_AU, cur_cost);
		}
		//for(int i = 0; i < r_dim * c_dim; i++) {
		//	assert(isfinite(A[i]));
		//}
	}
	for(int i = 0; i < r_dim * c_dim; i++) assert(isfinite(A[i]));
	//printf("iter = %d\n", n_iter);
	_iter = n_iter;
	double res = 0;
	res = cblas_ddot (r_dim * c_dim, C0, 1, A, 1);
	//round(r_dim, c_dim, A, r, c);
	_res = res;
	_clock = (clock() - time0) / (double)CLOCKS_PER_SEC;
	memcpy(C, C0, sizeof(double) * r_dim * c_dim);
	return res;
	//printf("answer after rounding: %f\n", cblas_ddot (r_dim * c_dim, C0, 1, A, 1));
}

void round(int r_dim, int c_dim, D* A, D* r, D* c) {
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

	cblas_dgemv (CblasRowMajor, CblasNoTrans, r_dim, c_dim, 1, A, c_dim, onec, 1, 0, row_sum, 1);
	cblas_dgemv (CblasRowMajor, CblasTrans, r_dim, c_dim, 1, A, c_dim, oner, 1, 0, col_sum, 1);

	D * row_dif = new D[r_dim];
	memcpy(row_dif, row_sum, sizeof(D) * r_dim);
	cblas_daxpy(r_dim, -1, r, 1, row_dif, 1);
	D * col_dif = new D[c_dim];
	memcpy(col_dif, col_sum, sizeof(D) * c_dim);
	cblas_daxpy(c_dim, -1, c, 1, col_dif, 1);

	D sa = cblas_dasum(r_dim, row_dif, 1);
	for(int i(0); i < r_dim; i++)
		for(int j(0); j < c_dim; j++)
			A[i * c_dim + j] = A[i * c_dim + j] + row_dif[i] * col_dif[j] / sa;
}

int main(int argc, char ** argv) {
	int n, m;
	//printf("%d %s %s %s\n", argc, argv[0], argv[1], argv[2]);
	if(argc >= 2) {
		sscanf(argv[1], "%lf", &_eta);
	//	printf("%f\n", _eta);
	}
	if(argc >= 3) {
		sscanf(argv[2], "%lf", &_eps);
	}
	double std = -1;
	if(argc >= 4) {
		sscanf(argv[3], "%lf", &std);
	}
	D _acc = -1;
	if(argc >= 5) {
		sscanf(argv[4], "%lf", &_acc);
	}
	scanf("%d%d", &n, &m);
	D * r = new D[n], * c = new D[m];
	capsum = 0;
	for(int i(0); i < n; i++) {
		double _r;
		scanf("%lf", &_r);
		r[i] = _r;
		capsum += r[i];
	}
	capsum = 0;
	for(int i(0); i < m; i++) {
		double _r;
		scanf("%lf", &_r);
		c[i] = _r;
		capsum += c[i];
	}

	//auto r_min = min_(r, n);
	//auto c_min = min_(c, m);
	//printf("%f  %f ", r_min, c_min);
	//exit(0);
	//for(int i(0); i < m ; i++) c[i] /= capsum;

	capsum = 0;
	int r_dim = n, c_dim = m;

	oner = new D[r_dim];
	onec = new D[c_dim];
	for(int i(0); i < r_dim; i++) oner[i] = 1;
	for(int i(0); i < c_dim; i++) onec[i] = 1;

	D * C = new D[n * m], * A = new D[n * m];
	int C_min = 0;
	for(int i(0); i < r_dim * c_dim; i++) {
		int a;
		scanf("%d", &a);
		//a /= 1000;
		C[i] = a;// / DIV;
		if (a > 0 && a > C_min) C_min = a;
		//C[i] = i % n == i / n ? 1 : 2;
	}
	//printf("C_min %d \n", C_min);


	double eta = _eta; //8000; //8000 for geo 300 for mnist, 500 for NLP
	//printf("eta = %f\n", eta);
	clock_t time1 = clock();
	/*above: initialization of r_dim, c_dim, r, c, C*/
	double leeta = 0, rieta = 2 * eta;
	do {
		double mid = (leeta + rieta) / 2;
		//printf("mid =%f\n", mid);
		int t1 = clock();
		eta = mid;
		double ans = greenkhorn(r_dim, c_dim, C, r, c, mid, A);
		//printf("time : %f\n", (clock() - t1) / (double)CLOCKS_PER_SEC);
		//printf("%lf %lf\n", abs(ans - std) / std, _acc);
		if(abs(ans - std) / std <= _acc) rieta = mid;
		else leeta = mid;
		if(std == -1) break;

	} while((rieta - leeta) > min(1., leeta * 0.01));
	//round(r_dim, c_dim, A, r, c);
	//printf("answer after rounding: %f\n", cblas_ddot (n * m, C, 1, A, 1));
	//printf("clock() = %f seconds\n", (clock()-time1)/(double)CLOCKS_PER_SEC );
	//std::cout<<"total time: "<< (clock()-time1)/(double)CLOCKS_PER_SEC <<'\n';
	//cout << "======================================================================" << endl;
	printf("%f %f %d %f\n", _res, _clock, _iter, eta);

}


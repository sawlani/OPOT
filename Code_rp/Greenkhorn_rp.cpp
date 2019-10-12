// Greenkhorn that implicitly tracks the dual variables
// (instead of keeping all exponentials around)
//
// eta is aggressiveness parameter, set it to smaller values for faster conv
//
// Code was optimized to take O(n+m) per iteration,
// but should still be an order of magnitude slower than
// hardware optimzied versions
//
//


#include <bits/stdc++.h>
#include <chrono>
using namespace std; 
typedef long long LL;
typedef double D;
typedef vector<int> VI;
typedef set<int> SI;
typedef map<int, int> MII;
typedef pair<int,int> PII;


#define MP make_pair 
#define A first 
#define B second 

#define PB push_back 
#define FR(i, a, b) for(int i=(a); i<(b); i++) 
#define FOR(i, n) FR(i, 0, n) 

template <typename T> inline void SetMin(T &a, T b) {if(b < a) a = b;}
template <typename T> inline void SetMax(T &a, T b) {if(b > a) a = b;}

const int MAXNM = 5000;
int n, m;

double r_cap[MAXNM];
double c_cap[MAXNM];
double cost[MAXNM][MAXNM];

double temp[MAXNM];

double x[MAXNM][MAXNM];
double r_s[MAXNM], c_s[MAXNM];
double y_left[MAXNM], y_right[MAXNM];

double eta = 60;
double epsilon = 3000.0;

double TotalResiduals() {
	double d_r = 0;
	FOR(i, n) d_r += fabs(r_cap[i] - r_s[i]);
	double d_c = 0;
	FOR(j, m) d_c += fabs(c_cap[j] - c_s[j]);
	return d_r + d_c;
}

double PrimalObj() {
	double total_cost = 0;
	FOR(i, n) FOR(j, m) {
		total_cost += x[i][j] * cost[i][j];
	}
	return total_cost;
}

void CalcSumAll() {
	FOR(i, n) r_s[i] = 0;
	FOR(j, m) c_s[j] = 0;
	FOR(i, n) FOR(j, m) {
		r_s[i] += x[i][j];
		c_s[j] += x[i][j];
	}
}
	
void DualToPrimalAll() {
	FOR(i, n) FOR(j, m) {
		double reduced_cost = y_left[i] + y_right[j] - cost[i][j];
		x[i][j] = exp(eta * reduced_cost); //this is x_{ij}
	}
	CalcSumAll();
}

void DualToPrimalRow(int i) {
	r_s[i] = 0;
	FOR(j, m) {
		c_s[j] -= x[i][j];
		double reduced_cost = y_left[i] + y_right[j] - cost[i][j];
		x[i][j] = exp(eta * reduced_cost); //this is x_{ij}
		r_s[i] += x[i][j];
		c_s[j] += x[i][j];
	}
}

void DualToPrimalCol(int j) {
	c_s[j] = 0;
	FOR(i, n) {
		r_s[i] -= x[i][j];
		double reduced_cost = y_left[i] + y_right[j] - cost[i][j];
		x[i][j] = exp(eta * reduced_cost); //this is x_{ij}
		r_s[i] += x[i][j];
		c_s[j] += x[i][j];
	}
}


double DualObj() {
	double ans = 0;
	FOR(i, n) ans += y_left[i] * r_cap[i];
	FOR(j, m) ans += y_right[j] * c_cap[j];
	return ans;
}


double FindMultiplier(int t, double goal) {
//find x s.t. \sum_{0<=i<t} exp(temp[i] + x) == goal
	double x = -temp[0];
	FR(i, 1, t) {
		SetMin<double>(x, -temp[i]);
	}
	double s = 0;
	FOR(i, t) {
		double v = temp[i] + x;
		if(v > -100 / eta) { //so no underflow happens...
			s += exp(v * eta);
		}
	}
	x -= log(s / goal) / eta;
	return x;
}

void WorkRow(int i) {
	FOR(j, m) temp[j] = y_right[j] - cost[i][j];
	y_left[i] = FindMultiplier(m, r_cap[i]);
}

void WorkColumn(int j) {
	FOR(i, n) temp[i] = y_left[i] - cost[i][j];
	y_right[j] = FindMultiplier(n, c_cap[j]);
}

void Show() {
	printf("%10.6lf...", -1);
	FOR(j, m) printf("%10.6lf ", y_right[j]);
	printf("\n");
	FOR(i, n) {
		printf("%10.6lf:::", y_left[i]);
		FOR(j, m) {
			printf("%10.6lf ", y_left[i] + y_right[j] - cost[i][j]);
		}
		printf("\n");
	}
}

int main() {
	scanf("%d%d", &n, &m);
	FOR(i, n) assert(scanf("%lf", &r_cap[i]) == 1);
	FOR(j, m) assert(scanf("%lf", &c_cap[j]) == 1);
	FOR(i, n) FOR(j, m) assert(scanf("%lf", &cost[i][j]) == 1);

double max_entry = 0;
	FOR(i, n) FOR(j, m) SetMax<double>(max_entry, cost[i][j]);
printf("max entry = %lf\n", max_entry);
eta /= max_entry;

	cerr << "DONE WITH IO" << endl;
	FOR(i, n) y_left[i] = 0;
	FOR(j, m) y_right[j] = 0;

	srand(time(0));
	int iter = 0;
	auto t1 = std::chrono::high_resolution_clock::now();
	while(1) {
//if(iter > 100) exit(0);
		if(iter % n == 0) {
			DualToPrimalAll();
		}

		iter++;
		double r_gain_max = -1;
		int i_rec;
		double c_gain_max = -1;
		int j_rec;
		FOR(i, n) {
			double r_gain = r_s[i] - r_cap[i] + r_cap[i] * log(r_cap[i] / r_s[i]);
			if(r_gain > r_gain_max) {
				r_gain_max = r_gain;
				i_rec = i;
			}
		}
		FOR(j, m) {
			double c_gain = c_s[j] - c_cap[j] + c_cap[j] * log(c_cap[j] / c_s[j]);
			if(c_gain > c_gain_max) {
				c_gain_max = c_gain;
				j_rec = j;
			}
		}
		if(r_gain_max > c_gain_max) {
			WorkRow(i_rec);
			DualToPrimalRow(i_rec);
		} else {
			WorkColumn(j_rec);
			DualToPrimalCol(j_rec);
		}

		double temp = TotalResiduals();
		
//printf("%d   %lf %d    %lf %d\n", iter, r_gain_max, i_rec, c_gain_max, j_rec);
		if(temp < epsilon) break;
		if(iter % 10000 == 0) {
			printf("Iteration %d:\n", iter);
			printf("TotalResiduals = %lf\n", temp);
			
			printf("Primal Objective Value = %g\n", PrimalObj());
			printf("Dual Objective Value = %g\n", DualObj());
		}
	}
	printf("%d iterations\n", iter);

//rounding stage, copied from round_transpoly.m from AWR17
	FOR(i, n) {
		double s = 0;
		FOR(j, m) s += x[i][j];
		double scale = min(1.0, r_cap[i] / s);
		FOR(j, m) x[i][j] *= scale;
	}
	FOR(j, m) {
		double s = 0;
		FOR(i, n) s += x[i][j];
		double scale = min(1.0, c_cap[j] / s);
		FOR(i, n) x[i][j] *= scale;
	}
	CalcSumAll();
	double temp = TotalResiduals();
	printf("Residuals before final low rank update = %lf\n", temp);
	FOR(i, n) FOR(j, m) {
		x[i][j] += (r_s[i] - r_cap[i]) * (c_s[j] - c_cap[j]) / temp * double(2.0);
	}

	CalcSumAll();
	printf("Residuals after rounding = %lf\n", TotalResiduals());

	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Time: " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << std::endl;
	printf("Primal Objective Value = %g\n", PrimalObj());
	printf("Dual Objective Value = %g\n", DualObj());
//FOR(i, n) printf("%lf\n", r_routed[i]);
//FOR(j, m) printf("%lf\n", c_routed[j]);
	return 0;
}


// Sinkhorn that implicitly tracks the dual variables
// (instead of keeping all exponentials around)
//
// eta is aggressiveness parameter, set it to smaller values for faster conv
//
// NOTE: this is unoptimized code strictly to measure convergence
// behavior vs iterations, do not try to time this
//
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
double y_left[MAXNM], y_right[MAXNM];

double eta = 5000; //100; //0.002
double epsilon = 1.0; 

double r_routed[MAXNM], c_routed[MAXNM];
int PRINT = 1;

double TotalResiduals() {
	FOR(i, n) r_routed[i] = r_cap[i];
	FOR(j, m) c_routed[j] = c_cap[j];

	FOR(i, n) FOR(j, m) {
		r_routed[i] -= x[i][j];
		c_routed[j] -= x[i][j];
	}
	
	double d_r = 0;
	FOR(i, n) d_r += fabs(r_routed[i]);
	double d_c = 0;
	FOR(j, m) d_c += fabs(c_routed[j]);
//printf("row diff = %lf col diff = %lf\n", d_r, d_c);
	return d_r + d_c;
}

double PrimalObj() {
	double total_cost = 0;
	FOR(i, n) FOR(j, m) {
		total_cost += x[i][j] * cost[i][j];
	}
	return total_cost;
}

void DualToPrimal() {
	double max_infeasibility = 0;
	FOR(i, n) FOR(j, m) {
		double reduced_cost = y_left[i] + y_right[j] - cost[i][j];
		SetMax<double>(max_infeasibility, reduced_cost);
		x[i][j] = exp(eta * reduced_cost); //this is x_{ij}
	}
//	if(PRINT == 1) printf("max infeasibility of dual solution = %lf\n", max_infeasibility * eta);
}

double DualObj() {
	double ans = 0;
	FOR(i, n) ans += y_left[i] * r_cap[i];
	FOR(j, m) ans += y_right[j] * c_cap[j];
	return ans;
}


double FindMultiplier(int t, double goal) {
//find x s.t.
//  \sum_{0<=i<t} exp(temp[i] + x) == goal
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
//printf("!@#$ %lf\n", log(s / goal) / eta);
	x -= log(s / goal) / eta;
//printf("%lf\n", x);

//printf("::%lf vs %lf  ===> %lf\n", s, goal, x);
//FR(i, 70, 75) printf("%d: %lf\n", i, temp[i] + x);
s = 0; FOR(i, t) s += exp(eta * (temp[i] + x));
assert(fabs(s - goal) < 0.01); 
/*
printf("%lf vs %lf\n", s, goal);
if(fabs(s - goal) > 0.01) {
	printf("%0.20lf\n", s - goal);
	FOR(i, t) temp[i] = (temp[i] + x) * eta;
	sort(temp, temp + t);
	reverse(temp, temp + t);
	FOR(i, 10) printf("%lf %lf\n", temp[i], exp(temp[i]));	
}
*/
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

//double s = 0; FOR(i, n) s += r_cap[i]; FOR(j, m) s -= c_cap[j]; assert(fabs(s) < 0.1);
//printf("%lf --> %lf\n", cost[0][0], cost[0][0] * eta); exit(0);
	cerr << "DONE WITH IO" << endl;
	FOR(i, n) y_left[i] = 0;
	FOR(j, m) y_right[j] = 0;

//puts("Init");
//Show();

//double sr = 0; FOR(i, n) sr += r_cap[i];
//double sc = 0; FOR(j, m) sc += c_cap[j];
//printf("%lf %lf\n", sr, sc);

//WorkColumn(90); exit(0);
	srand(time(0));
	int iter = 0;
	auto t1 = std::chrono::high_resolution_clock::now();
	while(1) {
//if(iter > 100) exit(0);
		iter++;
		FOR(i, n) {
			WorkRow(i);
//double s = 0;
//FOR(j, m) s += exp(eta * (y_left[i] + y_right[j] - cost[i][j]));
//printf("%d: %lf %lf\n", i, s, r_cap[i]); assert(fabs(s - r_cap[i]) < 0.1);
		}
		DualToPrimal();
//double s = 0;
//FOR(j, m) s += exp(eta * (y_left[0] + y_right[j] - cost[0][j]));
//printf("0: %lf %lf\n", s, r_cap[0]); assert(fabs(s - r_cap[0]) < 0.1);

//printf("%lf\n", TotalResiduals());
//FOR(i, n) {printf("%d: %lf/%lf\n", i, r_routed[i], r_cap[i]); assert(fabs(r_routed[i]) < 1e-9); }
		FOR(j, m) {
			WorkColumn(j);
		}
//printf("iteration %d:::\n", blah);
//Show();
		DualToPrimal();
		double temp = TotalResiduals();
//FOR(ij, max(n, m)) printf("%d %lf %lf\n", ij, r_routed[ij], c_routed[ij]);
		if(temp < epsilon) break;
		if(PRINT == 1 && iter % 10 == 0) {
			printf("Iteration %d:\n", iter);
			printf("TotalResiduals = %lf\n", temp);
			
			printf("Primal Objective Value = %lf\n", PrimalObj());
			printf("Dual Objective Value = %lf\n", DualObj());
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
	double temp = TotalResiduals();
//FOR(i, n) printf("%lf %lf\n", r_routed[i], c_routed[i]); //debug statement assuming n == m
	printf("Residuals before final low rank update = %lf\n", temp);
	FOR(i, n) FOR(j, m) {
//double temp1 = r_routed[i] * c_routed[j] / temp;
//if(fabs(temp1) > 1e-3) printf("%d %d %g\n", i, j, temp1);
		x[i][j] += r_routed[i] * c_routed[j] / temp * double(2);
	}
	printf("Residuals after rounding = %lf\n", TotalResiduals());
//FOR(i, n) printf("%lf %lf\n", r_routed[i], c_routed[i]); //debug statement assuming n == m
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Time: " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << std::endl;
	printf("Primal Objective Value = %lf\n", PrimalObj());
	printf("Dual Objective Value = %lf\n", DualObj());
//FOR(i, n) printf("%lf\n", r_routed[i]);
//FOR(j, m) printf("%lf\n", c_routed[j]);
	return 0;
}


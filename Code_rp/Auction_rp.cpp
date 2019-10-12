//
// Auction code for uncapacitied (matching) instnaces
//
//Runtime of this method depends heavily on how much
// the initial edge costs get divided by
//    It's set using the RESCALE parameter below.
//

#include <cstdio>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cassert>
using namespace std;

// Key rescaling parameter for tuning time / error tradeoffs
const int RESCALE = 1000; 

const int INF = 1<<30;
const int MAXNM = 5010;

int n, m;
int temp;
int a[MAXNM][MAXNM], p[MAXNM];
int to[MAXNM], from[MAXNM];

#define FR(i, a, b) for(int i=(a); i<(b); i++) 
#define FOR(i, n) FR(i, 0, n) 

template <typename T> inline void SetMin(T &a, T b) {if(b < a) a = b;}
template <typename T> inline void SetMax(T &a, T b) {if(b > a) a = b;}

double r[MAXNM], c[MAXNM];
double cost[MAXNM][MAXNM];

void Eval() {	
	int unmatched = 0;
	double total_cost = 0;
	FOR(i, n) {
		if(to[i] == -1) unmatched++;
		total_cost += cost[i][to[i]];
	}
	printf("%d UNMATCHED\n", unmatched);
	printf("total = %lf\n", total_cost);
	fflush(stdout);
}

void SanityCheck() {
	puts("---checking if is (partial) matching");
	FOR(i1, n) if(to[i1] != -1) assert(from[to[i1]] == i1);
	FOR(j1, m) if(from[j1] != -1) { if(to[from[j1]] != j1) printf("%d %d %d\n", from[j1], to[from[j1]], j1); assert(to[from[j1]] == j1);}
}

int main() {
	scanf("%d%d", &n, &m);
	FOR(i, n) scanf("%lf", &r[i]);
	FOR(j, m) scanf("%lf", &c[j]);
	FOR(i,n) FOR(j, m) scanf("%lf", &cost[i][j]);
	puts("DONE WITH IO");
	fflush(stdout);

	clock_t time1 = clock();
	FOR(i, n) FOR(j, m) a[i][j] = -(cost[i][j] / RESCALE);
	FOR(j, m) p[j] = 0;
	

	int matched = 0;
	FOR(i, n) to[i] = -1;
	FOR(j, m) from[j] = -1;
	
	int iter = 0;

	while(matched < n) {
		++iter;	
		FOR(i, n) {
			int j1, j2;
			int rec1 = -INF, rec2 = -INF;
			FOR(j, m) {
				int gain = a[i][j] - p[j];
				if(gain > rec1) {
					rec2 = rec1;
					j2 = j1;
					rec1 = gain;
					j1 = j;
				} else if(gain > rec2) {
					rec2 = gain;
					j2 = j;
				}
			}
//printf("%d (%d): %d (%d)   %d  (%d)\n", i, to[i], j1, from[j1], j2, from[j2]);
			if(to[i] != -1 && a[i][to[i]] - p[to[i]] < rec1 - 1) {
				from[to[i]] = -1;
				to[i] = -1;
			}
			if(to[i] == -1) {
				if(from[j1] == -1) {
					to[i] = j1;
					from[j1] = i;
					matched++;
				} else {
					to[from[j1]] = -1;
					to[i] = j1;
					from[j1] = i;
					p[j1] += rec1 - rec2 + 1;
				}
			}
//printf("%d (%d): %d (%d)   %d  (%d)\n", i, to[i], j1, from[j1], j2, from[j2]);
		}
		if(iter % 100 == 0) {
			printf("--iteration %d:::", iter);
			printf("%d matched, time = %lfs\n", matched, ((clock() - time1)/(double)CLOCKS_PER_SEC)); //break;
		}
//		Eval();
	}
	SanityCheck();
//FOR(i, n) printf("%d %d\n", to[i], from[i]);
	printf("%d iterations\n", iter);
	printf("clock() = %lf seconds\n", (clock()-time1)/(double)CLOCKS_PER_SEC );
	Eval(); 
	return 0;
}

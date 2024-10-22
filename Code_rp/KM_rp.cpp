//
//KM optimized for dense + small cost instances
//
//This version runs on DFS + implicit vertex potential tracking
//+ adjacency list to track of admissible edges
//
//This version uses CSR format for the sparse edges
//  in intermediate steps
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

const int INF = 1<<30;
const int MAXNM = 5010;

// Key rescaling parameter for tuning time / error tradeoffs
const int RESCALE = 1000; 

int n, m;
int temp;
int a[MAXNM][MAXNM], x[MAXNM], y[MAXNM];
int original[MAXNM][MAXNM];

int to[MAXNM], from[MAXNM];
int d[MAXNM];
char vis[MAXNM];

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

short e[MAXNM * MAXNM];
int estart[MAXNM], et;

int DFS(int i) {
	vis[i] = 1;
	FR(i1, estart[i], estart[i + 1]) {
		if(from[e[i1]] == -1) {
			to[i] = e[i1];
			from[e[i1]] = i;
			return 1;
		}
	}
	FR(i1, estart[i], estart[i + 1]) {
		if(!vis[from[e[i1]]] && DFS(from[e[i1]])) {
			to[i] = e[i1];
			from[e[i1]] = i;
			return 1;
		}
	}
	return 0;
}

int main() {
	scanf("%d%d", &n, &m);
	FOR(i, n) scanf("%lf", &r[i]);
	FOR(j, m) scanf("%lf", &c[j]);
	FOR(i,n) FOR(j, m) scanf("%lf", &cost[i][j]);
	puts("DONE WITH IO");

	clock_t time1 = clock();
	FOR(i, n) x[i] = 0;
	FOR(j, m) y[j] = 0;
	FOR(i, n) FOR(j, m) original[i][j] = cost[i][j] / RESCALE;
	FOR(i, n) FOR(j, m) a[i][j] = min(1<<20, original[i][j]); //min in case we want to use doubles or something as such

	FOR(i, n) to[i] = -1;
	FOR(j, m) from[j] = -1;
	int matched = 0;
	
	int iter = 0;

	while(matched < n) {
		printf("--iteration %d:::", ++iter);
		et = 0;
		FOR(i, n) {
			estart[i] = et;
			FOR(j, m) if(x[i] + y[j] == a[i][j]) {
				e[et++] = j;
			}
		}
		estart[n] = et;
printf("%d admissible arcs, ", et);
		int progress = 1;
		while(progress) {
			progress = 0;
			memset(vis, 0 ,sizeof(vis));
			FOR(i, n) if(to[i] == -1 && DFS(i)) {
				progress = 1;
				matched++;
			}
		}
		FOR(i, n) if(vis[i]) {
			x[i]++;
			if(to[i] != -1) y[to[i]]--;
		}
		printf("%d matched, time = %lfs\n", matched, ((clock() - time1)/(double)CLOCKS_PER_SEC)); //break;
	}
	printf("clock() = %lf seconds\n", (clock()-time1)/(double)CLOCKS_PER_SEC );
	Eval(); 
	return 0;
}

//
//MCMF based on KM,
// optimized for dense + small cost instances
//
// this version works general capacitated instances
//
// This version just does dinic on the tight (cost 0) dual edges,
//  and bounces flow forward and backwards using tables
//
//Runtime of this method depends heavily on how much
// the initial edge costs get divided by
//    It's set using the MULT_F/MULT_C parameters below.
//

#include <cstdio>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cassert>
using namespace std;

#define FR(i, a, b) for(int i=(a); i<(b); i++) 
#define FOR(i, n) FR(i, 0, n) 

template <typename T> inline void SetMin(T &a, T b) {if(b < a) a = b;}
template <typename T> inline void SetMax(T &a, T b) {if(b > a) a = b;}

const int INF = 1000000005;
const int MAXNM = 5010;

int n, m;
int temp;
int d_l[MAXNM], d_r[MAXNM];
int a[MAXNM][MAXNM], x[MAXNM], y[MAXNM];
int original[MAXNM][MAXNM];
int f[MAXNM][MAXNM];
int eid[MAXNM][MAXNM];

// Key rescaling parameters for tuning time / error tradeoffs
const double MULT_F = 100000;
const double MULT_C = 0.001;

double r[MAXNM], c[MAXNM];
double cost[MAXNM][MAXNM];

void Eval() {
	int m_used = 0;
	double total_cost = 0;
	vector<double> temp_l(r, r + n);
	vector<double> temp_r(c, c + m);
	FOR(i, n) {
		FOR(j, m) {
			if(f[i][j] > 0) m_used++;
			double f1 = double(f[i][j]) / MULT_F;
		 	temp_l[i] -= f1;			
		 	temp_r[j] -= f1;			
			total_cost += cost[i][j] * f1;;
		}
	}
	printf("%d edges with non-zero flows\n", m_used);

	double residue_l1 = 0;
	FOR(i, n) residue_l1 += fabs(temp_l[i]);
	FOR(j, m) residue_l1 += fabs(temp_r[j]);
	printf("difference in demand = %lf\n", residue_l1);
	printf("total cost = %0.6g\n", total_cost);
	fflush(stdout);
}

const int MAXV = MAXNM * 2 + 2;
const int MAXE = 1234567;
//WARNING: these are only meant to handle int capacities,
// need different values if we want 64-bit caps

int V, source, sink;
int eind;
int eadj[MAXE], eprev[MAXE], elast[MAXV], start[MAXV];
int ecap[MAXE];
int front, back, q[MAXV], dist[MAXV];

inline void AddEdge (int u, int v, int cap_uv, int cap_vu) {
//printf("adding edge %d-->%d, (%d, %d)\n", u, v, cap_uv, cap_vu);
	eadj[eind] = v;
	ecap[eind] = cap_uv;
	eprev[eind] = elast[u];
	elast[u] = eind++;
	
	eadj[eind] = u;
	ecap[eind] = cap_vu;
	eprev[eind] = elast[v];
	elast[v] = eind++;
}

bool BFS(){
	memset (dist, 63, V * sizeof (int));
	front = back = 0;
	q[back++] = source; dist [source] = 0;

	while (front < back) {
		int top = q[front++];
		for (int e = elast[top]; e != -1; e = eprev[e])
			if (ecap [e] > 0 &&
				dist[top] + 1 < dist[eadj[e]]) {
				dist[eadj[e]] = dist[top] + 1;
				q[back++] = eadj[e];
		}
	}
	return dist[sink] < INF;
}

int DFS(int u, int pcap) {
	if (u == sink) return pcap;
	int total = 0;

	for (int &e = start[u]; e != -1; e = eprev[e]) {
		if (ecap[e] > 0 && dist[u] + 1 == dist[eadj[e]]) {
			int p = DFS(eadj[e], min (pcap, ecap[e]));
			ecap[e] -= p;
			ecap[e ^ 1] += p;
   			pcap -= p;
			total += p;	
			if (pcap == 0) break;
		}
	}
	return total;
}

int main() {
	int cn = 0;
	while(scanf("%d%d", &n, &m) == 2) {
++cn;
printf("=====case %d====\n", cn);
		FOR(i, n) scanf("%lf", &r[i]);
		FOR(j, m) scanf("%lf", &c[j]);
		FOR(i,n) FOR(j, m) scanf("%lf", &cost[i][j]);
		puts("DONE WITH IO");

		clock_t time1 = clock();
		FOR(i, n) x[i] = 0;
		FOR(j, m) y[j] = 0;
		FOR(i, n) FOR(j, m) original[i][j] = int(cost[i][j] * MULT_C);
int max_cost = 0;
FOR(i, n) FOR(j, m) SetMax<int>(max_cost, original[i][j]);
printf("MAX COST = %d\n", max_cost);

		FOR(i, n) FOR(j, m) a[i][j] = min(1<<25, original[i][j]); //min in case we want to use doubles or something as such

		int s_l = 0;
		FOR(i, n) {
			d_l[i] = int(r[i] * MULT_F); 
			s_l += d_l[i];
		}
		int s_r = 0;
		FOR(j, m) {
			d_r[j] = int(c[j] * MULT_F);
			s_r += d_r[j];
		}

		int iter = 0;
		int goal = min(s_l, s_r);
		printf("TOTAL DEMAND = %d\n", goal);
		while(goal > 0 && iter < 200) {
			printf("--iteration %d::: %d units to go:::", ++iter, goal);
			V = n + m + 2;
			memset(elast, -1, V * sizeof(int));
			source = 0;
			sink = 1;
			eind = 0;

			FOR(i, n) if(d_l[i] > 0) {
				AddEdge(0, 2 + i, d_l[i], 0);
			}
			FOR(j, m) if(d_r[j] > 0) {
				AddEdge(2 + n + j, 1, d_r[j], 0);
			}

			FOR(i, n) FOR(j, m) if(x[i] + y[j] == a[i][j]) {
				AddEdge(2 + i, 2 + n + j, INF, f[i][j]);
				f[i][j] = eind - 1;
				//point capacity to the reverse capacity value
			} else {
				f[i][j] =-1;
			}
printf("%d admissible arcs ... ", eind); fflush(stdout);

			
int n_passes = 0;
			while (BFS()) {
n_passes++;
				memcpy (start, elast, V * sizeof (int));
				DFS(source, INF);
			}
printf("%d dinic passes", n_passes); fflush(stdout);

			FOR(i, n) FOR(j, m) {
				if(f[i][j] == -1) {
					f[i][j] = 0;
				} else {
					f[i][j] = ecap[f[i][j]];
					//bring the `reverse capacity' back
				}
			}
			
			for (int e = elast[source]; e != -1; e = eprev [e]) {
				goal -= ecap[e ^ 1];
				d_l[eadj[e] - 2] -= ecap[e ^ 1];
			}
			for (int e = elast[sink]; e != -1; e = eprev[e]) {
				d_r[eadj[e] - 2 - n] -= ecap[e];
			}

printf("time = %lfs\n", ((clock() - time1)/(double)CLOCKS_PER_SEC)); //break;

//KM-style dual adjustment
			FOR(i, n) if(dist[2 + i] < INF) x[i]++;
			FOR(j, m) if(dist[2 + n + j] < INF) y[j]--;
		}
//FOR(i, n) printf("%dL:::%d\n", i, d_l[i]);
//FOR(j, m) printf("%dR:::%d\n", j, d_r[j]);
		printf("TOTAL SOLVE TIME = %lf seconds\n", (clock()-time1)/(double)CLOCKS_PER_SEC );
		Eval(); 
	}
	return 0;
}

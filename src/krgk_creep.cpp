// krgk_creep: T -> 0+ creep dynamics of a directed elastic line in a random potential.
//
// Algorithm of A. B. Kolton, A. Rosso, T. Giamarchi and W. Krauth,
//   PRL 97, 057001 (2006)  and  PRB 79, 184207 (2009).
// Cleaned-up port of the creep driver ("-o") of the original code (original/vmc), with the
// same energy, moves, disorder and random numbers, so that both give the same events.
//
// Model: u_i integer, i = 0..L-1 periodic, u_i in a random potential periodic in u with period M,
//   E = sum_i [ V(i, u_i mod M) - F u_i + (c/2) (u_{i+1} - u_i)^2 ],  |u_{i+1} - u_i| <= mcut.
// Elementary moves: rigid displacement by +-1 of a block of consecutive sites (any length,
// including the whole line), allowed if the two border slopes stay within mcut.
// Creep step from a metastable state alpha: the minimal-barrier path, among the paths of
// forward elementary moves, to a metastable state gamma with E(gamma) < E(alpha). It is found
// by flooding: configurations are expanded in increasing order of energy, each one relaxed
// forward to its metastable state before its children are generated; the barrier is the
// highest energy (minus E(alpha)) of the configurations taken from the queue.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

static const double BIG = 1e300;
static const double PRECISION_FC = 0.005;  // step in F when locating F_c (as in the original)

struct Params {
	int L = 16, M = 64, mcut = 1;
	double disorder = 1.0, c = 2.0, f = 0.08;
	bool f_absolute = false;
	unsigned long seed = 12345;
	long events = 100;
	int limit = 5;            // search envelope: u_i <= alpha_i + limit*M
	int conf_every = 0;       // write the metastable configuration every n events (0: never)
	std::string out = ".";
};

// Mersenne Twister and Gaussian deviates exactly as in the original (mt.cpp):
// genrand_real3 in (0,1) and the polar Box-Muller of Numerical Recipes.
class Random {
	std::mt19937 mt;
	bool have = false;
	double saved = 0.0;
	double real3() { return ((double)mt() + 0.5) * (1.0 / 4294967296.0); }
public:
	explicit Random(unsigned long s) : mt((uint32_t)s) {}
	double gauss()
	{
		if (have) { have = false; return saved; }
		double v1, v2, r2;
		do {
			v1 = 2.0 * real3() - 1.0;
			v2 = 2.0 * real3() - 1.0;
			r2 = v1 * v1 + v2 * v2;
		} while (r2 >= 1.0 || r2 == 0.0);
		double fac = std::sqrt(-2.0 * std::log(r2) / r2);
		saved = v1 * fac; have = true;
		return v2 * fac;
	}
};

class Line {
public:
	int L, M, mcut;
	double c, F = 0.0;
	std::vector<double> V;  // V[i*M + j]
	std::vector<int> u;

	Line(const Params &p) : L(p.L), M(p.M), mcut(p.mcut), c(p.c), V((size_t)p.L * p.M), u(p.L, 0)
	{
		Random rng(p.seed);
		for (auto &v : V) v = rng.gauss();
		double mean = 0.0;
		for (double v : V) mean += v;
		mean /= V.size();
		for (auto &v : V) v = (v - mean) * p.disorder;
	}

	double pot(int i, int j) const { int m = j % M; if (m < 0) m += M; return V[(size_t)i * M + m]; }

	double energy(const std::vector<int> &x) const
	{
		double E = 0.0;
		for (int i = 0; i < L; i++) {
			int d = x[i + 1 == L ? 0 : i + 1] - x[i];
			E = E + pot(i, x[i]) - F * x[i] + c * 0.5 * d * d;
		}
		return E;
	}

	// Energy change of the block move (i, il): sites i+1 .. i+il-1 moved by dir, sites i and
	// i+il fixed (il = 2..L); il = L+1 is the rigid move of the whole line. Returns false if a
	// border slope would exceed mcut, or (with lim) if a moved site would pass lim.
	bool block(const std::vector<int> &x, int dir, int i, int il, double &d, const std::vector<int> *lim) const
	{
		if (il <= L) {
			int istop = (i + il) % L, a = x[(i + 1) % L], b = x[(istop - 1 + L) % L];
			if (std::abs(a + dir - x[i]) > mcut || std::abs(b + dir - x[istop]) > mcut) return false;
			d = c * (1.0 + dir * (a - x[i] + b - x[istop]));
		} else d = 0.0;
		for (int j = 1; j < il; j++) {
			int p = (i + j) % L;
			if (lim && x[p] + dir > (*lim)[p]) return false;
			d = d - pot(p, x[p]) + pot(p, x[p] + dir) - dir * F;
		}
		return true;
	}

	// One elementary move in direction dir: applies the first block move (shortest first, then
	// leftmost) that lowers the energy. Returns the smallest energy change per moved site among
	// the moves tried; it is negative if and only if a move was made.
	double move(int dir)
	{
		double dmin = BIG, d;
		for (int il = 2; il <= L + 1; il++) {
			int nstart = (il <= L) ? L : 1;
			for (int i = 0; i < nstart; i++) {
				if (!block(u, dir, i, il, d, nullptr)) continue;
				dmin = std::min(dmin, d / (il - 1));
				if (d < 0) {
					for (int j = 1; j < il; j++) u[(i + j) % L] += dir;
					return dmin;
				}
			}
		}
		return dmin;
	}

	void relax_forward() { while (move(+1) < 0.0) {} }

	// Critical force with the procedure of the original (find_critical): starting from u = 0 and
	// the initial F, relax backward; then advance forward by 2M, raising F by (dmin + PRECISION_FC)
	// whenever no forward move lowers the energy. F_c = F - PRECISION_FC at the end.
	double critical_force(double F0)
	{
		F = F0;
		std::fill(u.begin(), u.end(), 0);
		while (move(-1) < 0.0) {}
		F -= 5.0 * PRECISION_FC;
		int count = 0;
		while (count == 0) {
			std::vector<int> u0 = u;
			auto advanced = [&]() { int m = u[0] - u0[0]; for (int i = 1; i < L; i++) m = std::min(m, u[i] - u0[i]); return m; };
			while (advanced() < 2 * M) {
				double d = move(+1);
				if (d >= 0.0) { F += d + PRECISION_FC; count++; }
			}
			if (count == 0) { F -= 5.0 * PRECISION_FC; for (auto &x : u) x -= 2 * M; }
		}
		return F - PRECISION_FC;
	}

	// Compact key of a configuration: u_0 and the slopes, each in ceil(log2(2 mcut + 1)) bits.
	std::string key(const std::vector<int> &x) const
	{
		int bits = 1;
		while ((1 << bits) < 2 * mcut + 1) bits++;
		std::string k(sizeof(int) + ((size_t)(L - 1) * bits + 7) / 8, '\0');
		std::memcpy(&k[0], &x[0], sizeof(int));
		size_t pos = 0;
		for (int i = 1; i < L; i++) {
			unsigned s = (unsigned)(x[i] - x[i - 1] + mcut);
			for (int b = 0; b < bits; b++, pos++)
				if (s >> b & 1u) k[sizeof(int) + pos / 8] |= (char)(1 << (pos % 8));
		}
		return k;
	}

	std::vector<int> unkey(const std::string &k) const
	{
		int bits = 1;
		while ((1 << bits) < 2 * mcut + 1) bits++;
		std::vector<int> x(L);
		std::memcpy(&x[0], &k[0], sizeof(int));
		size_t pos = 0;
		for (int i = 1; i < L; i++) {
			unsigned s = 0;
			for (int b = 0; b < bits; b++, pos++)
				if (k[sizeof(int) + pos / 8] >> (pos % 8) & 1) s |= 1u << b;
			x[i] = x[i - 1] + (int)s - mcut;
		}
		return x;
	}

	struct Node {
		double E;
		std::string k;
		bool operator<(const Node &o) const { return E < o.E || (E == o.E && k < o.k); }
	};

	// Creep step from the current (metastable) configuration alpha. On return u is gamma, beta is
	// the last configuration taken from the queue (before its relaxation into gamma), and the
	// barrier max(E) - E(alpha) is returned; BIG if the envelope stops the search.
	double creep_step(std::vector<int> &beta, int limit, long &expanded)
	{
		const std::vector<int> alpha = u;
		std::vector<int> lim(L);
		for (int i = 0; i < L; i++) lim[i] = alpha[i] + limit * M;
		const double E0 = energy(u);
		double UB = 0.0, d;
		std::set<Node> queue;
		std::unordered_set<std::string> closed;
		expanded = 0;
		do {
			// children of the current metastable mother u
			bool inside = true;
			for (int i = 0; i < L; i++) inside = inside && u[i] <= lim[i];
			if (inside) {
				double Em = energy(u);
				for (int il = 2; il <= L + 1; il++) {
					int nstart = (il <= L) ? L : 1;
					for (int i = 0; i < nstart; i++) {
						if (!block(u, +1, i, il, d, &lim)) continue;
						std::vector<int> x = u;
						for (int j = 1; j < il; j++) x[(i + j) % L] += 1;
						std::string k = key(x);
						if (!closed.count(k)) queue.insert({Em + d, std::move(k)});
					}
				}
			}
			// lowest-energy configuration not expanded yet
			Node n;
			do {
				if (queue.empty()) { u = alpha; return BIG; }
				n = *queue.begin();
				queue.erase(queue.begin());
			} while (closed.count(n.k));
			closed.insert(n.k);
			expanded++;
			u = unkey(n.k);
			UB = std::max(UB, n.E - E0);
			beta = u;
			relax_forward();
		} while (energy(u) > E0);
		return UB;
	}
};

static void usage(const char *prog)
{
	std::cerr << "usage: " << prog << " [options]\n"
	          << "  --L N          line length (sites, periodic)              [16]\n"
	          << "  --M N          period of the disorder along u            [64]\n"
	          << "  --mcut N       maximal slope |u_{i+1}-u_i|               [1]\n"
	          << "  --disorder X   disorder strength (Gaussian, std X)       [1.0]\n"
	          << "  --c X          elastic constant                          [2.0]\n"
	          << "  --f X          drive F/F_c (or F with --absolute)        [0.08]\n"
	          << "  --absolute     --f is the absolute force F\n"
	          << "  --seed N       random seed                               [12345]\n"
	          << "  --events N     number of creep steps                     [100]\n"
	          << "  --limit N      search envelope u <= alpha + N*M          [5]\n"
	          << "  --conf-every N write the metastable state every N steps  [0: never]\n"
	          << "  --out DIR      output directory                          [.]\n";
}

int main(int argc, char **argv)
{
	Params p;
	for (int a = 1; a < argc; a++) {
		std::string o = argv[a];
		auto val = [&]() { if (a + 1 >= argc) { usage(argv[0]); std::exit(1); } return std::string(argv[++a]); };
		if (o == "--L") p.L = std::stoi(val());
		else if (o == "--M") p.M = std::stoi(val());
		else if (o == "--mcut") p.mcut = std::stoi(val());
		else if (o == "--disorder") p.disorder = std::stod(val());
		else if (o == "--c") p.c = std::stod(val());
		else if (o == "--f") p.f = std::stod(val());
		else if (o == "--absolute") p.f_absolute = true;
		else if (o == "--seed") p.seed = std::stoul(val());
		else if (o == "--events") p.events = std::stol(val());
		else if (o == "--limit") p.limit = std::stoi(val());
		else if (o == "--conf-every") p.conf_every = std::stoi(val());
		else if (o == "--out") p.out = val();
		else { usage(argv[0]); return 1; }
	}
	if (p.L < 3 || p.M < 2 || p.mcut < 1) { std::cerr << "need L >= 3, M >= 2, mcut >= 1\n"; return 1; }

	Line line(p);
	// F_c is located starting from F = --f, as the original does (it reads F/F_c as the initial F).
	double Fc = line.critical_force(p.f);
	line.F = p.f_absolute ? p.f : p.f * Fc;
	std::fill(line.u.begin(), line.u.end(), 0);
	line.relax_forward();

	std::ofstream log(p.out + "/params.txt");
	log << "L " << p.L << "\nM " << p.M << "\nmcut " << p.mcut << "\ndisorder " << p.disorder << "\nc " << p.c
	    << "\nf " << p.f << (p.f_absolute ? " (absolute)" : " (F/Fc)") << "\nseed " << p.seed << "\nlimit " << p.limit
	    << "\nFc " << Fc << "\nF " << line.F << "\n";
	std::ofstream ev(p.out + "/events.dat");
	ev << "# step barrier L_opt L_rel area mean_u E expanded\n";
	ev.precision(10);
	std::ofstream conf;
	if (p.conf_every > 0) conf.open(p.out + "/configs.bin", std::ios::binary);

	std::vector<int> alpha, beta(p.L);
	for (long s = 1; s <= p.events; s++) {
		alpha = line.u;
		long expanded;
		double UB = line.creep_step(beta, p.limit, expanded);
		if (UB >= BIG) { std::cerr << "step " << s << ": no lower state inside the envelope (raise --limit)\n"; return 2; }
		int Lopt = 0, Lrel = 0;
		long area = 0, sum = 0;
		for (int i = 0; i < p.L; i++) {
			Lopt += alpha[i] != beta[i];
			Lrel += alpha[i] != line.u[i];
			area += line.u[i] - alpha[i];
			sum += line.u[i];
		}
		ev << s << " " << UB << " " << Lopt << " " << Lrel << " " << area << " " << (double)sum / p.L << " "
		   << line.energy(line.u) << " " << expanded << "\n";
		if (p.conf_every > 0 && s % p.conf_every == 0) conf.write((const char *)line.u.data(), sizeof(int) * p.L);
	}
	return 0;
}

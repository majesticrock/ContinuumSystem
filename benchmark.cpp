#include <cmath>

static constexpr double __A = 4.;
static constexpr double __SINH_A = 27.289917197127752448908271590793818580289412485530296565528497015;

static constexpr double fermi_wavevector = 4.5;
static constexpr double K_MAX = 5;
static constexpr double K_MIN = 4;

static constexpr int N = 1000;
static constexpr double STEP = (K_MAX - K_MIN) / N;
static constexpr double scaling = 2. * __A / N;

static inline double index_sinh(int k_idx) {
	return fermi_wavevector + 0.5 * (K_MAX - K_MIN) * std::sinh(scaling * k_idx - __A) / __SINH_A;
}
static inline int momentum_sinh(double k) {
	k = (k - fermi_wavevector) / (K_MAX - K_MIN);
	return static_cast<int>(std::lround(
			(std::asinh(__SINH_A * k) + __A) / scaling
		));
}

static inline double index_lin(int k_idx) {
	return K_MIN + STEP * k_idx;
}
static inline int momentum_lin(double k) {
	return static_cast<int>(std::lround((k - K_MIN) / STEP));
}



static void __sinh(benchmark::State& state) {
	for (auto _ : state) {
		for (int i = 0; i < N; ++i) {
			double res = index_sinh(i);
			int pos = momentum_sinh(res);
		}
	}
}
BENCHMARK(__sinh);

static void __lin(benchmark::State& state) {
	for (auto _ : state) {
		for (int i = 0; i < N; ++i) {
			double res = index_lin(i);
			int pos = momentum_lin(res);
		}
	}
}
BENCHMARK(__lin);
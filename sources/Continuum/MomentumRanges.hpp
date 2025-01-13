#pragma once
#include "GlobalDefinitions.hpp"
#include <concepts>
#include <boost/math/quadrature/gauss.hpp>
#include <mrock/utility/Numerics/Integration/GeneralizedPrincipalValue.hpp>

namespace Continuum {
	struct MomentumRanges {
		static constexpr int n_gauss = 120;
#ifdef PHONON_SC_CHANNEL_ONLY
		static constexpr int n_dangerous_points = 3;
#else
		static constexpr int n_dangerous_points = 5;
#endif
		c_float K_MAX{};
		c_float K_MIN{};

		c_float INNER_K_MAX{};
		c_float INNER_K_MIN{};

		c_float LOWER_STEP{};
		c_float INNER_STEP{};
		c_float UPPER_STEP{};

		c_float const * K_F;
		std::array<c_float, 2> singularities_in_renormalized_dispersion;

		MomentumRanges(c_float const * k_F, const c_float omega_debye, c_float inner_offset);

		c_float index_to_momentum(int k_idx) const;
		int momentum_to_index(c_float k) const;
		int momentum_to_floor_index(c_float k) const;

		// Members that allow vector-like handling
		inline c_float operator[](int k_idx) const {
			return index_to_momentum(k_idx);
		}
		inline int size() const noexcept {
			return DISCRETIZATION;
		}

		std::vector<c_float> get_k_points() const;

		template<class Function>
		auto integrate(const Function& func, c_float begin, c_float end) const {
			decltype(func(begin)) value{ };
			if (is_zero(begin - end)) return value;
			if (begin < K_MIN) begin = K_MIN;
			if (end > K_MAX) end = K_MAX;

			int n{};
			while(n < n_dangerous_points && begin > dangerous_point(n)) {
				++n;
			}
			while(++n < n_dangerous_points && end > dangerous_point(n)) {
				value += __integrate(func, begin, dangerous_point(n - 1));
				begin = dangerous_point(n - 1);
			}
			value += __integrate(func, begin, end);
			return value;
		}

		template<class Function>
		auto cpv_integrate(const Function& func, c_float begin, c_float end, const c_float& singularity) const {
			auto distinction = [this, &singularity](c_float comp) {
				return ((std::abs(singularity - comp) < SINGULARITY_OFFSET) ? -1 : comp);
			};
			std::array<c_float, n_dangerous_points + 1> special_points;
			for (int i = 0; i < n_dangerous_points; ++i) {
				special_points[i] = distinction(dangerous_point(i));
			}
			special_points.back() = singularity;
			std::ranges::sort(special_points);
			return __cpv_integrate(func, begin, end, special_points);
		}

		template<class Function>
		inline auto integrate(const Function& func) const {
			return __integrate(func, K_MIN, INNER_K_MIN)
				+ __integrate(func, INNER_K_MIN, (*K_F))
				+ __integrate(func, (*K_F), INNER_K_MAX)
				+ __integrate(func, INNER_K_MAX, K_MAX);
		}

		template<class Function, class c_float>
		inline auto cpv_integrate(const Function& func, const c_float& singularity) const {
			auto distinction = [this, &singularity](c_float comp) {
				return ((std::abs(singularity - comp) < SINGULARITY_OFFSET) ? -1 : comp);
			};
			std::array<c_float, n_dangerous_points + 1> special_points;
			for (int i = 0; i < n_dangerous_points; ++i) {
				special_points[i] = distinction(dangerous_point(i));
			}
			special_points.back() = singularity;
			std::ranges::sort(special_points);
			return __cpv_integrate(func, K_MIN, K_MAX, special_points);
		}

	private:
		constexpr bool in_interval(c_float x, c_float lower, c_float upper) const {
			return ((x > lower) && (x < upper));
		}

		inline c_float dangerous_point(int i) const {
			assert(i < n_dangerous_points);
#ifdef PHONON_SC_CHANNEL_ONLY
			switch(i) {
				case 0:
					return INNER_K_MIN; break;
				case 1:
					return (*K_F); break;
				case 2:
					return INNER_K_MAX; break;
			}
#else
			switch(i) {
				case 0:
					return INNER_K_MIN; break;
				case 1:
					return singularities_in_renormalized_dispersion[0]; break;
				case 2:
					return (*K_F); break;
				case 3:
					return singularities_in_renormalized_dispersion[1]; break;
				case 4:
					return INNER_K_MAX; break;
			}
#endif
			throw std::runtime_error("Reached the end of MomentumRanges::dangerous_point");
		}

		template<class Function>
		inline auto __integrate(const Function& func, c_float begin, c_float end) const {
			if (is_zero(end - begin)) {
				return decltype(func(begin)){};
			}
			return boost::math::quadrature::gauss<c_float, n_gauss>::integrate(func, begin, end);
		}
		template<class Function, class Real_or_RandomAccessContainer>
		inline auto __cpv_integrate(const Function& func, c_float begin, c_float end, const Real_or_RandomAccessContainer& singularity) const {
			if (is_zero(end - begin)) {
				return decltype(func(begin)){};
			}
			return mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<c_float, n_gauss>::generalized_principal_value(func, begin, end, singularity);
		}
	};

	class MomentumIterator {
		MomentumRanges const* const _parent;
	public:
		c_float k{};
		int idx{};
		static inline int max_idx() { return DISCRETIZATION; }

		MomentumIterator(MomentumRanges const* const parent, int init = 0)
			: _parent(parent), k(_parent->index_to_momentum(init)), idx(init) {}

		inline c_float parent_step() const {
			if (k < _parent->INNER_K_MIN) return _parent->LOWER_STEP;
			if (k <= _parent->INNER_K_MAX) return _parent->INNER_STEP;
			return _parent->UPPER_STEP;
		}
		inline c_float max_k() const { return _parent->K_MAX; }
		inline c_float min_k() const { return _parent->K_MIN; }

		inline MomentumIterator& operator++() {
			++idx;
			k = _parent->index_to_momentum(idx);
			return *this;
		}
		inline MomentumIterator operator++(int) {
			MomentumIterator tmp = *this;
			++(*this);
			return tmp;
		}

		inline MomentumIterator& operator--() {
			--idx;
			k = _parent->index_to_momentum(idx);
			return *this;
		}
		inline MomentumIterator operator--(int) {
			MomentumIterator tmp = *this;
			--(*this);
			return tmp;
		}

		inline auto operator<=>(MomentumIterator const& other) const = default;
	};

	class InnerIterator {
		MomentumRanges const* const _parent;
		inline c_float index_to_momentum(int i) const {
			return (_parent->INNER_K_MIN + i * _parent->INNER_STEP);
		}
	public:
		c_float k{};
		int idx{};
		static inline int max_idx() { return _INNER_DISC + 1; }

		InnerIterator(MomentumRanges const* const parent, int init = 0)
			: _parent(parent), k(this->index_to_momentum(init)), idx(init) {}

		inline c_float parent_step() const {
			return _parent->INNER_STEP;
		}
		inline c_float max_k() const { return _parent->INNER_K_MAX; }
		inline c_float min_k() const { return _parent->INNER_K_MIN; }

		inline c_float operator[](int i) const { return this->index_to_momentum(i); }

		inline InnerIterator& operator++() {
			++idx;
			k = this->index_to_momentum(idx);
			return *this;
		}
		inline InnerIterator operator++(int) {
			InnerIterator tmp = *this;
			++(*this);
			return tmp;
		}

		inline InnerIterator& operator--() {
			--idx;
			k = this->index_to_momentum(idx);
			return *this;
		}
		inline InnerIterator operator--(int) {
			InnerIterator tmp = *this;
			--(*this);
			return tmp;
		}

		inline auto operator<=>(InnerIterator const& other) const = default;
	};

	class IEOMIterator {
		MomentumRanges const* const _parent;
		inline c_float index_to_momentum(int i) const {
			return _parent->index_to_momentum(i + 3 * _OUTER_DISC / 4);
		}
	public:
		c_float k{};
		int idx{};
		static inline int max_idx() { return _INNER_DISC + 45 * _OUTER_DISC / 100; }

		IEOMIterator(MomentumRanges const* const parent, int init = 0)
			: _parent(parent), k(this->index_to_momentum(init)), idx(init) {}

		inline c_float parent_step() const {
			if (k < _parent->INNER_K_MIN) return _parent->LOWER_STEP;
			if (k <= _parent->INNER_K_MAX) return _parent->INNER_STEP;
			return _parent->UPPER_STEP;
		}
		inline c_float max_k() const { return this->index_to_momentum(max_idx()); }
		inline c_float min_k() const { return this->index_to_momentum(0); }

		inline IEOMIterator& operator++() {
			++idx;
			k = this->index_to_momentum(idx);
			return *this;
		}
		inline IEOMIterator operator++(int) {
			IEOMIterator tmp = *this;
			++(*this);
			return tmp;
		}

		inline IEOMIterator& operator--() {
			--idx;
			k = this->index_to_momentum(idx);
			return *this;
		}
		inline IEOMIterator operator--(int) {
			IEOMIterator tmp = *this;
			--(*this);
			return tmp;
		}

		inline auto operator<=>(IEOMIterator const& other) const = default;
	};

	template<class T>
	concept is_momentum_iterator = std::same_as<T, MomentumIterator> || std::same_as<T, InnerIterator> || std::same_as<T, IEOMIterator>;

	template <class MomIt> requires is_momentum_iterator<MomIt>
	auto operator<=>(MomIt const& it, int i) { return it.idx <=> i; }
	template <class MomIt> requires is_momentum_iterator<MomIt>
	bool operator==(MomIt const& it, int i) { return it.idx == i; }
	template <class MomIt> requires is_momentum_iterator<MomIt>
	bool operator!=(MomIt const& it, int i) { return it.idx != i; }

	template <class MomIt> requires is_momentum_iterator<MomIt>
	std::ostream& operator<<(std::ostream& os, const MomIt& it) {
		os << it.idx << " -> " << it.k;
		return os;
	}
}
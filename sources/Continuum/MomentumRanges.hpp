#pragma once
#include "GlobalDefinitions.hpp"
#include <concepts>
#include <boost/math/quadrature/gauss.hpp>

namespace Continuum {
	struct MomentumRanges {
		static constexpr int n_gauss = 120;

		c_float K_MAX{};
		c_float K_MIN{};

		c_float INNER_K_MAX{};
		c_float INNER_K_MIN{};

		c_float LOWER_STEP{};
		c_float INNER_STEP{};
		c_float UPPER_STEP{};

		c_float const * K_F;

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

			if (begin <= INNER_K_MIN) {
				value += __integrate(func, begin, std::min(end, INNER_K_MIN));
				begin = INNER_K_MIN;
			}

			if (begin <= (*K_F) && end >= INNER_K_MIN) {
				value += __integrate(func, std::max(begin, INNER_K_MIN), std::min(end, (*K_F)));
				begin = (*K_F);
			}

			if (begin <= INNER_K_MAX && end >= (*K_F)) {
				value += __integrate(func, std::max(begin, (*K_F)), std::min(end, INNER_K_MAX));
				begin = INNER_K_MAX;
			}

			if (end >= INNER_K_MAX) {
				value += __integrate(func, std::max(begin, INNER_K_MAX), end);
			}

			return value;
		}

		template<class Function>
		inline auto integrate(const Function& func) const {
			return __integrate(func, K_MIN, INNER_K_MIN)
				+ __integrate(func, INNER_K_MIN, (*K_F))
				+ __integrate(func, (*K_F), INNER_K_MAX)
				+ __integrate(func, INNER_K_MAX, K_MAX);
		}

	private:
		template<class Function>
		inline auto __integrate(const Function& func, c_float begin, c_float end) const {
			if (is_zero(end - begin)) {
				return decltype(func(begin)){ };
			}
			return boost::math::quadrature::gauss<c_float, n_gauss>::integrate(func, begin, end);
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
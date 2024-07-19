#pragma once
#include "GlobalDefinitions.hpp"
#include <boost/math/quadrature/gauss.hpp>

namespace Continuum {
    class MomentumIterator;
    struct MomentumRanges {
        static constexpr int n_gauss = 60;

		c_float K_MAX{};
		c_float K_MIN{};

		c_float INNER_K_MAX{};
		c_float INNER_K_MIN{};

		c_float LOWER_STEP{};
		c_float INNER_STEP{};
        c_float UPPER_STEP{};

        c_float* K_F;

        MomentumRanges(c_float* k_F, const c_float omega_debye);

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
            if(is_zero(begin - end)) return value;

            if(begin <= INNER_K_MIN) {
                value += __integrate(func, begin, std::min(end, INNER_K_MIN));
                begin = INNER_K_MIN;
            }

            if(begin <= (*K_F) && end >= INNER_K_MIN) {
                value += __integrate(func, std::max(begin, INNER_K_MIN), std::min(end, (*K_F)));
                begin = (*K_F);
            }

            if(begin <= INNER_K_MAX && end >= (*K_F)){
                value += __integrate(func, std::max(begin, (*K_F)), std::min(end, INNER_K_MAX));
                begin = INNER_K_MAX;
            }

            if(end >= INNER_K_MAX) {
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
            if(is_zero(end - begin)) {
                return decltype(func(begin)){ };
            }
            return boost::math::quadrature::gauss<c_float, n_gauss>::integrate(func, begin, end);
        }
	};

    class MomentumIterator {
        MomentumRanges const * const _parent;
    public:
        c_float k{};
        int idx{};

        MomentumIterator(MomentumRanges const * const parent, int init = 0) 
            : _parent(parent), k(_parent->index_to_momentum(init)), idx(init) {}

        inline MomentumIterator& operator++() {
            ++idx;
            k = _parent->index_to_momentum(idx);
            return *this;
        }

        inline auto operator<=>(int other) { return idx <=> other; }
        inline auto operator==(int other) { return idx == other; }
        inline auto operator!=(int other) { return idx != other; }
    };

    class InnerIterator {
        MomentumRanges const * const _parent;
    public:
        c_float k{};
        int idx{};

        InnerIterator(MomentumRanges const * const parent, int init = 0) 
            : _parent(parent), k(_parent->index_to_momentum(init + _OUTER_DISC)), idx(init) {}

        inline InnerIterator& operator++() {
            ++idx;
            k = _parent->index_to_momentum(idx + _OUTER_DISC);
            return *this;
        }

        inline auto operator<=>(int other) { return idx <=> other; }
        inline auto operator==(int other) { return idx == other; }
        inline auto operator!=(int other) { return idx != other; }
    };
}
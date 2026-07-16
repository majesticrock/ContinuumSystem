#pragma once
#include "GlobalDefinitions.hpp"
#include "SCModel.hpp"

#include <mrock/symbolic_operators/TermLoader.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/utility/better_to_string.hpp>

#include <memory>
#include <map>

#ifndef _complex
#define _XP
#endif

#ifndef _XP
#include <mrock/iEoM/GeneralResolvent.hpp>
#define __iEoM_algorithm__ mrock::iEoM::GeneralResolvent<ModeHelper, c_complex>
#else
#include <mrock/iEoM/XPResolvent.hpp>
#define __iEoM_algorithm__ mrock::iEoM::XPResolvent<c_float, 10, false>
#endif

namespace Continuum {
	typedef __iEoM_algorithm__ iEoM_algorithm;

	class ModeHelper : public iEoM_algorithm
	{
	private:
		using m_iterator = InnerIterator;

		c_float compute_momentum(mrock::symbolic_operators::Momentum const& momentum, c_float k, c_float l, c_float q = 0) const;
		c_complex get_expectation_value(mrock::symbolic_operators::WickOperator const& op, c_float momentum) const;

		c_complex compute_phonon_sum(const mrock::symbolic_operators::WickTerm& term, c_float k) const;
		c_complex compute_em_sum(const mrock::symbolic_operators::WickTerm& term, c_float k) const;
	protected:
		constexpr static int hermitian_size = 2;
		constexpr static int antihermitian_size = 1;
		constexpr static int number_of_basis_terms = hermitian_size + antihermitian_size;

		static int hermitian_discretization;
		static int antihermitian_discretization;
		static int total_matrix_size;

		std::unique_ptr<SCModel> model;
		mrock::symbolic_operators::TermLoader wicks;

		void create_starting_states();
		void fill_matrices();
		void fill_M();

		void fill_block_M(int i, int j);
		void fill_block_N(int i, int j);

		c_complex computeTerm(const mrock::symbolic_operators::WickTerm& term, c_float k, c_float l) const;
	public:
		std::vector<c_float> continuum_boundaries() const;

		SCModel& getModel() {
			return *model;
		};
		const SCModel& getModel() const {
			return *model;
		};
		ModeHelper(ModelInitializer const& init);
	};
}

#undef __iEoM_algorithm__
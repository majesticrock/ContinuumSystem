#pragma once
#include "GlobalDefinitions.hpp"
#include "../../../FermionCommute/sources/TermLoader.hpp"
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include "../../../Utility/sources/better_to_string.hpp"
#include "../../../Utility/sources/Numerics/iEoM/GeneralResolvent.hpp"
#include "../../../Utility/sources/Numerics/iEoM/XPResolvent.hpp"
#include "SCModel.hpp"
#include <memory>
#include <map>

namespace Continuum {
	class ModeHelper : public Utility::Numerics::iEoM::XPResolvent<ModeHelper, c_float>
	{
	private:
		friend struct Utility::Numerics::iEoM::XPResolvent<ModeHelper, c_float>;
		using _parent = Utility::Numerics::iEoM::XPResolvent<ModeHelper, c_float>;

		std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> expectation_values;

		c_float compute_momentum(SymbolicOperators::Momentum const& momentum, c_float k, c_float l, c_float q = 0) const;
		c_complex get_expectation_value(SymbolicOperators::WickOperator const& op, c_float momentum) const;

		c_complex compute_phonon_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const;
		c_complex compute_em_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const;
	protected:
		SymbolicOperators::TermLoader wicks;
		//size_t TOTAL_BASIS{};
		constexpr static int hermitian_size = 3;
		constexpr static int antihermitian_size = 1;
		constexpr static int number_of_basis_terms = hermitian_size + antihermitian_size;

		static int hermitian_discretization;
		static int antihermitian_discretization;

		std::unique_ptr<SCModel> model;

		void createStartingStates();
		void fillMatrices();
		void fill_M();

		void fill_block_M(int i, int j);
		void fill_block_N(int i, int j);

		c_complex computeTerm(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const;
	public:
		const SCModel& getModel() const {
			return *model;
		};
		ModeHelper(Utility::InputFileReader& input);
	};
}
#pragma once
#include "GlobalDefinitions.hpp"
#include "../../../../FermionCommute/sources/TermLoader.hpp"
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include "../../../Utility/sources/better_to_string.hpp"
#include "../../../Utility/sources/Numerics/iEoM/GeneralResolvent.hpp"
#include "SCModel.hpp"
#include <memory>

namespace Continuum {
	class ModeHelper : public Utility::Numerics::iEoM::GeneralResolvent<ModeHelper, c_float>
	{
		friend class Utility::Numerics::iEoM::GeneralResolvent<ModeHelper, c_float>;
		using _parent = Utility::Numerics::iEoM::GeneralResolvent<ModeHelper, c_float>;
	protected:
		SymbolicOperators::TermLoader wicks;
		size_t TOTAL_BASIS{};

		std::unique_ptr<SCModel> model;

		int number_of_basis_terms{};
		int start_basis_at{};

		void createStartingStates();
		void fillMatrices();
		void fill_M();

		void fill_block_M(int i, int j);
		void fill_block_N(int i, int j);

		c_complex computeTerm(const SymbolicOperators::WickTerm&, c_float k, c_float l);
	public:
		ModeHelper(Utility::InputFileReader& input);
	};
}
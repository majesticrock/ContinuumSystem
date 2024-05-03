#include "ModeHelper.hpp"
#include <memory>
#include <cassert>
#include "../../../Utility/sources/Numerics/Interpolation.hpp"
#include "../../../Utility/sources/Numerics/TrapezoidalRule.hpp"
#include "../../../Utility/sources/Selfconsistency/IterativeSolver.hpp"

namespace Continuum {
	c_float ModeHelper::compute_momentum(SymbolicOperators::Momentum const& momentum, c_float k, c_float l, c_float q /*=0*/) const
	{
		// For now we do not encounter sums of momenta, e.g., k+l
		// Should that occur, we need to implement a logic for that here
		assert(momentum.momentum_list.size() == 1U);
		switch (momentum.momentum_list.front().second) {
		case 'k':
			return momentum.momentum_list.front().first * k;
		case 'l':
			return momentum.momentum_list.front().first * l;
		case 'q':
			return momentum.momentum_list.front().first * q;
		default:
			throw std::runtime_error("Momentum not recognized! " + momentum.momentum_list.front().second);
		}
	}
	
	c_complex ModeHelper::get_expectation_value(SymbolicOperators::WickOperator const& op, c_float momentum) const
	{
		const int index = momentum_to_index(momentum);
		return Utility::Numerics::linearly_interpolate(momentum, index_to_momentum(index), index_to_momentum(index + 1),
			expectation_values.at(op.type)[index], expectation_values.at(op.type)[index + 1]);
	}

	void ModeHelper::createStartingStates()
	{
		starting_states.resize(2, _parent::Vector::Zero(TOTAL_BASIS));
		std::fill(starting_states[0].begin(), starting_states[0].begin() + DISCRETIZATION, 1. / sqrt(DISCRETIZATION));
		std::fill(starting_states[1].begin() + 3 * DISCRETIZATION, starting_states[1].end(), 1. / sqrt(DISCRETIZATION));
	}

	void ModeHelper::fillMatrices()
	{
		M = _parent::Matrix::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = _parent::Matrix::Zero(TOTAL_BASIS, TOTAL_BASIS);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				fill_block_N(i, j);
				fill_block_M(i, j);
			}
		}
	}

	void ModeHelper::fill_M()
	{
		M = _parent::Matrix::Zero(TOTAL_BASIS, TOTAL_BASIS);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				fill_block_M(i, j);
			}
		}
	}

	void ModeHelper::fill_block_M(int i, int j)
	{
		for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < DISCRETIZATION; ++k) {
				if (!term.delta_momenta.empty()) {
					// only k=l and k=-l should occur. Additionally, only the magntiude should matter
					M(i + k * number_of_basis_terms, j + k * number_of_basis_terms) += computeTerm(term, k, k).real();
				}
				else {
					for (int l = 0; l < DISCRETIZATION; ++l) {
						M(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
							+= computeTerm(term, k, l).real();
					}
				}
			}
		}
	}

	void ModeHelper::fill_block_N(int i, int j)
	{
		for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < DISCRETIZATION; ++k) {
				if (!term.delta_momenta.empty()) {
					// only k=l and k=-l should occur. Additionally, only the magntiude should matter
					M(i + k * number_of_basis_terms, j + k * number_of_basis_terms) += computeTerm(term, k, k).real();
				}
				else {
					for (int l = 0; l < DISCRETIZATION; ++l) {
						M(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
							+= computeTerm(term, k, l).real();
					}
				}
			}
		}
	}

	c_complex ModeHelper::computeTerm(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
		if (term.sums.momenta.empty()) {
			c_complex value;
			if (term.coefficients.empty()) {
				value = term.multiplicity;
			}
			else {
				const SymbolicOperators::Coefficient* coeff_ptr = &term.coefficients.front();
				if (coeff_ptr->momenta.size() == 2U) {
					value = term.multiplicity * model->computeCoefficient(*coeff_ptr,
						compute_momentum(coeff_ptr->momenta[0], k, l), compute_momentum(coeff_ptr->momenta[1], k, l));
				}
				else if (coeff_ptr->momenta.size() == 1U) {
					value = term.multiplicity * model->computeCoefficient(*coeff_ptr,
						compute_momentum(coeff_ptr->momenta[0], k, l));
				}
				else {
					throw std::runtime_error("Number of momenta of coefficient is not handled! "
						+ std::to_string(coeff_ptr->momenta.size()));
				}
			}
			if (term.operators.empty()) return value;

			for (const auto& op : term.operators) {
				value *= this->get_expectation_value(op, this->compute_momentum(op.momentum, k, l));
			}
			return value * static_cast<c_float>(term.sums.spins.size() + 1U);
		}
		// For now, sums can only ever occur with U
		assert(term.coefficients.size() == 1U && term.coefficients.front().name == "U");

		auto integrand = [&](c_float q) {
			c_complex value{ this->get_expectation_value(term.operators.front(), 
				this->compute_momentum(term.operators.front().momentum, k, l, q))};
			
			for (auto it = term.operators.begin() + 1; it != term.operators.end(); ++it) {
				value *= this->get_expectation_value(*it, this->compute_momentum(it->momentum, k, l, q));
			}
			return value;
			};

		return term.multiplicity * this->model->computeCoefficient(term.coefficients.front(), c_float{}) *
			Utility::Numerics::Integration::trapezoidal_rule(integrand, model->u_lower_bound(k), model->u_upper_bound(k), DISCRETIZATION);
	}

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
		: _parent(this, SQRT_PRECISION<c_float>),
		number_of_basis_terms{ input.getInt("number_of_basis_terms") }, start_basis_at{ input.getInt("start_basis_at") }
	{
		model = std::make_unique<SCModel>(ModelInitializer(input));
		TOTAL_BASIS = DISCRETIZATION * this->number_of_basis_terms;
		wicks.load("../commutators/continuum/", true, this->number_of_basis_terms, start_basis_at);

		Utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(model.get(), &model->Delta);
		solver.computePhases();

		this->expectation_values = model->get_expectation_values();
	}
}
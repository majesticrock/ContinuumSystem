#include "ModeHelper.hpp"
#include <memory>
#include <cassert>
#include <Utility/Numerics/Interpolation.hpp>
#include <Utility/Numerics/Integration/TrapezoidalRule.hpp>
#include <Utility/Selfconsistency/IterativeSolver.hpp>
#include <Utility/Selfconsistency/BroydenSolver.hpp>
#include <Utility/ConstexprPower.hpp>

#include <boost/math/quadrature/gauss.hpp>

#define ieom_diag(k) k * k / (2 * PI * PI)
#define ieom_offdiag(k, l) k * k * l * l * model->momentumRanges.INNER_STEP / (4 * PI * PI * PI * PI)

#ifdef _complex
#define __conj(z) std::conj(z)
#else
#define __conj(z) z
#endif

namespace Continuum {
	c_float ModeHelper::compute_momentum(SymbolicOperators::Momentum const& momentum, c_float k, c_float l, c_float q /*=0*/) const
	{
		c_float momentum_value{};
		for (const auto& momentum_pair : momentum.momentum_list) {
			switch (momentum_pair.second) {
			case 'k':
				momentum_value += momentum_pair.first * k;
				break;
			case 'l':
				momentum_value += momentum_pair.first * l;
				break;
			case 'q':
				momentum_value += momentum_pair.first * q;
				break;
			default:
				throw std::runtime_error("Momentum not recognized! " + momentum_pair.second);
			}
		}
		return momentum_value;
	}

	c_complex ModeHelper::get_expectation_value(SymbolicOperators::WickOperator const& op, c_float momentum) const
	{
		if (op.type == SymbolicOperators::Number_Type) {
			return model->occupation(momentum);
		}
		else if (op.type == SymbolicOperators::SC_Type) {
			if (op.isDaggered) {
				return __conj(model->sc_expectation_value(momentum));
			}
			return model->sc_expectation_value(momentum);
		}
		throw std::runtime_error("Expectation value not recognized!");
	}

	void ModeHelper::createStartingStates()
	{
#ifdef _complex
		starting_states.resize(2, _parent::Vector::Zero(total_matrix_size));
		std::fill(starting_states[0].begin(), starting_states[0].begin() + MODE_DISC, sqrt(model->momentumRanges.INNER_STEP));
		std::fill(starting_states[1].begin() + 3 * MODE_DISC, starting_states[1].end(), sqrt(model->momentumRanges.INNER_STEP));
#else
		starting_states.push_back({ _parent::Vector::Zero(antihermitian_discretization), _parent::Vector::Zero(hermitian_discretization), "SC" });
		std::fill(starting_states[0][0].begin(), starting_states[0][0].begin() + MODE_DISC, sqrt(model->momentumRanges.INNER_STEP));
		std::fill(starting_states[0][1].begin(), starting_states[0][1].begin() + MODE_DISC, sqrt(model->momentumRanges.INNER_STEP));
#endif
	}

	void ModeHelper::fillMatrices()
	{
#ifdef _complex
		M.setZero(total_matrix_size, total_matrix_size);
		N.setZero(total_matrix_size, total_matrix_size);
#else
		K_plus.setZero(hermitian_discretization, hermitian_discretization);
		K_minus.setZero(antihermitian_discretization, antihermitian_discretization);
		L.setZero(hermitian_discretization, antihermitian_discretization);
#endif

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
#ifdef _complex
				fill_block_M(i, j);
				fill_block_N(i, j);
#else
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
				// N only contains offdiagonal blocks
				else if (i < hermitian_size && j >= hermitian_size) {
					fill_block_N(i, j);
				}
#endif
			}
		}

#ifndef _complex
		std::cout << "||K_+ - K_+^+|| = " << (K_plus - K_plus.adjoint()).norm() << std::endl;
		std::cout << "||K_- - K_-^+|| = " << (K_minus - K_minus.adjoint()).norm() << std::endl;
		/* constexpr int block_size = 2;
		//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
		for(int i = MODE_DISC / 2 - 10; i < MODE_DISC / 2 + 10; ++i) {
			//solver.compute(K_minus.block(i, i, block_size, block_size), false);
			//for(int j = 0; j < block_size; ++j){
				//if(solver.eigenvalues()(j) < 0) {
				//	std::cerr << i << ", " << j << " -> " << solver.eigenvalues()(j) << std::endl;
					double k = model->momentumRanges.index_to_momentum(_OUTER_DISC + i);
					double l = model->momentumRanges.index_to_momentum(_OUTER_DISC + i + 1);
					std::cout << i << "  " << k / model->fermi_wavevector << ": \n" << K_minus.block(i, i, block_size, block_size) << std::endl;

					double V0 = model->coulomb_scaling / (PhysicalConstants::vacuum_permitivity * (_screening * _screening));
					double V = model->coulomb_scaling / (PhysicalConstants::vacuum_permitivity * ((k - l) * (k - l) + _screening * _screening));
					std::cout << 8 * model->sc_expectation_value(k) * model->sc_expectation_value(l) * (2 * V0 - V)
						+ 2 * V * (1 - 2 * model->occupation(k)) * (1 - 2 * model->occupation(l)) << std::endl;
					std::cout << model->occupation(k) << "  " << model->sc_expectation_value(k) << std::endl;
					std::cout << "##################\n";
				//}
			//}
		} */
#else
		std::cout << "||M - M^+|| = " << (M - M.adjoint()).norm() << std::endl;
		std::cout << "||N - N^+|| = " << (N - N.adjoint()).norm() << std::endl;
#endif
	}

	void ModeHelper::fill_M()
	{
#ifdef _complex
		M.setZero(total_matrix_size, total_matrix_size);
#else
		K_plus.setZero(hermitian_discretization, hermitian_discretization);
		K_minus.setZero(antihermitian_discretization, antihermitian_discretization);
#endif

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
#ifndef _complex
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size))
#endif
					fill_block_M(i, j);
			}
		}
	}

	void ModeHelper::fill_block_M(int i, int j)
	{
		for (InnerIterator it(&model->momentumRanges); it < _INNER_DISC; ++it) {
			for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
				if (!term.delta_momenta.empty()) {
					if (term.sums.momenta.empty()) {
						if (term.coefficients.front().name == "g" || term.coefficients.front().name == "V") {
							// These kinds of terms scale as 1/V -> 0
							continue;
						}
					}
					M(i * MODE_DISC + it.idx, j * MODE_DISC + it.idx) += ieom_diag(it.k) * computeTerm(term, it.k, it.k);
				}
				else {
					for (InnerIterator jt(&model->momentumRanges); jt < _INNER_DISC; ++jt) {
						M(i * MODE_DISC + it.idx, j * MODE_DISC + jt.idx) += ieom_offdiag(it.k, jt.k) * computeTerm(term, it.k, jt.k);
					}
				}
			}
		}
	}

	void ModeHelper::fill_block_N(int i, int j)
	{
		for (InnerIterator it(&model->momentumRanges); it < _INNER_DISC; ++it) {
			for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
				if (!term.delta_momenta.empty()) {
					// only k=l and k=-l should occur. Additionally, only the magntitude should matter
					N(i * MODE_DISC + it.idx, j * MODE_DISC + it.idx) += ieom_diag(it.k) * computeTerm(term, it.k, it.k);
				}
				else {
					for (InnerIterator jt(&model->momentumRanges); jt < _INNER_DISC; ++jt) {
						N(i * MODE_DISC + it.idx, j * MODE_DISC + jt.idx) += ieom_offdiag(it.k, jt.k) * computeTerm(term, it.k, jt.k);
					}
				}
			}
		}
	}

	c_complex ModeHelper::compute_phonon_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
		const int q_dependend = term.whichOperatorDependsOn('q');
		SymbolicOperators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		SymbolicOperators::WickOperator const* const other_op = term.isBilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		c_complex value{};
		if (summed_op->type == SymbolicOperators::Number_Type) {
			value = -static_cast<c_float>(term.multiplicity) * model->integral_phonon(model->occupation, k);
		}
		else {
			value = -static_cast<c_float>(term.multiplicity) * model->integral_phonon(model->sc_expectation_value, k);
#ifdef _complex
			if (summed_op->isDaggered) value = std::conj(value);
#endif
		}
		if (other_op) {
			value *= this->get_expectation_value(*other_op, k);
		}
		return value;
	}

	c_complex ModeHelper::compute_em_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
		const int q_dependend = term.whichOperatorDependsOn('q');
		SymbolicOperators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		SymbolicOperators::WickOperator const* const other_op = term.isBilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		c_complex value{};
		if (summed_op->type == SymbolicOperators::Number_Type) {
			value = -static_cast<c_float>(term.multiplicity) * (model->fock_energy(k) + model->interpolate_delta_n(k));
		}
		else {
			value = static_cast<c_float>(term.multiplicity) * model->integral_screening(model->sc_expectation_value, k);
#ifdef _complex
			if (summed_op->isDaggered) value = std::conj(value);
#endif
		}
		if (other_op) {
			value *= this->get_expectation_value(*other_op, k);
		}
		return value;
	}

	c_complex ModeHelper::computeTerm(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
		if (term.sums.momenta.empty()) {
			c_complex value{ static_cast<c_float>(term.multiplicity) };

			if (!term.coefficients.empty()) {
				const SymbolicOperators::Coefficient* coeff_ptr = &term.coefficients.front();
				if (coeff_ptr->momenta.size() == 2U) {
					value *= model->computeCoefficient(*coeff_ptr, compute_momentum(coeff_ptr->momenta[0], k, l), compute_momentum(coeff_ptr->momenta[1], k, l));
				}
				else if (coeff_ptr->momenta.size() == 1U) {
					value *= model->computeCoefficient(*coeff_ptr, k, l);
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
		assert(term.coefficients.size() == 1U);

		if (term.coefficients.front().name == "g")
		{
			return compute_phonon_sum(term, k, l);
		}
		else if (term.coefficients.front().name == "V") {
			return compute_em_sum(term, k, l);
		}
		throw std::runtime_error("Something went wrong while computing terms...");
	}

	int ModeHelper::hermitian_discretization = 0;
	int ModeHelper::antihermitian_discretization = 0;
	int ModeHelper::total_matrix_size = 0;

	ModeHelper::ModeHelper(ModelInitializer const& init)
		: _parent(this, PRECISION,
#ifndef _complex
			MODE_DISC* hermitian_size, MODE_DISC* antihermitian_size, false,
#endif
			false)
	{
		hermitian_discretization = MODE_DISC * hermitian_size;
		antihermitian_discretization = MODE_DISC * antihermitian_size;
		total_matrix_size = MODE_DISC * number_of_basis_terms;

		model = std::make_unique<SCModel>(init);
		wicks.load("../commutators/continuum/", true, number_of_basis_terms, 0);

		//auto solver = Utility::Selfconsistency::make_iterative<c_complex>(model.get(), &model->Delta);
		auto solver = Utility::Selfconsistency::make_broyden<c_complex>(model.get(), &model->Delta, 200);
		solver.compute(true);
	}
}
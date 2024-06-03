#include "ModeHelper.hpp"
#include <memory>
#include <cassert>
#include "../../../Utility/sources/Numerics/Interpolation.hpp"
#include "../../../Utility/sources/Numerics/Integration/TrapezoidalRule.hpp"

#include "../../../Utility/sources/Selfconsistency/IterativeSolver.hpp"
#include "../../../Utility/sources/Selfconsistency/BroydenSolver.hpp"

#include <boost/math/quadrature/gauss.hpp>

#define ieom_diag(k) 4. * PI * k * k
#define ieom_offdiag(k, l) (2. / PI) * k * k * l * l * model->STEP

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
		if (op.type == SymbolicOperators::Number_Type) {
			return model->occupation(momentum);
		}
		else if (op.type == SymbolicOperators::SC_Type) {
			return model->sc_expectation_value(momentum);
		}
		assert(false && "Expectation value not recognized!");
	}

	void ModeHelper::createStartingStates()
	{
		starting_states.resize(1, { _parent::Vector::Zero(antihermitian_discretization), _parent::Vector::Zero(hermitian_discretization) });
		std::fill(starting_states[0][0].begin(), starting_states[0][0].begin() + DISCRETIZATION, sqrt(model->STEP));
		std::fill(starting_states[0][1].begin(), starting_states[0][1].begin() + DISCRETIZATION, sqrt(model->STEP));
	}

	void ModeHelper::fillMatrices()
	{
		K_plus.setZero(hermitian_discretization, hermitian_discretization);
		K_minus.setZero(antihermitian_discretization, antihermitian_discretization);
		L.setZero(hermitian_discretization, antihermitian_discretization);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
				// N only contains offdiagonal blocks
				else if (i < hermitian_size && j >= hermitian_size) {
					fill_block_N(i, j);
				}
			}
		}

		//std::cout << "||K_+ - K_+^+|| = " << (K_plus - K_plus.adjoint()).norm()   << " && ||K_+|| = " << K_plus.norm() << std::endl;
		//std::cout << "||K_- - K_-^+|| = " << (K_minus - K_minus.adjoint()).norm() << " && ||K_-|| = " << K_minus.norm() << std::endl;
		//for (int i = 0; i < K_plus.diagonal().real().size(); ++i) {
		//	if (K_plus.diagonal().real()(i) < -PRECISION) {
		//		std::cout << i << "+: " << K_plus.diagonal().real()(i) << "\n";
		//	}
		//}
		//for (int i = 0; i < K_minus.diagonal().real().size(); ++i) {
		//	if (K_minus.diagonal().real()(i) < -PRECISION) {
		//		std::cout << i << "-: " << K_minus.diagonal().real()(i) << "\n";
		//	}
		//}
	}

	void ModeHelper::fill_M()
	{
		K_plus = _parent::Matrix::Zero(hermitian_discretization, hermitian_discretization);
		K_minus = _parent::Matrix::Zero(antihermitian_discretization, antihermitian_discretization);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
			}
		}
	}

	void ModeHelper::fill_block_M(int i, int j)
	{
		for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
			for (int k_idx = 0; k_idx < DISCRETIZATION; ++k_idx) {
				const c_float k = this->model->index_to_momentum(k_idx);
				if (!term.delta_momenta.empty()) {
					//if (i == j && i == 0) {
					//	if(k_idx < 20)
					//		std::cout << "k=" <<k <<"\t\t" << term << " = " << computeTerm(term, k, k).real() << std::endl;
					//}
					// only k=l and k=-l should occur. Additionally, only the magntiude should matter

					if (term.sums.momenta.empty() && term.coefficients.front().name == "g") {
						// These kinds of terms scale as 1/N -> 0
						continue;
					}
					if (i < hermitian_size) {
						K_plus(i * DISCRETIZATION + k_idx, j * DISCRETIZATION + k_idx)
							+= ieom_diag(k) * std::real(computeTerm(term, k, k));
					}
					else {
						K_minus((i - hermitian_size) * DISCRETIZATION + k_idx, (j - hermitian_size) * DISCRETIZATION + k_idx)
							+= ieom_diag(k) * std::real(computeTerm(term, k, k));
					}
				}
				else {
					for (int l_idx = 0; l_idx < DISCRETIZATION; ++l_idx) {
						const c_float l = this->model->index_to_momentum(l_idx);
						if (i < hermitian_size) {
							K_plus(i * DISCRETIZATION + k_idx, j * DISCRETIZATION + l_idx)
								+= ieom_offdiag(k, l) * std::real(computeTerm(term, k, l));
						}
						else {
							K_minus((i - hermitian_size) * DISCRETIZATION + k_idx, (j - hermitian_size) * DISCRETIZATION + l_idx)
								+= ieom_offdiag(k, l) * std::real(computeTerm(term, k, l));
						}
					}
				}
			}
		}
	}

	void ModeHelper::fill_block_N(int i, int j)
	{
		for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
			for (int k_idx = 0; k_idx < DISCRETIZATION; ++k_idx) {
				const c_float k = this->model->index_to_momentum(k_idx);
				if (!term.delta_momenta.empty()) {
					// only k=l and k=-l should occur. Additionally, only the magntiude should matter
					L(i * DISCRETIZATION + k_idx, (j - hermitian_size) * DISCRETIZATION + k_idx)
						+= ieom_diag(k) * std::real(computeTerm(term, k, k));
				}
				else {
					for (int l_idx = 0; l_idx < DISCRETIZATION; ++l_idx) {
						const c_float l = this->model->index_to_momentum(l_idx);
						L(i * DISCRETIZATION + k_idx, (j - hermitian_size) * DISCRETIZATION + l_idx)
							+= ieom_offdiag(k, l) * std::real(computeTerm(term, k, l));
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
				value = static_cast<c_float>(term.multiplicity);
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
		assert(term.coefficients.size() == 1U && term.coefficients.front().name == "g");

		auto integrand = [&](c_float q) {
			c_complex value{ this->get_expectation_value(term.operators.front(),
				this->compute_momentum(term.operators.front().momentum, k, l, q)) };

			for (auto it = term.operators.begin() + 1; it != term.operators.end(); ++it) {
				value *= this->get_expectation_value(*it, this->compute_momentum(it->momentum, k, l, q));
			}
			return q * q * value;
			};

#ifdef approximate_theta
		if (k > this->model->g_upper_bound(k)) return 0;
#endif
		//c_float error;
		return (term.multiplicity / (2.0 * PI * PI)) * model->computeCoefficient(term.coefficients.front(), model->fermi_wavevector)
			//* Utility::Numerics::Integration::trapezoidal_rule(integrand, model->g_lower_bound(k), model->g_upper_bound(k), DISCRETIZATION);
			* boost::math::quadrature::gauss<double, 30>::integrate(integrand, model->g_lower_bound(k), model->g_upper_bound(k));
	}

	int ModeHelper::hermitian_discretization = 0;
	int ModeHelper::antihermitian_discretization = 0;

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
		: _parent(this, 1e-5, false, false) //SQRT_PRECISION
	{
		hermitian_discretization = DISCRETIZATION * hermitian_size;
		antihermitian_discretization = DISCRETIZATION * antihermitian_size;

		model = std::make_unique<SCModel>(ModelInitializer(input));
		wicks.load("../commutators/continuum/", true, number_of_basis_terms, 0);

		//Utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(model.get(), &model->Delta);
		Utility::Selfconsistency::BroydenSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(model.get(), &model->Delta, 200);
		solver.compute(true);

		this->expectation_values = model->get_expectation_values();
	}
}
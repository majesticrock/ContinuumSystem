#include "ModeHelper.hpp"
#include <memory>
#include <cassert>
#include <Utility/Numerics/Interpolation.hpp>
#include <Utility/Numerics/Integration/TrapezoidalRule.hpp>
#include <Utility/Selfconsistency/IterativeSolver.hpp>
#include <Utility/Selfconsistency/BroydenSolver.hpp>
#include <Utility/ConstexprPower.hpp>

#include <boost/math/quadrature/gauss.hpp>

#define ieom_diag(k) 4. * PI * k * k
#define ieom_offdiag(k, l) (2. / PI) * k * k * l * l * model->STEP

#ifdef _complex
#define __conj(z) std::conj(z)
#else
#define __conj(z) z
#endif

namespace Continuum {
	c_float ModeHelper::compute_momentum(SymbolicOperators::Momentum const& momentum, c_float k, c_float l, c_float q /*=0*/) const
	{
		c_float momentum_value{};
		for(const auto& momentum_pair : momentum.momentum_list) {
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
			if(op.isDaggered) {
				return __conj(model->sc_expectation_value(momentum));
			}
			return model->sc_expectation_value(momentum);
		}
		assert(false && "Expectation value not recognized!");
	}

	void ModeHelper::createStartingStates()
	{
#ifdef _complex
		starting_states.resize(2, _parent::Vector::Zero(total_matrix_size));
		std::fill(starting_states[0].begin(), starting_states[0].begin() + DISCRETIZATION, sqrt(model->STEP));
		std::fill(starting_states[1].begin() + 3 * DISCRETIZATION, starting_states[1].end(), sqrt(model->STEP));
#else
		starting_states.resize(1, { _parent::Vector::Zero(antihermitian_discretization), _parent::Vector::Zero(hermitian_discretization) });
		std::fill(starting_states[0][0].begin(), starting_states[0][0].begin() + DISCRETIZATION, sqrt(model->STEP));
		std::fill(starting_states[0][1].begin(), starting_states[0][1].begin() + DISCRETIZATION, sqrt(model->STEP));
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
		for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
			for (int k_idx = 0; k_idx < DISCRETIZATION; ++k_idx) {
				const c_float k = this->model->index_to_momentum(k_idx);
				if (!term.delta_momenta.empty()) {
					if (term.sums.momenta.empty()) {
						if (term.coefficients.front().name == "g" || term.coefficients.front().name == "V") {
							// These kinds of terms scale as 1/N -> 0
							continue;
						}
					}
					M(i * DISCRETIZATION + k_idx, j * DISCRETIZATION + k_idx)
							+= ieom_diag(k) * computeTerm(term, k, k);
				}
				else {
					for (int l_idx = 0; l_idx < DISCRETIZATION; ++l_idx) {
						const c_float l = this->model->index_to_momentum(l_idx);
						M(i * DISCRETIZATION + k_idx, j * DISCRETIZATION + l_idx)
								+= ieom_offdiag(k, l) * computeTerm(term, k, l);
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
					// only k=l and k=-l should occur. Additionally, only the magntitude should matter
					N(i * DISCRETIZATION + k_idx, j * DISCRETIZATION + k_idx)
						+= ieom_diag(k) * computeTerm(term, k, k);
				}
				else {
					for (int l_idx = 0; l_idx < DISCRETIZATION; ++l_idx) {
						const c_float l = this->model->index_to_momentum(l_idx);
						N(i * DISCRETIZATION + k_idx, j * DISCRETIZATION + l_idx)
							+= ieom_offdiag(k, l) * computeTerm(term, k, l);
					}
				}
			}
		}
	}

	c_complex ModeHelper::compute_phonon_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
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
		return (term.multiplicity / (2.0 * PI * PI)) * model->computeCoefficient(term.coefficients.front(), model->fermi_wavevector)
			* boost::math::quadrature::gauss<double, 30>::integrate(integrand, model->g_lower_bound(k), model->g_upper_bound(k));
	}

	c_complex ModeHelper::compute_em_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
		if(term.coefficients.front().momenta[0].momentum_list.size() == 1)
		{ // V(k, k-q)
			auto integrand = [&](c_float q) {
				c_complex value{ this->get_expectation_value(term.operators.front(),
				this->compute_momentum(term.operators.front().momentum, k, l, q)) };

				for (auto it = term.operators.begin() + 1; it != term.operators.end(); ++it) {
					value *= this->get_expectation_value(*it, this->compute_momentum(it->momentum, k, l, q));
				}
				const c_float weight = Utility::constexprPower<2>(std::min(k + q, model->g_upper_bound(q)))
					- Utility::constexprPower<2>(std::max(std::abs(k - q), model->g_lower_bound(q)));
				return (weight / q) * value;
				};

			return 0.5 * term.multiplicity * PhysicalConstants::em_factor
				* boost::math::quadrature::gauss<double, 30>::integrate(integrand, model->K_MIN, model->K_MAX);
		} 
		else if(term.coefficients.front().momenta[0].momentum_list.size() == 2) 
		{ // V(k-q, k)
			auto integrand = [&](c_float q) {
				c_complex value{ this->get_expectation_value(term.operators.front(),
				this->compute_momentum(term.operators.front().momentum, k, l, q)) };

				for (auto it = term.operators.begin() + 1; it != term.operators.end(); ++it) {
					value *= this->get_expectation_value(*it, this->compute_momentum(it->momentum, k, l, q));
				}
				const c_float weight = std::log(std::min(k + q, model->g_upper_bound(q)) / std::max(std::abs(k - q), model->g_lower_bound(q)));
				return (q * weight) * value;
				};

			return term.multiplicity * PhysicalConstants::em_factor
				* boost::math::quadrature::gauss<double, 30>::integrate(integrand, model->K_MIN, model->K_MAX);
		}
		throw std::runtime_error("Something went wrong while computing an em sum...");
	}

	c_complex ModeHelper::computeTerm(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
#ifndef _use_coulomb
		if(!term.coefficients.empty()){
			if(term.coefficients.front().name == "V") return c_complex{};
		}
#endif
		if (term.sums.momenta.empty()) {
			c_complex value { static_cast<c_float>(term.multiplicity) };
			
			if (!term.coefficients.empty()) {
				const SymbolicOperators::Coefficient* coeff_ptr = &term.coefficients.front();
				if (coeff_ptr->momenta.size() == 2U) {
					value *= model->computeCoefficient(*coeff_ptr,
						compute_momentum(coeff_ptr->momenta[0], k, l), compute_momentum(coeff_ptr->momenta[1], k, l));
				}
				else if (coeff_ptr->momenta.size() == 1U) {
					value *= model->computeCoefficient(*coeff_ptr,
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
		assert(term.coefficients.size() == 1U);

		if(term.coefficients.front().name == "g")
		{
			return compute_phonon_sum(term, k, l);
		}
		else if(term.coefficients.front().name == "V"){
			return compute_em_sum(term, k, l);
		}
		throw std::runtime_error("Something went wrong while computing terms...");
	}

	int ModeHelper::hermitian_discretization = 0;
	int ModeHelper::antihermitian_discretization = 0;
	int ModeHelper::total_matrix_size = 0;

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
		: _parent(this, SQRT_PRECISION, //SQRT_PRECISION
#ifndef _complex
		DISCRETIZATION * hermitian_size, DISCRETIZATION * antihermitian_size, false, 
#endif
		false)
	{
		hermitian_discretization = DISCRETIZATION * hermitian_size;
		antihermitian_discretization = DISCRETIZATION * antihermitian_size;
		total_matrix_size = DISCRETIZATION * number_of_basis_terms;

		model = std::make_unique<SCModel>(ModelInitializer(input));
		wicks.load("../commutators/continuum/", true, number_of_basis_terms, 0);

		Utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(model.get(), &model->Delta);
		//Utility::Selfconsistency::BroydenSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(model.get(), &model->Delta, 200);
		solver.compute(true, 1500);

		this->expectation_values = model->get_expectation_values();
	}
}
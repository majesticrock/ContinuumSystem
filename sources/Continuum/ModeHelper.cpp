#include "ModeHelper.hpp"
#include <memory>
#include <cassert>
#include <mrock/utility/Numerics/Interpolation.hpp>
#include <mrock/utility/Numerics/Integration/TrapezoidalRule.hpp>
#include <mrock/utility/Selfconsistency/IterativeSolver.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/ConstexprPower.hpp>
#include <mrock/utility/Numerics/Minimization/Bisection.hpp>
#include <boost/math/quadrature/gauss.hpp>

#define ieom_diag(k) k * k * (DISCRETIZATION * it.parent_step() / (2 * PI * PI))
#define ieom_offdiag(k, l) k * k * l * l * (DISCRETIZATION * it.parent_step() * jt.parent_step() / (4 * PI * PI * PI * PI))

#ifdef _complex
#define __conj(z) std::conj(z)
#else
#define __conj(z) z
#endif

namespace Continuum {
	c_float ModeHelper::compute_momentum(mrock::symbolic_operators::Momentum const& momentum, c_float k, c_float l, c_float q /*=0*/) const
	{
		c_float momentum_value{};
		for (const auto& momentum_symbol : momentum.momentum_list) {
			switch (static_cast<char>(momentum_symbol.name)) {
			case 'k':
				momentum_value += momentum_symbol.factor * k;
				break;
			case 'l':
				momentum_value += momentum_symbol.factor * l;
				break;
			case 'q':
				momentum_value += momentum_symbol.factor * q;
				break;
			default:
				throw std::runtime_error("Momentum not recognized! " + momentum_symbol.name);
			}
		}
		return momentum_value;
	}

	c_complex ModeHelper::get_expectation_value(mrock::symbolic_operators::WickOperator const& op, c_float momentum) const
	{
		if (op.type == mrock::symbolic_operators::Number_Type) {
			return model->occupation(momentum);
		}
		else if (op.type == mrock::symbolic_operators::SC_Type) {
			if (op.is_daggered) {
				return __conj(model->sc_expectation_value(momentum));
			}
			return model->sc_expectation_value(momentum);
		}
		throw std::runtime_error("Expectation value not recognized!");
	}

	void ModeHelper::createStartingStates()
	{
#ifndef _XP
		starting_states.resize(2, _parent::Vector::Zero(total_matrix_size));
		std::fill(starting_states[0].begin(), starting_states[0].begin() + m_iterator::max_idx(), sqrt(model->momentumRanges.INNER_STEP));
		std::fill(starting_states[1].begin() + 3 * m_iterator::max_idx(), starting_states[1].end(), sqrt(model->momentumRanges.INNER_STEP));
#else
		starting_states.push_back({ _parent::Vector::Zero(antihermitian_discretization), _parent::Vector::Zero(hermitian_discretization), "SC" });
		for (m_iterator it(&model->momentumRanges); it < m_iterator::max_idx(); ++it) {
			starting_states[0][0](it.idx) = 1. / sqrt((c_float)DISCRETIZATION);
			starting_states[0][1](it.idx) = 1. / sqrt((c_float)DISCRETIZATION);
		}
#endif
	}

	void ModeHelper::fillMatrices()
	{
#ifndef _XP
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
#ifndef _XP
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

		/* for (m_iterator it(&model->momentumRanges); it < m_iterator::max_idx(); ++it) {
			std::cout << it.k / model->fermi_wavevector << " -> " << K_plus(it.idx, m_iterator::max_idx() + it.idx) << "  " << K_plus(m_iterator::max_idx() + it.idx, it.idx) << " || "
				<< - 2 * model->interpolate_delta(it.k) * model->dispersion_to_fermi_level(it.k) / model->energy(it.k)  << std::endl;
		} */
#ifdef _XP
		std::cout << "||K_+ - K_+^+|| = " << (K_plus - K_plus.adjoint()).norm() << std::endl;
		std::cout << "||K_- - K_-^+|| = " << (K_minus - K_minus.adjoint()).norm() << std::endl;

		/* const auto ROD = residual_offdiagonality();
		std::cout << "ROD(+) = " << ROD.first << "    ROD(-) = " << ROD.second << std::endl;
		std::cout << model->dispersion_to_fermi_level(model->momentumRanges.index_to_momentum(_OUTER_DISC)) << "  "
			<< model->dispersion_to_fermi_level(model->momentumRanges.index_to_momentum(_OUTER_DISC + 1)) << std::endl; */

			//std::cout << K_plus.norm() << "   " << K_minus.norm() << std::endl;
#else
		std::cout << "||M - M^+|| = " << (M - M.adjoint()).norm() << std::endl;
		std::cout << "||N - N^+|| = " << (N - N.adjoint()).norm() << std::endl;
#endif
	}

	void ModeHelper::fill_M()
	{
#ifndef _XP
		M.setZero(total_matrix_size, total_matrix_size);
#else
		K_plus.setZero(hermitian_discretization, hermitian_discretization);
		K_minus.setZero(antihermitian_discretization, antihermitian_discretization);
#endif

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
#ifdef _XP
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size))
#endif
					fill_block_M(i, j);
			}
		}
	}

	void ModeHelper::fill_block_M(int i, int j)
	{
		for (m_iterator it(&model->momentumRanges); it < m_iterator::max_idx(); ++it) {
			if (is_zero(it.k)) continue;

			c_complex diag_buffer{};
			for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
				if (!term.delta_momenta.empty()) {
					if (term.sums.momenta.empty()) {
						if (term.coefficients.front().name == "g" || term.coefficients.front().name == "V") {
							// These kinds of terms scale as 1/V -> 0
							continue;
						}
					}
					diag_buffer += computeTerm(term, it.k, it.k);
					//M(i * m_iterator::max_idx() + it.idx, j * m_iterator::max_idx() + it.idx) += ieom_diag(it.k) * computeTerm(term, it.k, it.k);
				}
				else {
					for (m_iterator jt(&model->momentumRanges); jt < m_iterator::max_idx(); ++jt) {
						if (is_zero(jt.k)) continue;
						M(i * m_iterator::max_idx() + it.idx, j * m_iterator::max_idx() + jt.idx) += ieom_offdiag(it.k, jt.k) * computeTerm(term, it.k, jt.k);
					}
				}
			}
			M(i * m_iterator::max_idx() + it.idx, j * m_iterator::max_idx() + it.idx) += ieom_diag(it.k) * diag_buffer;
		}
	}

	void ModeHelper::fill_block_N(int i, int j)
	{
		for (m_iterator it(&model->momentumRanges); it < m_iterator::max_idx(); ++it) {
			if (is_zero(it.k)) continue;
			for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
				if (!term.delta_momenta.empty()) {
					// only k=l and k=-l should occur. Additionally, only the magntitude should matter
					N(i * m_iterator::max_idx() + it.idx, j * m_iterator::max_idx() + it.idx) += computeTerm(term, it.k, it.k);
				}
				else {
					throw std::runtime_error("Offdiagonal term in N!");
					/* for (m_iterator jt(&model->momentumRanges); jt < m_iterator::max_idx(); ++jt) {
						N(i * m_iterator::max_idx() + it.idx, j * m_iterator::max_idx() + jt.idx) += ieom_offdiag(it.k, jt.k) * computeTerm(term, it.k, jt.k);
						if(is_zero(jt.k)) continue;
					} */
				}
			}
			N(i * m_iterator::max_idx() + it.idx, j * m_iterator::max_idx() + it.idx) *= ieom_diag(it.k);
		}
	}

	c_complex ModeHelper::compute_phonon_sum(const mrock::symbolic_operators::WickTerm& term, c_float k, c_float l) const
	{
		const int q_dependend = term.which_operator_depends_on('q');
		mrock::symbolic_operators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		mrock::symbolic_operators::WickOperator const* const other_op = term.is_bilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		c_complex value{};
		if (summed_op->type == mrock::symbolic_operators::Number_Type) {
			value = model->phononInteraction.sc_channel_integral(model->occupation, k);
		}
		else {
			value = model->phononInteraction.sc_channel_integral(model->sc_expectation_value, k);
#ifdef _complex
			if (summed_op.is_daggered) value = std::conj(value);
#endif
		}
		if (other_op) {
			value *= this->get_expectation_value(*other_op, k);
		}
		value *= static_cast<c_float>(term.multiplicity);
		return value;
	}

	c_complex ModeHelper::compute_em_sum(const mrock::symbolic_operators::WickTerm& term, c_float k, c_float l) const
	{
		const int q_dependend = term.which_operator_depends_on('q');
		mrock::symbolic_operators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		mrock::symbolic_operators::WickOperator const* const other_op = term.is_bilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		c_complex value{};
		if (summed_op->type == mrock::symbolic_operators::Number_Type) {
			value = -(model->fock_coulomb(k) + model->interpolate_delta_n(k));
		}
		else {
			value = model->integral_screening(model->sc_expectation_value, k);
#ifdef _complex
			if (summed_op.is_daggered) value = std::conj(value);
#endif
		}
		if (other_op) {
			value *= this->get_expectation_value(*other_op, k);
		}
		value *= static_cast<c_float>(term.multiplicity);
		return value;
	}

	c_complex ModeHelper::computeTerm(const mrock::symbolic_operators::WickTerm& term, c_float k, c_float l) const
	{
		if (term.sums.momenta.empty()) {
			c_complex value{ static_cast<c_float>(term.multiplicity) };

			if (!term.coefficients.empty()) {
				const mrock::symbolic_operators::Coefficient* coeff_ptr = &term.coefficients.front();
				if (coeff_ptr->momenta.size() == 2U) {
					value *= model->compute_coefficient(*coeff_ptr, compute_momentum(coeff_ptr->momenta[0], k, l), compute_momentum(coeff_ptr->momenta[1], k, l));
				}
				else if (coeff_ptr->momenta.size() == 1U || coeff_ptr->momenta.empty()) {
					value *= model->compute_coefficient(*coeff_ptr, k, l);
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

	std::vector<c_float> ModeHelper::continuum_boundaries() const
	{
		const m_iterator buf(&model->momentumRanges);
		c_float prev_min { DBL_MAX };
		c_float min { DBL_MAX };
		for (InnerIterator it(&model->momentumRanges); it < InnerIterator::max_idx(); ++it) {
			if (this->model->energy(it.k) < prev_min)
			{
				prev_min = this->model->energy(it.k);
				min = it.k;
			}
		}

		std::cout << "Found min at " << min / model->fermi_wavevector 
			<< "k_F   E(min)=" << model->energy(min) 
			<< "   E(k_F)=" << model->energy(model->fermi_wavevector) << std::endl;
		return { 2 * model->energy(min),
			2 * std::max(model->energy(buf.max_k()), model->energy(buf.min_k())) };
	}

	int ModeHelper::hermitian_discretization = 0;
	int ModeHelper::antihermitian_discretization = 0;
	int ModeHelper::total_matrix_size = 0;

	ModeHelper::ModeHelper(ModelInitializer const& init)
		: _parent(this, SQRT_PRECISION,
#ifdef _XP
			m_iterator::max_idx()* hermitian_size, m_iterator::max_idx()* antihermitian_size, false,
#endif
			false)
	{
		hermitian_discretization = m_iterator::max_idx() * hermitian_size;
		antihermitian_discretization = m_iterator::max_idx() * antihermitian_size;
		total_matrix_size = m_iterator::max_idx() * number_of_basis_terms;

		model = std::make_unique<SCModel>(init);
		std::cout << "Working on " << model->info() << std::endl;
		wicks.load("../commutators/continuum/", true, number_of_basis_terms, 0);

#ifdef _iterative_selfconsistency
		auto solver = mrock::utility::Selfconsistency::make_iterative<c_complex>(model.get(), &model->Delta);
#else
		auto solver = mrock::utility::Selfconsistency::make_broyden<c_complex>(model.get(), &model->Delta, 200);
#endif
		solver.compute(true);

		const m_iterator buf(&model->momentumRanges);
		std::cout << "LB: <f_k> = " << model->sc_expectation_value(buf.min_k()) << "  k_min = "
			<< buf.min_k() / model->fermi_wavevector << "  E_min = " << model->energy(buf.min_k())
			<< "\nUB: <f_k> = " << model->sc_expectation_value(buf.max_k()) << "  k_max = " << buf.max_k() / model->fermi_wavevector
			<< "  E_max = " << model->energy(buf.max_k()) << std::endl;
	}
}
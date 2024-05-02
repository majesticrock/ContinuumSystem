#include "ModeHelper.hpp"

namespace Continuum {
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
            for(int k = 0; k < DISCRETIZATION; ++k){
                if(! term.delta_momenta.empty()){
                    // only k=l and k=-l should occur. Additionally, only the magntiude should matter
                    M(i + k * number_of_basis_terms, j + k * number_of_basis_terms) += computeTerm(term, k, k);
                } else {
                    for(int l = 0; l < DISCRETIZATION; ++l){
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
            for(int k = 0; k < DISCRETIZATION; ++k){
                if(! term.delta_momenta.empty()){
                    // only k=l and k=-l should occur. Additionally, only the magntiude should matter
                    M(i + k * number_of_basis_terms, j + k * number_of_basis_terms) += computeTerm(term, k, k);
                } else {
                    for(int l = 0; l < DISCRETIZATION; ++l){
                        M(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
                            += computeTerm(term, k, l).real();
                    }
                }
            }
        }
    }

    c_complex ModeHelper::computeTerm(const SymbolicOperators::WickTerm&, int k, int l)
    {
        if(term.isIdentity()){
            if(term.hasSingleCoefficient()){
                return term.getFactor() * this->model->computeCoefficient();
            }
        }
    }

    ModeHelper::ModeHelper(Utility::InputFileReader& input)
        : _parent(this, SQRT_PRECISION<c_float>), 
        number_of_basis_terms{ input.getInt("number_of_basis_terms") }, start_basis_at{ input.getInt("start_basis_at") }
    {
        model = std::make_uniqe<SCModel>(ModelInitializer(input));
        TOTAL_BASIS = DISCRETIZATION * this->number_of_basis_terms;
        wicks.load("../commutators/continuum", true, this->number_of_basis_terms, int start_basis_at);
    }
}
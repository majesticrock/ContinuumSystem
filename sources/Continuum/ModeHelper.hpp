#pragma once
#include "GlobalDefinitions.hpp"
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include "../../../Utility/sources/better_to_string.hpp"

namespace Continuum {
	class MatrixIsNegativeException : public std::runtime_error {
	public:
		c_float negative_eigenvalue{};
		MatrixIsNegativeException(const c_float& _negative_eigenvalue)
			: std::runtime_error("The matrix M is negative! Most negative eigenvalue = " + Utility::better_to_string(_negative_eigenvalue, std::chars_format::scientific, 6)),
			negative_eigenvalue(_negative_eigenvalue)
		{};
	};

	template <class EigenMatrixType>
	inline EigenMatrixType removeNoise(EigenMatrixType const& matrix) {
		return matrix.unaryExpr([](typename EigenMatrixType::Scalar const& val) {
			return (abs(val) < PRECISION<Utility::UnderlyingFloatingPoint_t<EigenMatrixType::Scalar>> ? typename EigenMatrixType::Scalar{} : val);
			});
	};

	class ModeHelper
	{
	protected:
		std::vector<SymbolicOperators::WickTermCollector> wicks_M{}, wicks_N{};

		size_t TOTAL_BASIS{};
	};
}
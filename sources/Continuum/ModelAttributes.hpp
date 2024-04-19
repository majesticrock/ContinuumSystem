#pragma once
#include "../../../Utility/sources/IsComplex.hpp"
#include "GlobalDefinitions.hpp"
#include <vector>

namespace Continuum {
	enum ComplexAttributePolicy { Magnitude, SeperateRealAndImaginary };

	template <typename DataType>
	struct ModelAttributes {
	public:
		/* On a 2D System we choose the following:
		* 0 - Delta_CDW
		* 1 - Delta_AFM
		* 2 - Delta_SC
		* 3 - Gamma_SC
		* 4 - Xi_SC (1D: Tau_SC, 3D: None? p-wave?) | For DOS: trying PS
		* 5 - Delta_eta
		* 6 - Gamma_n_up
		* 7 - Gamma_n_down
		* 8 - Delta_PS
		*/
		std::vector<DataType> selfconsistency_values;
		bool converged{};

		~ModelAttributes() = default;
		ModelAttributes() = default;
		ModelAttributes(std::initializer_list<DataType> i_list) : selfconsistency_values(i_list) {};
		ModelAttributes(ModelAttributes&& other) = default;
		ModelAttributes(const ModelAttributes& other) = default;
		ModelAttributes& operator=(const ModelAttributes& other) = default;
		ModelAttributes& operator=(ModelAttributes&& other) = default;

		// Using this constructor constructs the attribute vector with a fixed value, default is 0
		explicit ModelAttributes(const size_t number_of_attributes, const DataType& default_value = DataType{})
			: selfconsistency_values(number_of_attributes, default_value) {};

		ModelAttributes(const ModelAttributes<std::complex<DataType>>& other, ComplexAttributePolicy complexAttributePolicy)
			: selfconsistency_values(complexAttributePolicy == Magnitude ? other.selfconsistency_values.size() : 2U * other.selfconsistency_values.size()),
			converged{ other.converged }
		{
			if (complexAttributePolicy == Magnitude) {
				for (size_t i = 0U; i < selfconsistency_values.size(); ++i)
				{
					selfconsistency_values[i] = abs(other.selfconsistency_values[i]);
				}
			}
			else if (complexAttributePolicy == SeperateRealAndImaginary) {
				for (size_t i = 0; i < other.selfconsistency_values.size(); i++)
				{
					selfconsistency_values[i] = real(other.selfconsistency_values[i]);
					selfconsistency_values[i + other.selfconsistency_values.size()] = imag(other.selfconsistency_values[i]);
				}
			}
			else {
				throw std::runtime_error("ComplexAttributePolicy not recognized!");
			}
		};

		/*
		* Utility functions
		*/

		inline DataType& operator[](size_t i) {
			assert(i < selfconsistency_values.size());
			return selfconsistency_values[i];
		};
		inline const DataType& operator[](size_t i) const {
			assert(i < selfconsistency_values.size());
			return selfconsistency_values[i];
		};
		inline size_t size() const noexcept {
			return selfconsistency_values.size();
		}
		inline void push_back(const DataType& value) {
			selfconsistency_values.push_back(value);
		};
		inline void push_back(DataType&& value) {
			selfconsistency_values.push_back(std::move(value));
		};
		inline auto begin() {
			return selfconsistency_values.begin();
		}
		inline auto begin() const {
			return selfconsistency_values.begin();
		}
		inline auto end() {
			return selfconsistency_values.end();
		}
		inline auto end() const {
			return selfconsistency_values.end();
		}
		inline auto back() const {
			return selfconsistency_values.back();
		}
		template<class Vector>
		inline void fill_with(const Vector& vector) {
			this->selfconsistency_values.resize(vector.size());
			std::copy(vector.begin(), vector.end(), this->selfconsistency_values.begin());
		}
		inline void reset() {
			converged = false;
			std::fill(begin(), end(), DataType{});
		};
		inline void setZero() {
			this->reset();
		};

		inline bool isOrdered() const {
			for (const auto& value : selfconsistency_values)
			{
				if (!is_zero(value)) {
					return true;
				}
			}
			return false;
		};

		inline bool isFinite(const size_t i) const {
			return !is_zero(value);
		}

		/*
		* Arithmetric operators
		*/

		inline ModelAttributes& operator+=(const ModelAttributes& rhs) {
			for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
			{
				this->selfconsistency_values[i] += rhs.selfconsistency_values[i];
			}
			return *this;
		};
		inline ModelAttributes& operator-=(const ModelAttributes& rhs) {
			for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
			{
				this->selfconsistency_values[i] -= rhs.selfconsistency_values[i];
			}
			return *this;
		};
		template <class RealType>
		inline ModelAttributes& operator*=(const RealType rhs) {
			for (auto& value : this->selfconsistency_values)
			{
				value *= rhs;
			}
			return *this;
		};
		template <class RealType>
		inline ModelAttributes& operator/=(const RealType rhs) {
			for (auto& value : this->selfconsistency_values)
			{
				value /= rhs;
			}
			return *this;
		};
	};

	template <typename DataType>
	inline ModelAttributes<DataType> operator+(ModelAttributes<DataType> lhs, const ModelAttributes<DataType>& rhs) {
		return lhs += rhs;
	};
	template <typename DataType>
	inline ModelAttributes<DataType> operator-(ModelAttributes<DataType> lhs, const ModelAttributes<DataType>& rhs) {
		return lhs -= rhs;
	};

	template <typename DataType, class RealType>
	inline ModelAttributes<DataType> operator*(ModelAttributes<DataType> lhs, RealType rhs) {
		return lhs *= rhs;
	};

	template <typename DataType, class RealType>
	inline ModelAttributes<DataType> operator*(RealType lhs, ModelAttributes<DataType> rhs) {
		return rhs *= lhs;
	};

	template <typename DataType, class RealType>
	inline ModelAttributes<DataType> operator/(ModelAttributes<DataType> lhs, RealType rhs) {
		return lhs /= rhs;
	};

	template <typename DataType>
	inline std::ostream& operator<<(std::ostream& os, ModelAttributes<DataType> const& attributes) {
		for (const auto& value : attributes) {
			os << value << "  ";
		}
		return os;
	};
}
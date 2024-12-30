#pragma once
#include <mrock/utility/IsComplex.hpp>
#include "GlobalDefinitions.hpp"
#include <vector>
#include <complex>
#include <random>

namespace Continuum {
	enum ComplexAttributePolicy { Magnitude, SeperateRealAndImaginary };

	template <typename DataType>
	struct ModelAttributes {
	public:
		std::vector<DataType> selfconsistency_values;
		bool converged{};

		using value_type = DataType;
		using RealType = mrock::utility::UnderlyingFloatingPoint_t<DataType>;

		~ModelAttributes() = default;
		ModelAttributes() = default;
		ModelAttributes(std::initializer_list<DataType> i_list) : selfconsistency_values(i_list) {};
		ModelAttributes(ModelAttributes&& other) = default;
		ModelAttributes(const ModelAttributes& other) = default;
		ModelAttributes& operator=(const ModelAttributes& other) = default;
		ModelAttributes& operator=(ModelAttributes&& other) = default;

		inline static ModelAttributes Random(size_t size)
		{
			ModelAttributes<DataType> ret;
			ret.selfconsistency_values.resize(size);
			std::random_device dev;
			std::mt19937 rng(dev());
			std::uniform_real_distribution<> dis(0.0, 2.0 * 3.1415926);
			for (auto& value : ret) {
				if constexpr (mrock::utility::is_complex<DataType>()) {
					value = std::polar(dis(rng), dis(rng));
				}
				else {
					value = dis(rng);
				}
			}
			return ret;
		}

		inline static ModelAttributes Gaussian(int size, int center, int FWHM, double maximum) {
			ModelAttributes<DataType> ret;
			ret.selfconsistency_values.resize(size);
			for (int i = 0; i < size; ++i) {
				ret[i] = maximum * std::exp(-16. * 0.693147 * (i - center) * (i - center)
					/ static_cast<RealType>(FWHM * FWHM));
			}
			return ret;
		}

		template<class Allocator>
		inline static ModelAttributes FromAllocator(Allocator const& alloc, size_t size)
		{
			ModelAttributes<DataType> ret;
			ret.selfconsistency_values.resize(size);
			for (size_t i = 0U; i < size; ++i) {
				ret.selfconsistency_values[i] = alloc(i);
			}
			return ret;
		}

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
					selfconsistency_values[i] = std::abs(other.selfconsistency_values[i]);
				}
			}
			else if (complexAttributePolicy == SeperateRealAndImaginary) {
				for (size_t i = 0U; i < other.selfconsistency_values.size(); ++i)
				{
					selfconsistency_values[i] = std::real(other.selfconsistency_values[i]);
					selfconsistency_values[i + other.selfconsistency_values.size()] = std::imag(other.selfconsistency_values[i]);
				}
			}
			else {
				throw std::runtime_error("ComplexAttributePolicy not recognized!");
			}
		};

		/*
		* utility functions
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
		inline auto front() const {
			return selfconsistency_values.front();
		}
		template<class Vector>
		inline void fill_with(const Vector& vector) {
			this->selfconsistency_values.resize(vector.size());
			std::copy(vector.begin(), vector.end(), this->selfconsistency_values.begin());
		}
		template<class Vector>
		inline void fill_with(const Vector& vector, RealType weight) {
			assert(this->size() == vector.size());
			for (size_t i = 0U; i < this->size(); ++i) {
				this->selfconsistency_values[i] = (1. - weight) * this->selfconsistency_values[i] + weight * vector[i];
			}
		}
		inline void clear_noise(RealType precision) {
			for (auto& val : this->selfconsistency_values) {
				if (std::abs(val) < precision) val = DataType{};
			}
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
			return !is_zero(this->selfconsistency_values[i]);
		}

		inline ModelAttributes<RealType> real() const {
			if constexpr (mrock::utility::is_complex<DataType>()) {
				ModelAttributes<RealType> ret;
				ret.converged = this->converged;
				ret.selfconsistency_values.resize(this->size());
				for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
				{
					ret.selfconsistency_values[i] = std::real(this->selfconsistency_values[i]);
				}
				return ret;
			}
			else {
				return *this;
			}
		};
		inline ModelAttributes<RealType> imag() const {
			if constexpr (mrock::utility::is_complex<DataType>()) {
				ModelAttributes<RealType> ret;
				ret.converged = this->converged;
				ret.selfconsistency_values.resize(this->size());
				for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
				{
					ret.selfconsistency_values[i] = std::imag(this->selfconsistency_values[i]);
				}
				return ret;
			}
			else {
				ModelAttributes<DataType> ret;
				ret.converged = this->converged;
				ret.selfconsistency_values.resize(this->size());
				return ret;
			}
		};
		inline ModelAttributes<RealType> abs() const {
			if constexpr (mrock::utility::is_complex<DataType>()) {
				return ModelAttributes<RealType>(*this, Magnitude);
			}
			else {
				ModelAttributes<DataType> ret(*this);
				for (auto& val : ret) {
					val = std::abs(val);
				}
				return ret;
			}
		};
		inline const std::vector<DataType>& as_vector() const {
			return this->selfconsistency_values;
		};

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
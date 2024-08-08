#pragma once
#include "ModelInitializer.hpp"

namespace Continuum {
	struct Base_Incrementer {
		const c_float _Delta;
		Base_Incrementer(c_float Delta) : _Delta{ Delta } {}
		virtual ~Base_Incrementer() = default;
		virtual void increment(ModelInitializer& init, int n = 1) const = 0;
		virtual double current(ModelInitializer const& init) const = 0;
	};
	struct Temperature_Incrementer : public Base_Incrementer {
		Temperature_Incrementer(c_float Delta) : Base_Incrementer(Delta) {}
		void increment(ModelInitializer& init, int n = 1) const override;
		double current(ModelInitializer const& init) const override;
	};
	struct PhononCoupling_Incrementer : public Base_Incrementer {
		PhononCoupling_Incrementer(c_float Delta) : Base_Incrementer(Delta) {}
		void increment(ModelInitializer& init, int n = 1) const override;
		double current(ModelInitializer const& init) const override;
	};
	struct DebyeFrequency_Incrementer : public Base_Incrementer {
		DebyeFrequency_Incrementer(c_float Delta) : Base_Incrementer(Delta) {}
		void increment(ModelInitializer& init, int n = 1) const override;
		double current(ModelInitializer const& init) const override;
	};
	struct FermiEnergy_Incrementer : public Base_Incrementer {
		FermiEnergy_Incrementer(c_float Delta) : Base_Incrementer(Delta) {}
		void increment(ModelInitializer& init, int n = 1) const override;
		double current(ModelInitializer const& init) const override;
	};
	struct CoulombScaling_Incrementer : public Base_Incrementer {
		CoulombScaling_Incrementer(c_float Delta) : Base_Incrementer(Delta) {}
		void increment(ModelInitializer& init, int n = 1) const override;
		double current(ModelInitializer const& init) const override;
	};
}
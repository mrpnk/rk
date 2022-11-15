#pragma once
#include <array>
#include <cmath>
#include <iostream>

template<int N>
struct vec : std::array<double, N>{
	vec(){}
	vec(std::initializer_list<double> parts) {
		for (int i = 0; auto & a:parts)
			(*this)[i++] = a;
	}

	template<int M1, int M2>
	vec(double v0, vec<M1> const& m1, vec<M2> const& m2) {
		static_assert(1+M1+M2 == N);

		int i = 0;
		(*this)[i++] = v0;
		for (; i < 1 + M1;)
			(*this)[i++] = m1[i - 1];
		for (; i < 1 + M1 + M2;)
			(*this)[i++] = m2[i - 1 - M1];
	}

	vec operator*(double f) const {
		vec ret = *this;
		for(auto& a : ret) a *= f;
		return ret;
	}
	vec operator+(double f) const {
		vec ret = *this;
		for(auto& a : ret) a += f;
		return ret;
	}
	vec operator+(vec const& v) const {
		vec ret;
		for(int i = 0; i < N; ++i) ret[i] = (*this)[i] + v[i];
		return ret;
	}
	vec operator-(vec const& v) const {
		vec ret;
		for(int i = 0; i < N; ++i) ret[i] = (*this)[i] - v[i];
		return ret;
	}
	vec operator/(vec const& v) const {
		vec ret;
		for(int i=0; i<N;++i) ret[i] = (*this)[i]/v[i];
		return ret;
	}
	bool isbad() const {
		for(int i = 0; i < N; ++i) if(std::isnan((*this)[i]))return true;
		return false;
	}
	friend vec<N> abs(vec<N> const& v) {
		vec ret;
		for(int i = 0; i < N; ++i) ret[i] = abs(v[i]);
		return ret;
	}
	friend double reduceMax(vec<N> const& v) {
		double ret=-std::numeric_limits<double>::max();
		for(int i = 0; i < N; ++i) ret = std::max(ret,v[i]);
		return ret;
	}
	friend std::ostream& operator<<(std::ostream& os, vec const& v) {
		for (int i = 0; i < N; ++i)os << v[i] << '\t';
		return os;
	}
};


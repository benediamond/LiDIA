// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Frank Lehmann (FL), Markus Maurer (MM)
//                Thorsten Rottschaefer (TR), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================

#ifndef LIDIA_DENSE_POWER_SERIES_GF_ELEMENT_H_GUARD_
#define LIDIA_DENSE_POWER_SERIES_GF_ELEMENT_H_GUARD_

#ifndef LIDIA_GF_ELEMENT_H_GUARD_
#include "LiDIA/gf_element.h"
#endif
#ifndef LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_
#include "LiDIA/finite_fields/base_dense_power_series.h"
#endif
#ifndef LIDIA_DENSE_POWER_SERIES_H_GUARD_
#include "LiDIA/dense_power_series.h"
#endif
#ifndef LIDIA_SPARSE_POWER_SERIES_H_GUARD_
#include "LiDIA/sparse_power_series.h"
#endif
#ifndef LIDIA_FFT_PRIME_H_GUARD_ // ?
#include "LiDIA/fft_prime.h"
#endif

#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#define IN_NAMESPACE_LIDIA
#endif

#ifndef NO_PSR_GF_ELEMENT

//******************************************************************************
//************* Specialization : dense_power_series < gf_element > *************
//******************************************************************************

template <>
class dense_power_series<gf_element>
		: public base_dense_power_series<gf_element>
{

protected:
	//
	// ***** protected member functions *****
	//

	void square(const dense_power_series<gf_element> &a);
	void multiply(const dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);
	void invert(const dense_power_series<gf_element> &a);
	void power(const dense_power_series<gf_element> &a, long n);
	void divide(const dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);
	void divide(const gf_element &b, const dense_power_series<gf_element> &a);

public:
	//
	// ***** constructors / destructor *****
	//

	dense_power_series();
	dense_power_series(const gf_element &a, lidia_size_t l);
	dense_power_series(const base_vector<gf_element> &a, lidia_size_t f);
	dense_power_series(const base_dense_power_series<gf_element> &a);
	~dense_power_series() {}

	//
	// ***** assignment - operator *****
	//

	dense_power_series<gf_element> &operator=(
			const dense_power_series<gf_element> &a);
	dense_power_series<gf_element> &operator=(
			const sparse_power_series<gf_element> &a);

	//
	// ***** arithmetic via functions *****
	//

	friend void square(dense_power_series<gf_element> &c,
			const dense_power_series<gf_element> &a);

	friend void multiply(dense_power_series<gf_element> &c,
			const dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);

	friend void invert(dense_power_series<gf_element> &c,
			const dense_power_series<gf_element> &a);

	friend void power(dense_power_series<gf_element> &c,
			const dense_power_series<gf_element> &a, long n);

	friend void divide(dense_power_series<gf_element> &c,
			dense_power_series<gf_element> &a, dense_power_series<gf_element> &b);

	friend void divide(dense_power_series<gf_element> &c, const gf_element &b,
			const dense_power_series<gf_element> &a);

	//
	// ***** arithmetic via operators *****
	//

	friend dense_power_series<gf_element> operator*(
			const dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);

	friend dense_power_series<gf_element> &operator*=(
			dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);

	friend dense_power_series<gf_element> operator/(
			const dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);

	friend dense_power_series<gf_element> operator/(
			const gf_element &b, const dense_power_series<gf_element> &a);

	friend dense_power_series<gf_element> &operator/=(
			dense_power_series<gf_element> &a,
			const dense_power_series<gf_element> &b);
};

inline void square(
		dense_power_series<gf_element> &c, const dense_power_series<gf_element> &a)
{
	c.square(a);
}

inline void multiply(dense_power_series<gf_element> &c,
		const dense_power_series<gf_element> &a,
		const dense_power_series<gf_element> &b)
{
	c.multiply(a, b);
}

inline void invert(
		dense_power_series<gf_element> &c, const dense_power_series<gf_element> &a)
{
	c.invert(a);
}

inline void power(dense_power_series<gf_element> &c,
		const dense_power_series<gf_element> &a, long n)
{
	c.power(a, n);
}

inline void divide(dense_power_series<gf_element> &c,
		dense_power_series<gf_element> &a, dense_power_series<gf_element> &b)
{
	c.divide(a, b);
}

inline void divide(dense_power_series<gf_element> &c, const gf_element &b,
		const dense_power_series<gf_element> &a)
{
	c.divide(b, a);
}

//
// ***** arithmetic via operators *****
//

inline dense_power_series<gf_element> operator*(
		const dense_power_series<gf_element> &a,
		const dense_power_series<gf_element> &b)
{
	dense_power_series<gf_element> c;

	c.multiply(a, b);
	return c;
}

inline dense_power_series<gf_element> &operator*=(
		dense_power_series<gf_element> &a, const dense_power_series<gf_element> &b)
{
	a.multiply(a, b);
	return a;
}

inline dense_power_series<gf_element> operator/(
		const dense_power_series<gf_element> &a,
		const dense_power_series<gf_element> &b)
{
	dense_power_series<gf_element> c;

	c.divide(a, b);
	return c;
}

inline dense_power_series<gf_element> operator/(
		const gf_element &b, const dense_power_series<gf_element> &a)
{
	dense_power_series<gf_element> c;

	c.divide(b, a);
	return c;
}

inline dense_power_series<gf_element> &operator/=(
		dense_power_series<gf_element> &a, const dense_power_series<gf_element> &b)
{
	a.divide(a, b);
	return a;
}

#endif // NO_PSR_GF_ELEMENT

#ifdef LIDIA_NAMESPACE
} // end of namespace LiDIA
#undef IN_NAMESPACE_LIDIA
#endif

#endif // LIDIA_DENSE_POWER_SERIES_GF_ELEMENT_H_GUARD_

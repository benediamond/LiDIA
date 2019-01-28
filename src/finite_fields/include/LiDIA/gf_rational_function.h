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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================

//	Description: This class realizes rational functions and
//              rational functions modulo some modulus which can
//              be either given as Fp_polynomial or as Fp_poly_modulus.

#ifndef LIDIA_GF_RATIONAL_FUNCTION_H_GUARD_
#define LIDIA_GF_RATIONAL_FUNCTION_H_GUARD_

#ifndef LIDIA_LIDIA_H_GUARD_
#include "LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_GF_POLYNOMIAL_H_GUARD_
#include "LiDIA/gf_polynomial.h"
#endif

#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#define IN_NAMESPACE_LIDIA
#endif

class gf_rational_function
{
private:
	gf_polynomial *num;
	gf_polynomial *den;

public:
	gf_rational_function();
	gf_rational_function(const galois_field &K);
	gf_rational_function(const gf_polynomial &);
	gf_rational_function(const gf_polynomial &, const gf_polynomial &);
	gf_rational_function(const gf_rational_function &);
	~gf_rational_function();

	//******************************************************************
	// some simple assignment functions   etc.
	//******************************************************************

	// void kill();

	void set_field(const galois_field &K);
	const galois_field &get_field() const;

	lidia_size_t degree_numerator() const;
	lidia_size_t degree_denominator() const;

	void get_coefficient_numerator(gf_element &a, lidia_size_t i) const;
	void get_coefficient_denominator(gf_element &a, lidia_size_t i) const;

	void set_coefficient_numerator(const gf_element &a, lidia_size_t i);
	void set_coefficient_denominator(const gf_element &a, lidia_size_t i);

	void set_coefficient_numerator(lidia_size_t i);
	void set_coefficient_denominator(lidia_size_t i);

	const gf_element lead_coefficient_numerator() const;
	const gf_element lead_coefficient_denominator() const;

	const gf_element const_term_numerator() const;
	const gf_element const_term_denominator() const;

	//*************************************************************
	// Assignments
	//*************************************************************

	gf_rational_function &operator=(const gf_rational_function &);
	gf_rational_function &operator=(const gf_polynomial &);

	void assign(const gf_rational_function &f);

	void assign(const gf_polynomial &f);

	void assign(const gf_polynomial &f, const gf_polynomial &g);

	void assign_numerator(const gf_polynomial &a);

	void assign_denominator(const gf_polynomial &a);

	void assign_zero();
	void assign_one();
	void assign_x();

	// void randomize(lidia_size_t deg_num, lidia_size_t deg_denom = 0);

	gf_polynomial &numerator();
	gf_polynomial &denominator();

	const gf_polynomial &numerator() const;
	const gf_polynomial &denominator() const;

	bigint operator()(const bigint &a) const;

	//******* Comparisons ***************************************

	friend bool operator==(
			const gf_rational_function &, const gf_rational_function &);
	friend bool operator!=(
			const gf_rational_function &a, const gf_rational_function &b);

	friend bool equal_mod(const gf_rational_function &,
			const gf_rational_function &, const gf_polynomial &);

	friend bool equal(const gf_rational_function &, const gf_rational_function &,
			const gf_poly_modulus &);

	bool is_zero() const;
	bool is_one() const;

	//***************************************************************
	// procedural versions for arithmetics
	//***************************************************************

#ifndef HEADBANGER

	void negate() { num->negate(); }

	friend void negate(gf_rational_function &c, const gf_rational_function &a)
	{
		c.assign(a);
		c.negate();
	}

	friend void add(gf_rational_function &, const gf_rational_function &,
			const gf_rational_function &);

	friend void subtract(gf_rational_function &, const gf_rational_function &,
			const gf_rational_function &);

	friend void multiply(gf_rational_function &c, const gf_rational_function &a,
			const gf_rational_function &b)
	{
		multiply(*c.num, *a.num, *b.num);
		multiply(*c.den, *a.den, *b.den);
	}

	friend void multiply(gf_rational_function &x, const gf_rational_function &a,
			const gf_polynomial &b)
	{
		multiply(*x.num, *a.num, b);
		x.den->assign(*a.den);
	}

	friend void multiply(gf_rational_function &x, const gf_polynomial &a,
			const gf_rational_function &b)
	{
		multiply(*x.num, a, *b.num);
		x.den->assign(*b.den);
	}

	friend void multiply(gf_rational_function &x, const gf_element &a,
			const gf_rational_function &b)
	{
		multiply(*x.num, *b.num, a);
		x.den->assign(*b.den);
	}

	friend void multiply(gf_rational_function &c, const gf_rational_function &a,
			const gf_element &b)
	{
		multiply(*c.num, *a.num, b);
		c.den->assign(*a.den);
	}

	void multiply_by_2() { add(*num, *num, *num); }

	void normalize()
	{
		gf_polynomial c(gcd(*num, *den));
		divide(*num, *num, c);
		divide(*den, *den, c);
	}

	friend void square(gf_rational_function &x, const gf_rational_function &a)
	{
		square(*x.num, *a.num);
		square(*x.den, *a.den);
	}

#ifndef HEADBANGER
	friend void div_rem(gf_polynomial &q, gf_rational_function &f)
	{
		gf_polynomial rr;

		div_rem(q, rr, *f.num, *f.den);
		f.num->assign(rr);
	}

	friend void divide(gf_rational_function &q, const gf_rational_function &a,
			const gf_rational_function &b);

	friend void divide(gf_rational_function &q, const gf_polynomial &a,
			const gf_rational_function &b);

	friend void divide(gf_rational_function &q, const gf_rational_function &a,
			const gf_polynomial &b)
	{
		q.den->assign(*a.den);
		multiply(*q.num, *a.num, b);
	}
#endif

	void invert();

	friend void invert(gf_rational_function &c, const gf_rational_function &a);

	void reduce(const gf_polynomial &f)
	{
		remainder(*num, *num, f);
		remainder(*den, *den, f);
	}

	void reduce(const gf_poly_modulus &f)
	{
		remainder(*num, *num, f.modulus());
		remainder(*den, *den, f.modulus());
	}

	//***************************************************************
	// operators
	//***************************************************************

	gf_rational_function &operator+=(const gf_rational_function &a)
	{
		add(*this, *this, a);
		return *this;
	}

	gf_rational_function &operator-=(const gf_rational_function &a)
	{
		subtract(*this, *this, a);
		return *this;
	}

	gf_rational_function &operator*=(const gf_rational_function &a)
	{
		multiply(*this, *this, a);
		return *this;
	}

	gf_rational_function &operator/=(const gf_rational_function &a)
	{
		divide(*this, *this, a);
		return *this;
	}

	//***************************************************************
	// Miscellaneaous
	//***************************************************************

	friend void swap(gf_rational_function &x, gf_rational_function &y)
	{
		swap(*x.num, *y.num);
		swap(*x.den, *y.den);
	}

	// friend void shift(gf_rational_function &c, const gf_rational_function &a,
	//                   lidia_size_t n) // c = a * x^n

	friend void derivative(
			gf_rational_function &c, const gf_rational_function &a);
	// c = derivative of a

	friend gf_rational_function derivative(const gf_rational_function &a)
	{
		gf_rational_function x;

		derivative(x, a);
		return x;
	}

	//***************************************************************
	// operators
	//***************************************************************

	friend gf_rational_function operator-(const gf_rational_function &a)
	{
		gf_rational_function c(a);

		c.negate();
		return c;
	}

	friend gf_rational_function operator+(
			const gf_rational_function &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		add(c, a, b);
		return c;
	}

	friend gf_rational_function operator+(
			const gf_polynomial &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		add(c, gf_rational_function(a), b);
		return c;
	}

	friend gf_rational_function operator+(
			const gf_rational_function &a, const gf_polynomial &b)
	{
		gf_rational_function c(a.get_field());

		add(c, gf_rational_function(b), a);
		return c;
	}

	friend gf_rational_function operator-(
			const gf_rational_function &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		subtract(c, a, b);
		return c;
	}

	friend gf_rational_function operator-(
			const gf_polynomial &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		subtract(c, gf_rational_function(a), b);
		return c;
	}

	friend gf_rational_function operator-(
			const gf_rational_function &a, const gf_polynomial &b)
	{
		gf_rational_function c(a.get_field());

		subtract(c, a, gf_rational_function(b));
		return c;
	}

	friend gf_rational_function operator*(
			const gf_rational_function &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		multiply(c, b, a);
		return c;
	}

	friend gf_rational_function operator*(
			const gf_polynomial &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		multiply(c, b, a);
		return c;
	}

	friend gf_rational_function operator*(
			const gf_rational_function &a, const gf_polynomial &b)
	{
		gf_rational_function c(a.get_field());

		multiply(c, b, a);
		return c;
	}

	friend gf_rational_function operator*(
			const gf_element &a, const gf_rational_function &b)
	{
		gf_rational_function c(b.get_field());

		multiply(c, a, b);
		return c;
	}

	friend gf_rational_function operator/(
			const gf_rational_function &a, const gf_rational_function &b)
	{
		gf_rational_function c(a.get_field());

		divide(c, a, b);
		return c;
	}

	friend gf_rational_function operator/(
			const gf_polynomial &a, const gf_rational_function &b)
	{
		gf_rational_function c(b.get_field());

		divide(c, a, b);
		return c;
	}

	friend gf_rational_function operator/(
			const gf_rational_function &a, const gf_polynomial &b)
	{
		gf_rational_function c(a.get_field());

		divide(c, a, b);
		return c;
	}
	//***************************************************************
	//
	//	Modular Arithmetic without pre-conditioning
	//
	//***************************************************************

#ifndef HEADBANGER
	// arithmetic mod f.
	// all inputs and outputs are polynomials of degree less than deg(f).
	// Note that numerator and denominator are reduced modulo f, but
	// no conversion to gp_polynomial is done.
	// ASSUMPTION: f is assumed monic, and deg(f) > 0.
	// NOTE: if you want to do many computations with a fixed f,
	//       use the gp_poly_modulus data structure and associated
	//       routines below.

	friend void add_mod(gf_rational_function &c, const gf_rational_function &a,
			const gf_rational_function &b, const gf_polynomial &f);

	friend void subtract_mod(gf_rational_function &c,
			const gf_rational_function &a, const gf_rational_function &b,
			const gf_polynomial &f);

	friend void multiply_mod(gf_rational_function &c,
			const gf_rational_function &a, const gf_rational_function &b,
			const gf_polynomial &f)
	{
		multiply_mod(*c.num, *a.num, *b.num, f);
		multiply_mod(*c.den, *a.den, *b.den, f);
	}

	friend void square_mod(gf_rational_function &c, const gf_rational_function &a,
			const gf_polynomial &f)
	{
		square_mod(*c.num, *a.num, f);
		square_mod(*c.den, *a.den, f);
	}

	friend void convert_mod(
			gf_polynomial &c, const gf_rational_function &a, const gf_polynomial &f)
	// c = a % f, error if *a.den is not invertible
	{
		invert_mod(c, *a.den, f);
		multiply_mod(c, c, *a.num, f);
	}

	friend bool convert_mod_status(
			gf_polynomial &c, const gf_rational_function &a, const gf_polynomial &f)
	// if (a.den, f) = 1, returns 1 and sets c = a % f
	// otherwise, returns 0 and sets c = (a.den, f)
	{
		bool t = invert_mod_status(c, *a.den, f);
		if (t)
			multiply_mod(c, c, *a.num, f);
		return t;
	}

	friend void divide_mod(gf_rational_function &q, const gf_rational_function &a,
			const gf_rational_function &b, const gf_polynomial &f)
	{
		gf_rational_function h(b);

		h.invert();
		multiply_mod(q, a, h, f);
	}
#endif

	//***************************************************************
	// now the gp_poly_modulus functions for faster arithmetic
	//***************************************************************

	// If you need to do a lot of arithmetic modulo a fixed f,
	// build gp_poly_modulus F for f.  This pre-computes information about f
	// that speeds up the computation a great deal.
	// f should be monic, and deg(f) > 0.

	friend void add(gf_rational_function &x, const gf_rational_function &a,
			const gf_rational_function &b, const gf_poly_modulus &F);

	friend void subtract(gf_rational_function &x, const gf_rational_function &a,
			const gf_rational_function &b, const gf_poly_modulus &F);

	friend void multiply(gf_rational_function &x, const gf_rational_function &a,
			const gf_rational_function &b, const gf_poly_modulus &F)
	// x = (a * b) % f, deg(a), deg(b) < f.degree()
	{
		multiply(*x.num, *a.num, *b.num, F);
		multiply(*x.den, *a.den, *b.den, F);
	}

	friend void square(gf_rational_function &x, const gf_rational_function &a,
			const gf_poly_modulus &F)
	// x = a^2 % f, deg(a) < f.degree()
	{
		square(*x.num, *a.num, F);
		square(*x.den, *a.den, F);
	}

	friend void divide(gf_rational_function &x, const gf_rational_function &a,
			const gf_rational_function &b, const gf_poly_modulus &F)
	{
		gf_rational_function h(b);

		h.invert();
		multiply(x, a, h, F);
	}

#endif // HEADBANGER // what in god's name is headbanger?
};

#ifdef LIDIA_NAMESPACE
} // end of namespace LiDIA
#undef IN_NAMESPACE_LIDIA
#endif

#endif // LIDIA_FP_RATIONAL_FUNCTION_H_GUARD_

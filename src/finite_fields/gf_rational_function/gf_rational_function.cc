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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "LiDIA/gf_rational_function.h"
#include <cctype>

#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#endif

// ---------- CONSTRUCTORS / DESTRUCTORS ----------

gf_rational_function::gf_rational_function()
{
	num = new polynomial<gf_element>;
	den = new polynomial<gf_element>;
	set_field(galois_field(3)); // ? set_modulus(3);
	den->assign_one();
	num->assign_zero();
}

gf_rational_function::gf_rational_function(const galois_field &K)
{
	num = new gf_polynomial;
	den = new gf_polynomial;
	if (K.characteristic().is_even())
		lidia_error_handler("gf_rational_function", "ct::no odd prime number");
	set_field(K);
	den->assign_one();
	num->assign_zero();
}

gf_rational_function::gf_rational_function(const gf_polynomial &f)
{
	num = new gf_polynomial;
	den = new gf_polynomial;
	set_field(f.get_field());
	num->assign(f);
	den->assign_one();
}

gf_rational_function::gf_rational_function(
		const gf_polynomial &f, const gf_polynomial &g)
{
	if (g.get_field() != f.get_field()) // sketchy. special method for this?
		lidia_error_handler("gf_rational_function", "ct::different moduli");
	num = new gf_polynomial;
	den = new gf_polynomial;
	set_field(f.get_field());
	num->assign(f);
	den->assign(g);
}

gf_rational_function::gf_rational_function(const gf_rational_function &f)
{
	num = new gf_polynomial;
	den = new gf_polynomial;
	set_field(f.get_field());
	num->assign(*f.num);
	den->assign(*f.den);
}

gf_rational_function::~gf_rational_function()
{
	if (num != NULL)
		delete num;

	if (den != NULL)
		delete den;
}

//******************************************************************
// some simple assignment functions   etc.
//******************************************************************

void gf_rational_function::set_field(const galois_field &K)
{
	num->set_field(K);
	den->set_field(K);
}

const galois_field &gf_rational_function::get_field() const
{
	return num->get_field();
}

lidia_size_t gf_rational_function::degree_numerator() const
{
	return num->degree();
}

lidia_size_t gf_rational_function::degree_denominator() const
{
	if (num->is_zero())
		return -1;
	else
		return den->degree();
}

void gf_rational_function::get_coefficient_numerator(
		gf_element &a, lidia_size_t i) const
{
	num->get_coefficient(a, i);
}

void gf_rational_function::get_coefficient_denominator(
		gf_element &a, lidia_size_t i) const
{
	den->get_coefficient(a, i);
}

void gf_rational_function::set_coefficient_numerator(
		const gf_element &a, lidia_size_t i)
{
	num->set_coefficient(a, i);
}

void gf_rational_function::set_coefficient_denominator(
		const gf_element &a, lidia_size_t i)
{
	den->set_coefficient(a, i);
}

void gf_rational_function::set_coefficient_numerator(lidia_size_t i)
{
	num->set_coefficient(i);
}

void gf_rational_function::set_coefficient_denominator(lidia_size_t i)
{
	den->set_coefficient(i);
}

const gf_element gf_rational_function::lead_coefficient_numerator() const
{
	return num->lead_coeff();
}

const gf_element gf_rational_function::lead_coefficient_denominator() const
{
	if (num->is_zero())
	{
		lidia_error_handler("gf_rational_function",
				"lead_coeff_denominator::zero rational function");
		return num->get_field();
	}
	else
		return den->lead_coeff();
}

const gf_element gf_rational_function::const_term_numerator() const
{
	return num->const_term();
}

const gf_element gf_rational_function::const_term_denominator() const
{
	if (num->is_zero())
	{
		lidia_error_handler("gf_rational_function",
				"const_term_denominator::zero rational function");
		return num->get_field();
	}
	else
		return den->const_term();
}

//********** Assignments **************************************

gf_rational_function &gf_rational_function::operator=(
		const gf_rational_function &f)
{
	if (this != &f)
		this->assign(*f.num, *f.den);
	return *this;
}

gf_rational_function &gf_rational_function::operator=(const gf_polynomial &f)
{
	set_field(f.get_field());
	num->assign(f);
	den->assign_one();
	return *this;
}

void gf_rational_function::assign(const gf_rational_function &f)
{
	num->assign(*f.num);
	den->assign(*f.den);
}

void gf_rational_function::assign(
		const gf_polynomial &f, const gf_polynomial &g)
{
	if (f.get_field() != g.get_field())
		lidia_error_handler("gf_rational_function", "assign::different moduli");

	num->assign(f);
	den->assign(g);
}

void gf_rational_function::assign(const gf_polynomial &f)
{
	set_field(f.get_field());
	num->assign(f);
	den->assign_one();
}

void gf_rational_function::assign_numerator(const gf_polynomial &a)
{
	if (a.get_field() != num->get_field())
		lidia_error_handler(
				"gf_rational_function", "assign_numerator::different moduli");
	num->assign(a);
}

void gf_rational_function::assign_denominator(const gf_polynomial &a)
{
	if (a.get_field() != num->get_field())
		lidia_error_handler(
				"gf_rational_function", "assign_denominator::different moduli");
	den->assign(a);
}

void gf_rational_function::assign_zero()
{
	num->assign_zero();
	den->assign_one();
}

void gf_rational_function::assign_one()
{
	num->assign_one();
	den->assign_one();
}

void gf_rational_function::assign_x()
{
	num->assign_x();
	den->assign_one();
}

// void gf_rational_function::randomize(lidia_size_t deg_num,
//                                      lidia_size_t deg_denom) {

gf_polynomial &gf_rational_function::numerator()
{
	return static_cast<gf_polynomial &>(*num);
}

gf_polynomial &gf_rational_function::denominator()
{
	return static_cast<gf_polynomial &>(*den);
}

const gf_polynomial &gf_rational_function::numerator() const
{
	return static_cast<const gf_polynomial &>(*num);
}

const gf_polynomial &gf_rational_function::denominator() const
{
	return static_cast<const gf_polynomial &>(*den);
}

//************ Comparisons *************************************

bool operator==(const gf_rational_function &a, const gf_rational_function &b)
{
	gf_polynomial h1, h2;

	multiply(h1, *a.num, *b.den);
	multiply(h2, *a.den, *b.num);
	return (h1 == h2);
}

bool operator!=(const gf_rational_function &a, const gf_rational_function &b)
{
	return (!(a == b));
}

bool gf_rational_function::is_zero() const { return num->is_zero(); }

bool gf_rational_function::is_one() const { return ((*num) == (*den)); }

bool equal_mod(const gf_rational_function &a, const gf_rational_function &b,
		const gf_polynomial &f)
{
	gf_polynomial h1, h2;

	if ((b.is_zero() && !a.is_zero()) || (a.is_zero() && !b.is_zero()))
		return false;

	if (a.is_zero() && b.is_zero())
		return true;

	multiply_mod(h1, *a.num, *b.den, f);
	multiply_mod(h2, *a.den, *b.num, f);
	return (h1 == h2);
}

bool equal(const gf_rational_function &a, const gf_rational_function &b,
		const gf_poly_modulus &f)
{
	gf_polynomial h1, h2;

	if (a.is_zero() != b.is_zero())
		return false;

	if (a.is_zero() && b.is_zero())
		return true;

	multiply(h1, *a.num, *b.den, f);
	multiply(h2, *a.den, *b.num, f);
	return (h1 == h2);
}

// gf_element gf_rational_function::operator()(const gf_element &a) const {

//******** procedural versions for arithmetic *****************

void add(gf_rational_function &x, const gf_rational_function &a,
		const gf_rational_function &b)
{
	if (&a == &b)
	{
		add(*x.num, *a.num, *a.num);
		x.den->assign(*a.den);
	}
	else if ((*a.den) == (*b.den))
	{
		x.den->assign(*a.den);
		add(*x.num, *a.num, *b.num);
	}
	else
	{
		gf_polynomial h1, h2;

		multiply(h1, *a.num, *b.den);
		multiply(h2, *a.den, *b.num);

		multiply(*x.den, *a.den, *b.den);
		add(*x.num, h1, h2);
	}
}

void add_mod(gf_rational_function &x, const gf_rational_function &a,
		const gf_rational_function &b, const gf_polynomial &f)
{
	if (&a == &b)
	{
		add(*x.num, *a.num, *a.num);
		x.den->assign(*a.den);
	}
	else if ((*a.den) == (*b.den))
	{
		x.den->assign(*a.den);
		add(*x.num, *a.num, *b.num);
	}
	else
	{
		gf_polynomial h1, h2;
		multiply_mod(h1, *a.num, *b.den, f);
		multiply_mod(h2, *a.den, *b.num, f);
		multiply_mod(*x.den, *a.den, *b.den, f);
		add(*x.num, h1, h2);
	}
}

void add(gf_rational_function &x, const gf_rational_function &a,
		const gf_rational_function &b, const gf_poly_modulus &F)
{
	if (&a == &b)
	{
		add(*x.num, *a.num, *a.num);
		x.den->assign(*a.den);
	}
	else if ((*a.den) == (*b.den))
	{
		x.den->assign(*a.den);
		add(*x.num, *a.num, *b.num);
	}
	else
	{
		gf_polynomial h1, h2;
		multiply(h1, *a.num, *b.den, F);
		multiply(h2, *a.den, *b.num, F);
		multiply(*x.den, *a.den, *b.den, F);
		add(*x.num, h1, h2);
	}
}

void subtract(gf_rational_function &x, const gf_rational_function &a,
		const gf_rational_function &b)
{
	if (&a == &b)
	{
		x.num->assign_zero();
		x.den->assign_one();
	}
	else if ((*a.den) == (*b.den))
	{
		x.den->assign(*a.den);
		subtract(*x.num, *a.num, *b.num);
	}
	else
	{
		gf_polynomial h1, h2;

		multiply(h1, *a.num, *b.den);
		multiply(h2, *a.den, *b.num);
		multiply(*x.den, *a.den, *b.den);
		subtract(*x.num, h1, h2);
	}
}

void subtract_mod(gf_rational_function &x, const gf_rational_function &a,
		const gf_rational_function &b, const gf_polynomial &f)
{
	if (&a == &b)
	{
		x.num->assign_zero();
		x.den->assign_one();
	}
	else if ((*a.den) == (*b.den))
	{
		x.den->assign(*a.den);
		subtract(*x.num, *a.num, *b.num);
	}
	else
	{
		gf_polynomial h1, h2;
		multiply_mod(h1, *a.num, *b.den, f);
		multiply_mod(h2, *a.den, *b.num, f);
		multiply_mod(*x.den, *a.den, *b.den, f);
		subtract(*x.num, h1, h2);
	}
}

void subtract(gf_rational_function &x, const gf_rational_function &a,
		const gf_rational_function &b, const gf_poly_modulus &F)
{
	if (&a == &b)
	{
		x.num->assign_zero();
		x.den->assign_one();
	}
	else if ((*a.den) == (*b.den))
	{
		x.den->assign(*a.den);
		subtract(*x.num, *a.num, *b.num);
	}
	else
	{
		gf_polynomial h1, h2;
		multiply(h1, *a.num, *b.den, F);
		multiply(h2, *a.den, *b.num, F);
		multiply(*x.den, *a.den, *b.den, F);
		subtract(*x.num, h1, h2);
	}
}

void divide(gf_rational_function &q, const gf_rational_function &a,
		const gf_rational_function &b)
{
	if (&a == &b)
	{
		q.num->assign_one();
		q.den->assign_one();
	}
	else if (&q == &b)
	{
		gf_rational_function r(*b.den, *b.num);
		multiply(q, a, r);
	}
	else
	{
		multiply(*q.num, *a.num, *b.den);
		multiply(*q.den, *a.den, *b.num);
	}
}

void divide(gf_rational_function &q, const gf_polynomial &a,
		const gf_rational_function &b)
{
	if (&q == &b)
	{
		gf_rational_function r(*b.den, *b.num);
		multiply(q, r, a);
	}
	else
	{
		q.den->assign(*b.num);
		multiply(*q.num, a, *b.den);
	}
}

void invert(gf_rational_function &c, const gf_rational_function &a)
{
	if (&a == &c)
	{
		gf_rational_function r(*a.den, *a.num);
		c.assign(r);
	}
	else
	{
		c.num->assign(*a.den);
		c.den->assign(*a.num);
	}
}

void gf_rational_function::invert()
{
	gf_rational_function c(*this);

	num->assign(c.denominator());
	den->assign(c.numerator());
}

//************ Miscellaneous functions ************************

void derivative(gf_rational_function &f, const gf_rational_function &g)
{
	square(*f.den, *g.den);
	subtract(
			*f.num, derivative(*g.num) * (*g.den), derivative(*g.den) * (*g.num));
}

// not doing i/o for now

#ifdef LIDIA_NAMESPACE
} // end of namespace LiDIA
#endif

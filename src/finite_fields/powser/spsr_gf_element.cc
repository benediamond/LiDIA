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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "LiDIA/gf_element.h"
#include "LiDIA/sparse_power_series.h"

#include "LiDIA/finite_fields/base_sparse_power_series.h"
#include "LiDIA/finite_fields/base_sparse_power_series.cc"
#include "LiDIA/finite_fields/coeff_sparse_power_series.h"
#include "LiDIA/finite_fields/coeff_sparse_power_series.cc"

#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#endif

template class base_sparse_power_series<gf_element>;

template class spc<gf_element>;
template std::istream &operator>>(std::istream &in, spc<gf_element> &c);
template std::ostream &operator<<(std::ostream &out, const spc<gf_element> &c);
template void swap(spc<gf_element> &a, spc<gf_element> &b);
template int coeff_cmp_zero(const spc<gf_element> &a, const spc<gf_element> &b);

//
// ****  constructor/destructor functions    ******
//

sparse_power_series<gf_element>::sparse_power_series()
		: base_sparse_power_series<gf_element>()
{
	debug_handler("sparse_power_series< gf_element >",
			"sparse_power_series< gf_element > ()");
}

sparse_power_series<gf_element>::sparse_power_series(
		const gf_element &a, lidia_size_t l)
		: base_sparse_power_series<gf_element>(a, l)
{
	debug_handler("sparse_power_series< gf_element >",
			"sparse_power_series< gf_element > (const gf_element&, lidia_size_t)");
}

sparse_power_series<gf_element>::sparse_power_series(
		const base_vector<gf_element> &a, lidia_size_t f)
		: base_sparse_power_series<gf_element>(a, f)
{
	debug_handler("sparse_power_series< gf_element >",
			"sparse_power_series< gf_element > (const base_vector< gf_element > "
			"&, lidia_size_t)");
}

sparse_power_series<gf_element>::sparse_power_series(
		const base_sparse_power_series<gf_element> &x)
		: base_sparse_power_series<gf_element>(x)
{
	debug_handler("sparse_power_series< gf_element >",
			"sparse_power_series< gf_element > (const "
			"base_sparse_power_series< gf_element > &)");
}

#if 0
sparse_power_series<gf_element>::
~sparse_power_series ()
{
	debug_handler ("sparse_power_series< gf_element >",
		       "~sparse_power_series< gf_element > ()");
}
#endif

//
//  *****  assignment operator  *****
//

sparse_power_series<gf_element> &sparse_power_series<gf_element>::operator=(
		const base_sparse_power_series<gf_element> &x)
{
	debug_handler("sparse_power_series< gf_element >",
			"operator = (const base_sparse_power_series< gf_element > &");

	if (&x != this)
		base_sparse_power_series<gf_element>::operator=(
				static_cast<const base_sparse_power_series<gf_element> &>(x));

	return (*this);
}

//
// ***** arithmetical procedures *****
//

void sparse_power_series<gf_element>::multiply(
		const sparse_power_series<gf_element> &a,
		const sparse_power_series<gf_element> &b)
{
	debug_handler("sparse_power_series< gf_element >",
			"multiply(const sp_pow &, const sp_pow &)");

	a.init_test("multiply (3 x sp_pow)");
	b.init_test("multiply (3 x sp_pow)");

	a.sort_test();
	b.sort_test();

	lidia_size_t cf, cl, i, az, bz, ix, tmp_e;
	lidia_size_t max_size;
	lidia_size_t all, asz, bsz;
	gf_element tmp, tmp_c;

	sort_vector<spc<gf_element> > *ac, *bc;
	base_vector<gf_element> cc;

	cf = a.first + b.first;
	cl = (a.last + b.first) < (a.first + b.last) ? (a.last + b.first)
																							 : (a.first + b.last);

	if (a.is_zero() || b.is_zero())
	{
		assign_zero(cl);
	}
	else
	{
		// both series contain non-zero coefficients
		max_size = cl - cf + 1;
		asz = a.coeff->size();
		bsz = b.coeff->size();

		ac = a.coeff;
		bc = b.coeff;

		cc.set_capacity(max_size);		 // tmp. storage for coeff.
		for (i = 0; i < max_size; i++) // init. with 0
			cc[i].assign_zero();

		for (az = 0; az < asz; az++)
		{
			for (bz = 0; bz < bsz; bz++)
			{
				tmp_e = (*ac)[az].exp + (*bc)[bz].exp;

				if (tmp_e <= cl)
				{
					// coeff will be valid in result
					LiDIA::multiply(tmp_c, (*ac)[az].coeff, (*bc)[bz].coeff);
					LiDIA::add(cc[tmp_e - cf], cc[tmp_e - cf], tmp_c);
				}
				else
					bz = bsz;
			}
		}

		// --- counting number of non-zero elts. in product ---

		for (i = 0, all = 0; i < max_size; i++)
		{
			if (!cc[i].is_zero())
				all++;
		}

		if (all > 0)
		{
			coeff->set_capacity(all);
			coeff->set_size(all);

			for (i = 0, ix = 0; i <= (cl - cf); i++)
			{
				if (!cc[i].is_zero())
				{
					(*coeff)[ix].coeff = cc[i];
					(*coeff)[ix].exp = i + cf;

					ix++;
				}
			}

			first = (*coeff)[0].exp;
			last = cl;
		}
	}
}

void sparse_power_series<gf_element>::invert(
		const sparse_power_series<gf_element> &a)
{
	debug_handler("sparse_power_series< gf_element >",
			"invert (sp_pow & , const sp_pow &)");

	a.init_test("invert (2 x sp_pow)");
	a.sort_test();

	lidia_size_t i, j;
	lidia_size_t cf, cl, cx;
	lidia_size_t all, aalloc;

	gf_element bc0, ma0, tb, tmp;

	base_vector<spc<gf_element> > bc;
	sort_vector<spc<gf_element> > *ac;

	if (!a.is_zero())
	{
		cf = -a.first;
		cl = a.last - 2 * a.first;

		all = cl - cf + 1;
		aalloc = a.coeff->size();
		ac = a.coeff;

		bc.set_capacity(all);

		for (i = 0; i < all; i++)
			bc[i].coeff.assign_zero();

		// Compute the inverse of a in bc;
		// to avoid a lot of reduction steps, bc is
		// a vector of bigint; therefore, the necessary
		// reductions have to be done by hand.

		LiDIA::invert(bc0, (*ac)[0].coeff);
		LiDIA::negate(ma0, bc0); // ma0 == -1/a_0

		bc[0].coeff = bc0;
		bc[0].exp = cf;

		tb = bc[0].coeff;

		for (i = cf + 1, all = 1; i <= cl; i++)
		{
			// a * (b=1/a) is the one-approximation = 1 + \SUM {i=0} {cl-cf} {0 * X^i}.
			// The coefficient of this serie
			// that corresponds to X^(n) equals \SUM {a.first<=i<=n+a.first, j=n-i} {a_i * b_j}.
			// This sum is used to compute the b with the greatest index that appears in
			// the sum; this is b with index i+j-a.first = i+j+cf which is stored in
			// bc[i+j+cf-cf] = bc[i+j]. Therefore, the sum a[j].exp + (i-1) of the exponents of a_(a[j].exp)
			// and b_(i-1) gives us the index of the coefficient in bc which can be updated by
			// the product of a_(a[j].exp) and b_(i-1).

			// tb is the coefficient of b=1/a with exponent 'i-1'.

			if (!tb.is_zero())
			{
				for (j = 1; (j < aalloc) && ((cx = i - 1 + (*ac)[j].exp) <= cl - cf);
						 j++)
				{
					LiDIA::multiply(tmp, tb, (*ac)[j].coeff);
					LiDIA::add(bc[cx].coeff, bc[cx].coeff, tmp);
				}
			}

			// remainder(bc[i - cf].coeff, bc[i - cf].coeff, bigmod::modulus());

			if (!bc[i - cf].coeff.is_zero())
			{
				LiDIA::multiply(bc[all].coeff, bc[i - cf].coeff, ma0);
				// LiDIA::remainder(bc[all].coeff, bc[all].coeff,
				//                  bigmod::modulus());
				bc[all].exp = i;
				tb = bc[all].coeff;

				all++;
			}
			else
				tb.assign_zero();
		}

		first = cf;
		last = cl;
		coeff->set_capacity(all);

		ac = coeff;

		for (i = 0; i < all; i++)
		{
			(*ac)[i].coeff = bc[i].coeff;
			(*ac)[i].exp = bc[i].exp;
		}
	}
	else
	{
		lidia_error_handler(
				"sparse_power_series< gf_element >", "inverting of zero failed");
	}
}

void sparse_power_series<gf_element>::square(
		const sparse_power_series<gf_element> &a)
{
	debug_handler("sparse_power_series< gf_element >",
			"square (sp_pow & , const sp_pow &)");

	a.init_test("square (sp_pow & , const sp_pow &)");
	a.sort_test();
	this->multiply(a, a);
}

void sparse_power_series<gf_element>::power(
		const sparse_power_series<gf_element> &a, long n)
{
	debug_handler("sparse_power_series< gf_element >",
			"power (sp_pow & , const sp_pow & , long)");

	a.init_test("power (sp_pow & , const sp_pow & , long)");
	a.sort_test();

	sparse_power_series<gf_element> z;

	if (n >= 0)
	{
		z = a;
	}
	else
	{
		z.invert(a);
		n = -n;
	}

	set(gf_element(1), a.last - a.first); // res = 1 * X^0

	while (n > 1)
	{
		if (n & 1)
		{
			this->multiply(*this, z);
		}

		z.square(z);

		n >>= 1;
	}

	if (n == 1)
		this->multiply(*this, z);
}

void sparse_power_series<gf_element>::divide(
		const sparse_power_series<gf_element> &a,
		const sparse_power_series<gf_element> &b)
{
	debug_handler("sparse_power_series< gf_element >",
			"divide (sp_pow &, const sp_pow &, const sp_pow &)");

	a.init_test("divide(3 x sp_pow)");
	b.init_test("divide(3 x sp_pow)");
	a.sort_test();
	b.sort_test();

	sparse_power_series<gf_element> invb;
	invb.invert(b);
	this->multiply(a, invb);
}

void sparse_power_series<gf_element>::divide(
		const gf_element &b, const sparse_power_series<gf_element> &a)
{
	debug_handler("sparse_power_series< gf_element >",
			"divide (sp_pow & , const gf_element & , const sp_pow &)");

	a.init_test("divide (sp_pow , cont gf_element & , sp_pow)");
	this->invert(a);
	base_sparse_power_series<gf_element>::multiply(b, *this);
}

#ifdef LIDIA_NAMESPACE
} // end of namespace LiDIA
#endif

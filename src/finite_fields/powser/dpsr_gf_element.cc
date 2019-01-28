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
#include "LiDIA/math_vector.h"
#include "LiDIA/dense_power_series.h"

#include "LiDIA/finite_fields/base_dense_power_series.h"
#include "LiDIA/finite_fields/base_dense_power_series.cc"

#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#endif

template class base_dense_power_series<gf_element>;

//
// ************************************
// ** dense_power_series <gf_element> *
// ************************************
//

dense_power_series<gf_element>::dense_power_series()
		: base_dense_power_series<gf_element>()
{
	debug_handler("dense_power_series< gf_element >", "dense_power_series()");
}

dense_power_series<gf_element>::dense_power_series(
		const gf_element &a, lidia_size_t l)
		: base_dense_power_series<gf_element>(a, l)
{
	debug_handler("dense_power_series< gf_element >",
			"dense_power_series(const gf_element&, lidia_size_t)");
}

dense_power_series<gf_element>::dense_power_series(
		const base_vector<gf_element> &a, lidia_size_t f)
		: base_dense_power_series<gf_element>(a, f)
{
	debug_handler("dense_power_series< gf_element >",
			"dense_power_series(const base_vector< gf_element > &, lidia_size_t)");
}

dense_power_series<gf_element>::dense_power_series(
		const base_dense_power_series<gf_element> &a)
		: base_dense_power_series<gf_element>(a)
{
	debug_handler("dense_power_series< gf_element >",
			"dense_power_series(const base_dense_power_series< gf_element > &)");
}

#if 0
dense_power_series<gf_element>::
~dense_power_series ()
{
	debug_handler ("dense_power_series< gf_element >", "~dense_power_series()");
}
#endif

// ************************************************
// ************ assignment - operator *************
// ************************************************

dense_power_series<gf_element> &dense_power_series<gf_element>::operator=(
		const dense_power_series<gf_element> &a)
{
	debug_handler("dense_power_series< gf_element >",
			"operator = (const dense_power_series< gf_element > &)");

	if (&a != this)
	{
		base_dense_power_series<gf_element>::operator=(
				static_cast<const base_dense_power_series<gf_element> &>(a));
	}
	return *this;
}

dense_power_series<gf_element> &dense_power_series<gf_element>::operator=(
		const sparse_power_series<gf_element> &a)
{
	debug_handler("dense_power_series< gf_element >",
			"operator = (const sparse_power_series< gf_element > &)");
	base_dense_power_series<gf_element>::operator=(
			static_cast<const base_sparse_power_series<gf_element> &>(a));
	return *this;
}

// ************************************************
// ************** friend - functions **************
// ************************************************

// ************************************************
// ********** arithmetic via functions ************
// ************************************************

void dense_power_series<gf_element>::square(
		const dense_power_series<gf_element> &a)
{
	debug_handler("dense_power_series< gf_element >",
			"square(dense_power_series< gf_element > &, const "
			"dense_power_series< gf_element > &)");

	lidia_size_t nc;
	lidia_size_t pc;
	lidia_size_t i, j;
	lidia_size_t ed, n;
	lidia_size_t non_zero_index_a;
	int ident = 0;

	math_vector<gf_element> *C;
	math_vector<gf_element> *A = a.coeff;

	gf_element x, tmp1, tmp2; // purpose of x?

	if (A->size() == 0)
		lidia_error_handler("dense_power_series< gf_element "
												">::multiply(dense_power_series< gf_element "
												"> &, dense_power_series< gf_element > &, "
												"dense_power_series< gf_element > "
												"&",
				"Argument not initialized.");

	if (a.is_zero(non_zero_index_a))
		assign_zero(2 * a.last + 1);

	else
	{
		// precision and first of c

		pc = A->size() - non_zero_index_a;
		nc = 2 * (a.first + non_zero_index_a);

		// &a == this ?

		if (this != &a)
		{
			C = coeff;
		}
		else
		{
			ident = 1;
			C = new math_vector<gf_element>;
		}

		C->set_capacity(pc);

		// square a

		for (i = 0, n = 2 * non_zero_index_a; i < pc; i++, n++)
		{
			if (n & 1)
				ed = (n - 1) >> 1;
			else
				ed = (n >> 1) - 1;

			// tmp2.assign_zero();

			for (j = non_zero_index_a; j <= ed; j++)
			{
				LiDIA::multiply(tmp1, (*A)[j], (*A)[n - j]);
				LiDIA::add(tmp2, tmp2, tmp1);
			}

			(*C)[i] = tmp2;

			LiDIA::add((*C)[i], (*C)[i], (*C)[i]);

			if (!(n & 1))
			{
				LiDIA::square(x, (*A)[n >> 1]);
				LiDIA::add((*C)[i], (*C)[i], x);
			}
		}

		first = nc;
		last = nc + pc - 1;

		// copy result if necessary

		if (ident)
		{
			delete coeff;
			coeff = C;
		}

	} // end else a.is_zero (...)
}

void dense_power_series<gf_element>::multiply(
		const dense_power_series<gf_element> &a,
		const dense_power_series<gf_element> &b)
{
	debug_handler("dense_power_series< gf_element >",
			"multiply(dense_power_series< gf_element > &, "
			"dense_power_series< gf_element > "
			"&, dense_power_series< gf_element > &");

	if (&a == &b)
	{
		this->square(a);
	}
	else
	{
		lidia_size_t nc, pc;
		lidia_size_t i, j, n;
		lidia_size_t non_zero_index_a;
		lidia_size_t non_zero_index_b;
		bool zero_a;
		bool zero_b;
		int ident = 0;

		math_vector<gf_element> *A = a.coeff;
		math_vector<gf_element> *B = b.coeff;
		math_vector<gf_element> *C;

		gf_element tmp1, tmp2;

		if (A->size() == 0 || B->size() == 0)
			lidia_error_handler("dense_power_series< gf_element "
													">::multiply(dense_power_series< "
													"gf_element > &, dense_power_series< gf_element > &, "
													"dense_power_series< gf_element > &",
					"Arguments not initialized.");

		zero_a = a.is_zero(non_zero_index_a);
		zero_b = b.is_zero(non_zero_index_b);

		if (zero_a && zero_b)
			assign_zero(a.last + b.last + 1);

		else if (zero_a)
			assign_zero(a.last + b.first + non_zero_index_b);

		else if (zero_b)
			assign_zero(b.last + a.first + non_zero_index_a);

		else
		{
			// precision and first of c

			pc = A->size() - non_zero_index_a < B->size() - non_zero_index_b
							 ? A->size() - non_zero_index_a
							 : B->size() - non_zero_index_b;

			nc = a.first + b.first + non_zero_index_a + non_zero_index_b;

			// c == a or c == b ?

			if ((this != &a) && (this != &b))
			{
				C = coeff;
			}
			else
			{
				ident = 1;
				C = new math_vector<gf_element>;
			}

			C->set_capacity(pc);

			// multiply a and b

			for (i = 0, n = non_zero_index_a + non_zero_index_b; i < pc; i++, n++)
			{
				tmp2.assign_zero();

				for (j = non_zero_index_a; j <= n - non_zero_index_b; j++)
				{
					LiDIA::multiply(tmp1, (*A)[j], (*B)[n - j]);
					LiDIA::add(tmp2, tmp2, tmp1);
				}

				(*C)[i] = tmp2;
			}

			first = nc;
			last = nc + pc - 1;

			// copy result if necessary

			if (ident)
			{
				delete coeff;
				coeff = C;
			}

		} // end - else if (zero_a && zero_b)

	} // end - else if (&a == &b)
}

void dense_power_series<gf_element>::invert(
		const dense_power_series<gf_element> &a)
{
	debug_handler("dense_power_series< gf_element >",
			"invert(dense_power_series< gf_element > &, const "
			"dense_power_series< gf_element > &)");

	lidia_size_t non_zero_index_a;
	lidia_size_t nc;
	lidia_size_t pc;
	lidia_size_t i, j, k;
	int ident = 0;

	gf_element y, tmp1, tmp2;

	math_vector<gf_element> *C;
	math_vector<gf_element> *A = a.coeff;

	// check for invalid argument

	if (A->size() == 0)
	{
		lidia_error_handler("dense_power_series< gf_element >::invert"
												"(dense_power_series< gf_element > &, "
												"dense_power_series< gf_element > &)",
				"Argument not initialized.");
	}

	// division by zero ?

	if (a.is_zero(non_zero_index_a))
	{
		lidia_error_handler("dense_power_series< gf_element >::invert"
												"(dense_power_series< gf_element > &, "
												"dense_power_series< gf_element > &)",
				"Division by zero.");
	}
	else
	{
		// precision and first of c

		pc = A->size() - non_zero_index_a;
		nc = -(a.first + non_zero_index_a);

		// &a == this ?

		if (this != &a)
			C = coeff;
		else
		{
			ident = 1;
			C = new math_vector<gf_element>;
		}

		C->set_capacity(pc);

		// invert a

		LiDIA::invert((*C)[0], (*A)[non_zero_index_a]);
		LiDIA::negate(y, (*C)[0]);

		for (i = 1; i < pc; i++)
		{
			// tmp2.assign_zero();

			for (j = 1, k = 1 + non_zero_index_a; j <= i; j++, k++)
			{
				LiDIA::multiply(tmp1, (*A)[k], (*C)[i - j]);
				LiDIA::add(tmp2, tmp2, tmp1);
			}

			(*C)[i] = tmp2;
			LiDIA::multiply((*C)[i], (*C)[i], y);
		}

		first = nc;
		last = nc + pc - 1;

		// copy C if necessary

		if (ident)
		{
			delete coeff;
			coeff = C;
		}

	} // end else a.is_zero (...)
}

void dense_power_series<gf_element>::power(
		const dense_power_series<gf_element> &a, long n)
{
	debug_handler("dense_power_series< gf_element >",
			"power(dense_power_series< gf_element > &, const "
			"dense_power_series< gf_element > &, long)");

	dense_power_series<gf_element> z;
	lidia_size_t non_zero_index_a;
	bool zero_a;

	// check for invalid argument

	if ((a.coeff)->size() == 0)
	{
		lidia_error_handler("dense_power_series< gf_element "
												">::power(dense_power_series< gf_element > "
												"&, dense_power_series< gf_element > &, long)",
				"Argument not initialized.");
	}

	zero_a = a.is_zero(non_zero_index_a);

	if (n == 0)
	{
		if (zero_a)
			assign_one(0);
		else
			assign_one(a.last - (a.first + non_zero_index_a));
	}
	else if (zero_a)
		assign_zero(static_cast<lidia_size_t>(n * a.last));
	else
	{
		// initialize z which holds the squares

		if (n < 0)
		{
			z.invert(a);
			n = -n;
		}
		else
			z = a;

		assign_one((a.coeff)->size() - 1 - non_zero_index_a);

		// repeated squaring

		while (n > 1)
		{
			// n odd

			if (n & 1)
				this->multiply(*this, z);

			z.square(z);

			// divide n by 2

			n = n >> 1;
		}

		if (n == 1)
			this->multiply(*this, z);
	}
}

void dense_power_series<gf_element>::divide(
		const dense_power_series<gf_element> &a,
		const dense_power_series<gf_element> &b)
{
	debug_handler("dense_power_series< gf_element >",
			"divide(dense_power_series< gf_element > &, "
			"const dense_power_series< gf_element > &, const "
			"dense_power_series< gf_element > &");

	dense_power_series<gf_element> inv_b;
	inv_b.invert(b);
	multiply(a, inv_b);
}

void dense_power_series<gf_element>::divide(
		const gf_element &b, const dense_power_series<gf_element> &a)
{
	debug_handler("dense_power_series< gf_element >",
			"divide"
			"(dense_power_series< gf_element > &, gf_element&, "
			"dense_power_series< gf_element > &)");

	dense_power_series<gf_element> d;
	d.invert(a);
	base_dense_power_series<gf_element>::multiply(b, d);
}

#ifdef LIDIA_NAMESPACE
} // end of namespace LiDIA
#endif

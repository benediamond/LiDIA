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
//	Author	: Thomas Pfahler (TPf), Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include "LiDIA/gf_polynomial.h"
#include "LiDIA/finite_fields/Fp_polynomial_fft.h" // why?
#include "LiDIA/finite_fields/Fp_polynomial_util.h"
#include "LiDIA/factorization.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

//
// HB: computes a root of f over GF(p^m)
//
gf_element
find_root( const gf_polynomial & f )
{
// f is monic, and has deg(f) distinct roots.
// returns  roots

    if (f.degree() == 0)
		lidia_error_handler( "gf_pol_util::find_root(const gf_polynomial)", "polynomial is constant." );
	
    if (f.degree() == 1)
		return( - f[ 0 ] );

    gf_polynomial h, local_f( f );

    const galois_field & K = f.get_field();
    const bigint &p = K.characteristic();

    if( p != 2 )
    {
		gf_element one(K), r(K);
		one.assign_one();
		const bigint & q = K.number_of_elements();
		bigint q1( q );
		q1.divide_by_2();
		
		while( local_f.degree() > 1 )
		{
			do
			{
				r.randomize();
				power_x_plus_a_mod(h, r, q1, local_f);
				subtract(h, h, one);
				gcd(h, h, local_f);
			}
			while( h.degree() <= 0 || h.degree() == local_f.degree() );

			if( h.degree() <= local_f.degree() / 2 )
				local_f.assign( h );
			else
				divide( local_f, local_f, h );
		}
    }
    else
		lidia_error_handler( "gf_pol_util::find_root(const gf_polynomial)", "case p = 2 not implemented." );

	return( - local_f[ 0 ] );
}


void
find_roots(base_vector< gf_element > &x, const gf_polynomial & f)
{
// f is monic, and has deg(f) distinct roots.
// returns the list of roots
	debug_handler("gf_pol_util.cc", "find_roots(base_vector< gf_element > &, gf_polynomial&)");

	if (f.degree() == 0)
		return;

	if (f.degree() == 1) {
		lidia_size_t k = x.size();

		if (x.capacity() < k + 1)
			x.set_capacity(k + 1);
		x.set_size(k + 1);

		negate(x[k], f.const_term());
		return;
	}

	gf_polynomial h;

	const galois_field &K = f.get_field();
	const bigint &p = K.characteristic();

	if (p != 2) {
		gf_element one(K), r(K);
		one.assign_one();
		const bigint &q = K.number_of_elements();
		bigint q1(q);
		q1.divide_by_2();

		{
			do {
				r.randomize();
				power_x_plus_a_mod(h, r, q1, f);
				subtract(h, h, one);
				gcd(h, h, f);
			} while (h.degree() <= 0 || h.degree() == f.degree());
		}
	}
	else {
		// p == 2
		do {
			lidia_size_t i, k = K.degree(), n = f.degree();
			gf_polynomial g;
			g = randomize(K, n-1);
			h = g;
			for (i = 1; i < k; i++) {
				square(g, g, f);
				add(h, h, g);
			}
			gcd(h, h, f);
		} while (h.degree() <= 0 || h.degree() == f.degree());
	}

	find_roots(x, h);
	divide(h, f, h);
	find_roots(x, h);
}


//
// Task:	Returns the list of roots of f (without multiplicities)
//
// Conditions:	always: f mustn't be the zero polynomial
//		if (flag == 0): no assumptions [default]
//		if (flag != 0): f monic with deg(f) distinct roots.
//

base_vector<gf_element> find_roots(const gf_polynomial &f, int flag)
{
	base_vector< gf_element > v;
	if (flag != 0) { // this is what berl and can-zass call. the original code
		find_roots(v, f);
		lidia_size_t i;
		gf_polynomial h, g;
		h.assign_x(f.get_field());
		g.assign_one(f.get_field());
		for (i = 0; i < v.size(); i++) {
			h[0] = -v[i];
			multiply(g, g, h);
		}
		if (g != f)
			std::cout << "Mist." << v << std::endl << f << std::endl << g << std::endl;
	} else { // adapted from factoring.cc
		polynomial<gf_element> x_to_the_q, x;
		polynomial<gf_element> g(f);
		divide(g, g, g.lead_coeff()); // g.make_monic();

		x.assign_x(g.get_field());
		gf_poly_modulus F(g);

		// this is the costly step we want to avoid if flag != 0
		power_x(x_to_the_q, g.get_field().number_of_elements(), F);

		gcd(g, g, x_to_the_q - x); // split off linear factors

		factorization<polynomial<gf_element> > sqfr;
		square_free_decomp(sqfr, g); // decompose acc. to multiplicities

		// any component of sqfr is a product
		// of distinct linear factors
		lidia_size_t i;
		for (i = 0; i < sqfr.no_of_composite_components(); i++)
			find_roots(v, sqfr.composite_base(i).base());
	}
	return v;
}



void compose(gf_polynomial& x, const gf_polynomial& g,
	     const gf_polynomial& h, const gf_poly_modulus& F)
// x = g(h) mod f
{
	lidia_size_t m = square_root(g.degree() + 1);
	if (m == 0) {
		x.assign_zero(h.get_field());
		return;
	}

	gf_poly_argument A;
	A.build(h, F, m);
	A.compose(x, g, F);
}



void compose2(gf_polynomial& x1, gf_polynomial& x2,
	      const gf_polynomial& g1, const gf_polynomial& g2,
	      const gf_polynomial& h, const gf_poly_modulus& F)
	// xi = gi(h) mod f (i=1,2)
	// ALIAS RESTRICTION:  xi may not alias gj, for i != j
{
	lidia_size_t m = square_root(g1.degree() + g2.degree() + 2);
	if (m == 0) {
		x1.assign_zero(h.get_field());
		x2.assign_zero(h.get_field());
		return;
	}

	gf_poly_argument A;
	A.build(h, F, m);
	A.compose(x1, g1, F);
	A.compose(x2, g2, F);
}



void compose3(gf_polynomial& x1, gf_polynomial& x2, gf_polynomial& x3,
	      const gf_polynomial& g1, const gf_polynomial& g2, const gf_polynomial& g3,
	      const gf_polynomial& h, const gf_poly_modulus& F)
	// xi = gi(h) mod f (i=1..3)
	// ALIAS RESTRICTION:  xi may not alias gj, for i != j
{
	lidia_size_t m = square_root(g1.degree() + g2.degree() + g3.degree() + 3);
	if (m == 0) {
		x1.assign_zero(h.get_field());
		x2.assign_zero(h.get_field());
		x3.assign_zero(h.get_field());
		return;
	}

	gf_poly_argument A;
	A.build(h, F, m);
	A.compose(x1, g1, F);
	A.compose(x2, g2, F);
	A.compose(x3, g3, F);
}



void
trace_map(gf_polynomial & w, const gf_polynomial & a, lidia_size_t d,
	  const gf_poly_modulus & F, const gf_polynomial & b)
	// w = a+a^q+...+^{q^{d-1}} mod f;
	// it is assumed that d >= 0, and b = X^q mod f, q a power of p
{
	debug_handler("gf_polynomial", "trace_map(gf_polynomial&, gf_polynomial&, lidia_size_t, gf_poly_modulus&, gf_polynomial&)");

	w.ffield = (gf_polynomial::common_field(a.ffield, b.ffield));
	gf_polynomial::build_frame(w.ffield);


	gf_polynomial z(b), y(a), t;
	w.assign_zero();

	while (d) {
		if (d == 1) {
			if (w.is_zero())
				w.assign(y);
			else {
				compose(w, w, z, F);
				add(w, w, y);
			}
		}
		else {
			if ((d & 1) == 0) {
				compose2(z, t, z, y, z, F);
				add(y, t, y);
			}
			else {
				if (w.is_zero()) {
					w.assign(y);
					compose2(z, t, z, y, z, F);
					add(y, t, y);
				}
				else {
					compose3(z, t, w, z, y, w, z, F);
					add(w, w, y);
					add(y, t, y);
				}
			}
		}
		d = d >> 1;
	}

	gf_polynomial::delete_frame();
}


void inner_product(gf_element& x, const base_vector< gf_element > & a,
		   const gf_polynomial &b, lidia_size_t offset)
{
	debug_handler("gf_pol_util.cc", "inner_product(gf_element&, base_vector< gf_element > &, gf_polynomial&, lidia_size_t)");

	lidia_size_t n = comparator<lidia_size_t>::min(a.size(), b.degree() + 1 + offset);
	lidia_size_t i;
	gf_element accum, t;

	for (i = offset; i < n; i++) {
		multiply(t, a[i], b[i - offset]);
		add(accum, accum, t);
	}
	x.assign(accum); // remainder(x, accum, b.modulus());
}

void power_compose(gf_polynomial &y, const gf_polynomial &h, lidia_size_t q,
                   const gf_poly_modulus &F) {
	// w = X^{q^d} mod f;
	// it is assumed that d >= 0, and b = X^q mod f, q a power of p
	debug_handler("gf_polynomial",
                "power_compose(gf_polynomial&, gf_polynomial&, lidia_size_t, "
                "gf_polynomial&)");

	// common_field(h.ffield, f.ffield)
	// gf_polynomial::build_frame(y.ffield);

	gf_polynomial z(h);
	lidia_size_t sw;

	y.assign_x();

	while (q) {
		sw = 0;

		if (q > 1)
			sw = 2;
		if (q & 1) {
			if (y.is_x())
				y.assign(z);
			else
				sw = sw | 1;
		}

		switch (sw) {
		case 0:
			break;

		case 1:
			compose(y, y, z, F);
			break;

		case 2:
			compose(z, z, z, F);
			break;

		case 3:
			compose2(y, z, y, z, z, F);
			break;
		}

		q = q >> 1;
	}
}

static void tandem_power_compose(gf_polynomial &y1, gf_polynomial &y2,
                                 const gf_polynomial &h, lidia_size_t q1,
                                 lidia_size_t q2, const gf_poly_modulus &F) {
    debug_handler(
        "gf_pol_util.cc",
        "tandem_power_compose(gf_polynomial&, gf_polynomial&, gf_polynomial& "
        "h, lidia_size_t, lidia_size_t, gf_poly_modulus&)");

    y1.assign_x(h.get_field());
    y2.assign_x(h.get_field());
    // gf_polynomial::build_frame(ffield);

    gf_polynomial z(h);
    lidia_size_t sw;

    while (q1 || q2) {
        sw = 0;

        if (q1 > 1 || q2 > 1)
            sw = 4;

        if (q1 & 1) {
            if (y1.is_x())
                y1.assign(z);
            else
                sw = sw | 2;
        }

        if (q2 & 1) {
            if (y2.is_x())
                y2.assign(z);
            else
                sw = sw | 1;
        }

        switch (sw) {
        case 0:
            break;

        case 1:
            compose(y2, y2, z, F);
            break;

        case 2:
            compose(y1, y1, z, F);
            break;

        case 3:
            compose2(y1, y2, y1, y2, z, F);
            break;

        case 4:
            compose(z, z, z, F);
            break;

        case 5:
            compose2(z, y2, z, y2, z, F);
            break;

        case 6:
            compose2(z, y1, z, y1, z, F);
            break;

        case 7:
            compose3(z, y1, y2, z, y1, y2, z, F);
            break;
        }

        q1 = q1 >> 1;
        q2 = q2 >> 1;
    }

    // gf_polynomial::delete_frame();
}

//***********************************************************************
//
//  Algorithms which are useful in counting points on elliptic curves:
//		compute_degree, prob_compute_degree
//
//***********************************************************************

static lidia_size_t base_case(const gf_polynomial &h, lidia_size_t q, int a,
                              const gf_poly_modulus &F) {
    debug_handler(
        "gf_pol_util.cc",
        "base_case(gf_polynomial&, lidia_size_t, int, gf_poly_modulus&)");

    lidia_size_t b;
    gf_polynomial lh(h);
    // lh.set_max_degree(F.modulus().degree() - 1);

    b = 0;
    while (b < a - 1 && !lh.is_x()) {
        b++; // no binary search? see Shoup, Frobenius, 1992, p. 23
        power_compose(lh, lh, q, F);
    }

    if (!lh.is_x())
        b++;

    return b;
}

static int compute_split(int lo, int hi, const fac_vec &fvec) {
    // duplicated code, but is it worth bothering trying to set up an import?
    debug_handler("gf_pol_util.cc", "compute_split(int, int, fac_vec&)");

    int mid, i;
    double total, sum;

    total = 0;
    for (i = lo; i <= hi; i++)
        total = total + fvec[i].len;

    mid = lo - 1;
    sum = 0;
    while (sum < total / 2) {
        mid++;
        sum = sum + fvec[mid].len;
    }

    if (mid == hi || (mid != lo && 2 * sum > total + fvec[mid].len))
        mid--;

    return mid;
}

static void rec_compute_degree(int lo, int hi, const gf_polynomial &h,
                               const gf_poly_modulus &F, fac_vec &fvec) {
    debug_handler("gf_pol_util.cc",
                  "rec_compute_degree(int, int, gf_polynomial&, "
                  "gf_poly_modulus&, fac_vec&)");

    int mid;
    lidia_size_t q1, q2;
    gf_polynomial h1, h2;

    if (h.is_x()) {
        fvec.clear(lo, hi);
        return;
    }

    if (lo == hi) {
        fvec[lo].b = base_case(h, fvec[lo].q, fvec[lo].a, F);
        return;
    }

    mid = compute_split(lo, hi, fvec);

    q1 = fvec.prod(lo, mid);
    q2 = fvec.prod(mid + 1, hi);

    tandem_power_compose(h1, h2, h, q1, q2, F);
    rec_compute_degree(lo, mid, h2, F, fvec);
    rec_compute_degree(mid + 1, hi, h1, F, fvec);
}

//
// Task:	The common degree of the irreducible factors of f is computed.
//		This routine is useful in counting points on elliptic curves.
//
// Conditions:	f = F.modulus() is assumed to be an "equal degree" polynomial,
//		h = x^p mod f,
//              d is multiple of common degree, if d == -1, such information
//              not known
//

lidia_size_t compute_degree(const gf_polynomial &h, const gf_poly_modulus &F,
                            lidia_size_t d) {
    debug_handler("gf_pol_util.cc",
                  "compute_degree(gf_polynomial&, gf_poly_modulus&)");

    // h.comp_modulus(F.modulus(), "compute_degree");

    if (h.is_x())
        return 1;

    lidia_size_t res;

    if (d == -1)
        d = F.modulus().degree();

    fac_vec fvec(d);

    int i, NumFactors = fvec.number_of_factors();

    rec_compute_degree(0, NumFactors - 1, h, F, fvec);

    res = 1;

    for (i = 0; i < NumFactors; i++) {
        res = res * static_cast<lidia_size_t>(
                        power(static_cast<udigit>(fvec[i].q),
                              static_cast<unsigned int>(fvec[i].b)));
    }

    return res;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

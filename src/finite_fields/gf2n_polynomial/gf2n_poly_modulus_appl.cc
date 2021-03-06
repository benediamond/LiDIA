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


#include	"LiDIA/gf2n_poly_modulus.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/timer.h"

#include	<cassert>
#include        <cstdlib>


#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	lidia_size_t n, d, tests, i;
	bigint q;
	timer t1, t2;
	t1.set_print_mode(HMS_MODE);
	t2.set_print_mode(HMS_MODE);

	std::cout << "\n Input Field Degree over GF(2):  ";
	std::cin >> n;

	gf2n_init(n);
	shift_left(q, bigint(1), n);

	gf2nIO::noprefix();

	std::cout << "\n Polynomial degree : ";
	std::cin >> d;

	std::cout << "\n Number of Tests : "; std::cin >> tests;

	gf2n_polynomial f;
	gf2n_poly_modulus fpm;
	gf2n_polynomial erg1, erg2, x, xq;

	std::cout << "\n\n";

	for (i = 1; i <= tests; i++) {
		f.randomize(d);
		f.make_monic();
		fpm.build(f);

		x.assign_x();
		erg1.assign_x();

		t1.start_timer();
		power(xq, x, q, fpm);
		t1.stop_timer();
		std::cout << "\nPower needed time " << t1 << std::flush;
		erg1.assign(x);

		t2.start_timer();
		Xq(erg1, fpm);
		t2.stop_timer();
		std::cout << "\nXq    needed time " << t2 << std::flush;

		if (xq != erg1)
			std::exit(1);
		std::cout << "\n Test " << i << " successful " << std::flush;
	}

	std::cout << "\n\n";
	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}

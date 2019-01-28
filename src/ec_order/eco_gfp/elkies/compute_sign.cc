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


//  Uses the Dewaghe idea top compute the sign of the eigenvalue



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void eco_prime::sign_dewaghe (lidia_size_t & alpha,
			      const ff_pol & fC)
{
	if (l % 4 == 1) {
		lidia_error_handler("eco_prime", "sign_dewaghe::Dewaghe's idea only "
				    "usable for l = 3 mod 4 !!");
		return;
	}

	if (alpha < 0)
		alpha = l - alpha;

	ff_pol curve(A.get_field());
	ff_element h;

	curve.set_coefficient(3);
	curve.set_coefficient(A, 1);
	curve.set_coefficient(B, 0);

	resultant(h, fC, curve);

	if (h.is_square() != (jacobi(static_cast<udigit>(alpha), static_cast<udigit>(l)) == 1)) {
		alpha = -alpha;
	}
	if (alpha < 0)
		alpha = l + alpha;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

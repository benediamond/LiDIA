#include <LiDIA/bigint.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/meq_prime.h>
#include <LiDIA/bigmod.h>

using namespace std;
using namespace LiDIA;

int main(int argc, const char * argv[]) {
  meq_prime meq;
  cout << "success? " << meq.set_prime((unsigned long) 7) << endl;

  galois_field K(17, 1);
  gf_element::set_uninitialized_field(K);
  gf_element c(K);
  c.assign(4);
  polynomial< gf_element > f(K);
  meq.build_poly_in_Y(f, c);
  cout << f << endl;

  return 0;
}

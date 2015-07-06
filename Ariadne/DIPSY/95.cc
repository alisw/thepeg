#include <iostream>
#include <map>
#include <cmath>

using namespace std;

int main() {

  const double dummy = sqrt(dummy);

  multimap<double,double> bdist;
  multimap<double,double> fdist;
  double w = 0.0;
  double b = 0.0;
  double n = 0.0;
  double f = 0.0;
  double m = 0.0;
  double sumw = 0.0;
  double sumnw = 0.0;
  double sumnw2 = 0.0;

  while ( cin >> w >> b >> n >> f) {
    sumw += w;
    sumnw += n*w;
    sumnw2 += n*n*w;
    bdist.insert(make_pair(b, w));
    fdist.insert(make_pair(f, w));
  }
  cerr << "read " << bdist.size() << " lines" << endl;
  cerr << "sumw " << sumw << endl;
  cerr << "<n>  " << sumnw/sumw
       << " +- " << sqrt((sumnw2/sumw - sumnw*sumnw/(sumw*sumw))/bdist.size()) << endl;
  double w5 = 0.05*sumw;

  multimap<double,double>::iterator next = bdist.begin();
  multimap<double,double>::iterator prev = next;
  while ( ++next != bdist.end() ) {
    next->second += prev->second;
    if ( next->second > w5&& prev->second <= w5 ) {
      double b5 = prev->first + (w5 - prev->second)*
	(next->first - prev->first)/(next->second - prev->second);
      cout << "5% centrality in b < " << b5 << endl;
      break;
    }
    prev = next;
  }

  multimap<double,double>::reverse_iterator rnext = fdist.rbegin();
  multimap<double,double>::reverse_iterator rprev = rnext;
  while ( ++rnext != fdist.rend() ) {
    rnext->second += rprev->second;
    if ( rnext->second > w5 && rprev->second <= w5 ) {
      double f5 = rprev->first + (w5 - rprev->second)*
	(rnext->first - rprev->first)/(rnext->second - rprev->second);
      cout << "5% centrality for f > " << f5 << endl;
      break;
    }
    rprev = rnext;
  }
}

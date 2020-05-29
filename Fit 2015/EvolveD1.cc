#include <apfel/apfelxx.h>

// Parameters
double pars[1] = {0};

// Evolution parameters
const double mc = 1.4;
const double mb = 4.5;
const int PerturbativeOrder = 0;   // 0 = LO , 1 = NLO
const double mu0 = 1;   // Q0
const double AlphaQCDRef = 0.13939;   // alpha_S(MZ) from MSTW08: LO 0.13939 , NLO 0.12018
const double MuAlphaQCDRef = 91.1876;

// Initial scale functions
extern"C"
{
  double zd1u_  (double *x, double *mass);
  double zd1d_  (double *x, double *mass);
  double zd1s_  (double *x, double *mass);
  double zd1c_  (double *x, double *mass);
}

//_________________________________________________________________________________
std::map<int, double> Dists(double const& x, double const&) 
{
  double xx = x;

  // Call all functions only once.
  const double dwn = zd1d_ (&xx, &pars[0]);
  const double up  = zd1u_ (&xx, &pars[0]);
  const double str = zd1s_ (&xx, &pars[0]);
  const double chm = zd1c_ (&xx, &pars[0]);

  // Construct QCD evolution basis conbinations.
  double const Singlet = 2*dwn + 2*up + 2*str + 2*chm;
  double const Valence = 0;
  double const T3      = 2*up - 2*dwn;
  double const T8      = 2*up + 2*dwn - 4*str;
  double const T15     = 2*up + 2*dwn + 2*str - 6*chm;
  double const T24     = 2*up + 2*dwn + 2*str + 2*chm;
  double const Gluon   = 0;   //  D1g(Q0)=0  case 1
//  double const Gluon   = 0;   //  D1g(Q0)=D1u(Q0)/4  case 2
//  double const Gluon   = 0;   //  D1g(Q0)=D1u(Q0)  case 3

  // Fill in map in the QCD evolution basis.
  std::map<int,double> QCDEvMap;
  QCDEvMap[0]  = Gluon;
  QCDEvMap[1]  = Singlet;
  QCDEvMap[2]  = Valence;
  QCDEvMap[3]  = T3;
  QCDEvMap[4]  = Valence;
  QCDEvMap[5]  = T8;
  QCDEvMap[6]  = Valence;
  QCDEvMap[7]  = T15;
  QCDEvMap[8]  = Valence;
  QCDEvMap[9]  = T24;
  QCDEvMap[10] = Valence;
  QCDEvMap[11] = T24;
  QCDEvMap[12] = Valence;
  return QCDEvMap;
}

extern "C" 
{
  // Running coupling object
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, {0, 0, 0, mc, mb}, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 101, 3};
  double as(double const& mu) { return Alphas.Evaluate(mu); };

  // Evolution object
  const apfel::Grid g{{{100, 1e-5, 3}, {60, 1e-1, 3}, {50, 6e-1, 3}, {50, 8e-1, 3}}};
//  const auto DglapObj = apfel::InitializeDglapObjectsQCDtrans(g, {0, 0, 0, mc, mb});  // chiral-odd Space-like (PDFs)
  const auto DglapObj = apfel::InitializeDglapObjectsQCDT(g, {0, 0, 0, mc, mb});  // chiral-even Time-like (FFs)
  std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> TabulatedPDFs;

  //_____________________________________________________________________________
  void initialiseevolution_(double* p)
  {
    TabulatedPDFs = NULL;

    // Set silent mode for APFEL++
    apfel::SetVerbosityLevel(0);

    // Allocate parameters
//    for (int i = 0; i < 10; i++)
//	  pars[i] = p[i];
    pars[0] = p[0];

    // Set initial-scale distributions
    std::unique_ptr<apfel::Dglap<apfel::Distribution>> EvolvedPDFs = apfel::BuildDglap(DglapObj, Dists, mu0, PerturbativeOrder, as);

    // Tabulate distributions
    TabulatedPDFs = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>(new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*EvolvedPDFs, 50, 1, 100, 3});
  }

  //_____________________________________________________________________________
  void evolvetransversities_(double const& x, double const& mu, double* xf)
  {
    // Evolve and rotate distributions
    for (auto const& f : apfel::QCDEvToPhys(TabulatedPDFs->EvaluateMapxQ(x, mu)))
      if (f.first != 21)
	xf[f.first + 6] = f.second;
  }
}

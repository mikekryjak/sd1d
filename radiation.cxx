
#include <globals.hxx> // for mesh object
#include <output.hxx>
#include <boutexception.hxx>
#include <utils.hxx>

#include "radiation.hxx"

#include <fstream>
#include <string>
#include <sstream>

using std::string;

const Field3D RadiatedPower::power(const Field3D &Te, const Field3D &Ne, const Field3D &Ni) {
  Field3D result;
  result.allocate();
  
  for(const auto& i : result) {
    result[i] = power(Te[i], Ne[i], Ni[i]);
  }
  
  return result;
}

InterpRadiatedPower::InterpRadiatedPower(const string &filename) {
  std::ifstream file(filename.c_str());
  
  output.write("Loading data from file: %s\n", filename.c_str());

  if(!file.is_open())
    throw BoutException("InterpRadiatedPower: Couldn't open file %s\n", filename.c_str());

  string line;
  int linenr = 1;
  while( std::getline(file, line) ) {
    // Expecting either a comment, blank line, or two numbers
    // Remove comments, then whitespace from left and right
    string strippedline = trim( trimComments( line ) );
    
    if(strippedline.length() == 0)
      continue;
    
    std::stringstream ss(strippedline);
    BoutReal t, p;
    if( !(ss >> t) )
      throw BoutException("InterpRadiatedPower: file '%s' line %d: %s\n",
                          filename.c_str(), linenr, line.c_str());
    
    if( !(ss >> p) )
      throw BoutException("InterpRadiatedPower: file '%s' line %d: %s\n",
                          filename.c_str(), linenr, line.c_str());
    
    te_array.push_back(t);
    p_array.push_back(p);
  
    linenr++;
  }
  

  file.close();
}

BoutReal InterpRadiatedPower::power(BoutReal Te, BoutReal ne, BoutReal ni) {
  return 0.0;
}

////////////////////////////////////////////////////////////////
// 

BoutReal HydrogenRadiatedPower::power(BoutReal Te, BoutReal ne, BoutReal ni) {
  
}

// Collision rate coefficient <sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::ionisation(BoutReal T) {
  double fION;	//collision rate coefficient <sigma*v> [m3/s]
  double TT,X,S;

  TT=T;
  
  if (TT<1.0) TT=1.0;
  X=log10(TT);
  
  if (TT>=20.0) 
    S=-0.5151*X-2.563/X-5.231;
  else
    S=-3.054*X-15.72*exp(-X)+1.603*exp(-X*X);
  
  fION=pow(10.0,S-6.0);           
  
  return fION;
}

//<sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::recombination(BoutReal n, BoutReal T) {
  double fREC;	//<sigma*v> [m3/s]
  double TT,RDNE,RTE,DNE,E,RN,RT,RNJ,RTI,suma;
  int i,j,i1,j1;
  
  if(n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;
  
  double MATA[9][9]= {
    {-2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03, \
     -5.386830982777E-03, -1.626039237665E-04,  6.080907650243E-06, \
     2.101102051942E-05, -2.770717597683E-06,  1.038235939800E-07,}, 
    { 3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03, \
      -9.532840484460E-04,  1.888048628708E-04, -1.014890683861E-05, \
      2.245676563601E-05, -4.695982369246E-06,  2.523166611507E-07,}, 
    {-2.799644392058E-02, -7.452514292790E-03,  6.958711963182E-03, \
     4.631753807534E-04,  1.288577690147E-04, -1.145028889459E-04, \
     -2.245624273814E-06,  3.250878872873E-06, -2.145390398476E-07,},  
    { 1.209545317879E-02,  2.709299760454E-03, -2.139257298118E-03, \
      -5.371179699661E-04, -1.634580516353E-05,  5.942193980802E-05, \
      -2.944873763540E-06, -9.387290785993E-07,  7.381435237585E-08,},
    {-2.436630799820E-03, -7.745129766167E-04,  4.603883706734E-04, \
     1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05, \
     1.002105099354E-06,  1.392391630459E-07, -1.299713684966E-08,},
    { 2.837893719800E-04,  1.142444698207E-04, -5.991636837395E-05, \
      -2.257565836876E-05,  3.425262385387E-06,  1.118965496365E-06, \
      -1.291320799814E-07, -1.139093288575E-08,  1.265189576423E-09,},
    {-1.886511169084E-05, -9.382783518064E-06,  4.729262545726E-06, \
     1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08, \
     7.786155463269E-09,  5.178505597480E-10, -6.854203970018E-11,},
    { 6.752155602894E-07,  3.902800099653E-07, -1.993485395689E-07, \
      -6.618240780594E-08,  2.042041097083E-08,  3.708616111085E-10, \
      -2.441127783437E-10, -9.452402157390E-12,  1.836615031798E-12,},
    {-1.005893858779E-08, -6.387411585521E-09,  3.352589865190E-09, \
     1.013364275013E-09, -3.707977721109E-10,  7.068450112690E-12, \
     3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,},
  };
	
  RDNE=n;
  RTE=T;
	
  DNE=RDNE;
  TT=RTE;
  E=DNE*1.0E-14;
	
  if (TT<1.0) TT=1.0;
	
  RN=log(E);
  RT=log(TT);
	
  suma=0.0;
  for (i=1;i<=9;i++)
    {
      i1=i-1;
      for (j=1;j<=9;j++)
        {
          j1=j-1;
          RNJ=pow(RN,j1);
          if ((RN==0.0) && (j1==0)) RNJ=1.0;
          RTI=pow(RT,i1);
          if ((RT==0.0) && (i1==0)) RTI=1.0;		
          suma=suma+MATA[j-1][i-1]*RNJ*RTI;
        }
    }

  fREC=exp(suma)*1.0E-6;

  return fREC;
}

// <sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::chargeExchange(BoutReal Te) {
  double fCX;		//<sigma*v> [m3/s]
  double TT,S;
  
  TT=Te;
  
  if (TT<1.0) TT=1.0;    
  S=-14.0+log10(TT)/3.0;
  
  fCX=pow(10.0,S); 
  
  return fCX;
}

// <sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::excitation(BoutReal Te) {
  double fEXC;	//<sigma*v> [m3/s]
  double TT,Y;
  
  TT=Te;
  
  if (TT<1.0) TT=1.0;
  Y=10.2/TT;
  
  fEXC=49.0E-14/(0.28+Y)*exp(-Y)*sqrt(Y*(1.0+Y));
  
  return fEXC;
}


/////////////////////////////////////////////////////////////////////////////


BoutReal UpdatedRadiatedPower::power(BoutReal Te, BoutReal ne, BoutReal ni) {
  throw BoutException("UpdatedRadiatedPower::power not implemented");
}


//<sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::recombination(BoutReal n, BoutReal T) {
  double TT, RDNE, RTE, DNE, E, RN, RT, RNJ, RTI, suma, fHAV, fRAD;
  int i, j, i1, j1;

  if (n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;

  if (T < 0.025) {
    T = 0.025; // 300K
  }

  double MATA[9][9] = {
      {
          -2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03,
          -5.386830982777E-03, -1.626039237665E-04, 6.080907650243E-06,
          2.101102051942E-05, -2.770717597683E-06, 1.038235939800E-07,
      },
      {
          3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03,
          -9.532840484460E-04, 1.888048628708E-04, -1.014890683861E-05,
          2.245676563601E-05, -4.695982369246E-06, 2.523166611507E-07,
      },
      {
          -2.799644392058E-02, -7.452514292790E-03, 6.958711963182E-03,
          4.631753807534E-04, 1.288577690147E-04, -1.145028889459E-04,
          -2.245624273814E-06, 3.250878872873E-06, -2.145390398476E-07,
      },
      {
          1.209545317879E-02, 2.709299760454E-03, -2.139257298118E-03,
          -5.371179699661E-04, -1.634580516353E-05, 5.942193980802E-05,
          -2.944873763540E-06, -9.387290785993E-07, 7.381435237585E-08,
      },
      {
          -2.436630799820E-03, -7.745129766167E-04, 4.603883706734E-04,
          1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05,
          1.002105099354E-06, 1.392391630459E-07, -1.299713684966E-08,
      },
      {
          2.837893719800E-04, 1.142444698207E-04, -5.991636837395E-05,
          -2.257565836876E-05, 3.425262385387E-06, 1.118965496365E-06,
          -1.291320799814E-07, -1.139093288575E-08, 1.265189576423E-09,
      },
      {
          -1.886511169084E-05, -9.382783518064E-06, 4.729262545726E-06,
          1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08,
          7.786155463269E-09, 5.178505597480E-10, -6.854203970018E-11,
      },
      {
          6.752155602894E-07, 3.902800099653E-07, -1.993485395689E-07,
          -6.618240780594E-08, 2.042041097083E-08, 3.708616111085E-10,
          -2.441127783437E-10, -9.452402157390E-12, 1.836615031798E-12,
      },
      {
          -1.005893858779E-08, -6.387411585521E-09, 3.352589865190E-09,
          1.013364275013E-09, -3.707977721109E-10, 7.068450112690E-12,
          3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,
      },
  };

  RDNE = n;
  RTE = T;

  DNE = RDNE;
  TT = RTE;
  E = DNE * 1.0E-14;

  RN = log(E);
  RT = log(TT);

  suma = 0.0;
  for (i = 1; i <= 9; i++) {
    i1 = i - 1;
    for (j = 1; j <= 9; j++) {
      j1 = j - 1;
      RNJ = pow(RN, j1);
      if ((RN == 0.0) && (j1 == 0))
        RNJ = 1.0;
      RTI = pow(RT, i1);
      if ((RT == 0.0) && (i1 == 0))
        RTI = 1.0;
      suma = suma + MATA[j - 1][i - 1] * RNJ * RTI;
    }
  }

  fHAV = exp(suma) * 1.0E-6 / (1.0 + 0.125 * TT);

  double A = 3.92E-20;
  double B = 3.0E-124 * pow(1.6E-19, -4.5);
  double Ry = 13.60569;
  double chi = 0.35;

  fRAD = A * pow(Ry, 1.5) * 1.0 / (sqrt(TT) * (Ry + chi * TT));

  return fHAV + fRAD + B * DNE * pow(TT, -5.0);
}

// <sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::chargeExchange(BoutReal T) {
  if (T < 0.025) {
    T = 0.025; // 300K
  }

  double E = 10.0;

  double cxcoeffs[9][9] = {
      {
          -1.829079582E1, 1.640252721E-1, 3.364564509E-2, 9.530225559E-3,
          -8.519413900E-4, -1.247583861E-3, 3.014307546E-4, -2.499323170E-5,
          6.932627238E-7,
      },
      {
          2.169137616E-1, -1.106722014E-1, -1.382158680E-3, 7.348786287E-3,
          -6.343059502E-4, -1.919569450E-4, 4.075019352E-5, -2.850044983E-6,
          6.966822400E-8,
      },
      {
          4.307131244E-2, 8.948693625E-3, -1.209480567E-2, -3.675019470E-4,
          1.039643391E-3, -1.553840718E-4, 2.670827249E-6, 7.695300598E-7,
          -3.783302282E-8,
      },
      {
          -5.754895093E-4, 6.062141761E-3, 1.075907882E-3, -8.119301728E-4,
          8.911036876E-6, 3.175388950E-5, -4.515123642E-6, 2.187439284E-7,
          -2.911233952E-9,
      },
      {
          -1.552077120E-3, -1.210431588E-3, 8.297212634E-4, 1.361661817E-4,
          -1.008928628E-4, 1.080693990E-5, 5.106059414E-7, -1.299275586E-7,
          5.117133050E-9,
      },
      {
          -1.876800283E-4, -4.052878752E-5, -1.907025663E-4, 1.141663042E-5,
          1.775681984E-5, -3.149286924E-6, 3.105491555E-8, 2.274394089E-8,
          -1.130988251E-9,
      },
      {
          1.125490271E-4, 2.875900436E-5, 1.338839629E-5, -4.340802793E-6,
          -7.003521917E-7, 2.318308730E-7, -6.030983538E-9, -1.755944926E-9,
          1.005189187E-10,
      },
      {
          -1.238982763E-5, -2.616998140E-6, -1.171762874E-7, 3.517971869E-7,
          -4.928692833E-8, 1.756388999E-10, -1.446756796E-10, 7.143183138E-11,
          -3.989884106E-12,
      },
      {
          4.163596197E-7, 7.558092849E-8, -1.328404104E-8, -9.170850254E-9,
          3.208853884E-9, -3.952740759E-10, 2.739558476E-11, -1.693040209E-12,
          6.388219930E-14,
      },
  };

  double lograte = 0.0;
  for (int i = 0; i <= 8; i++) {
    for (int j = 0; j <= 8; j++) {
      lograte = lograte + cxcoeffs[i][j] * pow(log(T), i) * pow(log(E), j);
    }
  }

  return 1.0E-6 * exp(lograte);
}

// <sigma*v> [m3/s]
// ORIGINAL FUNCTION
//BoutReal UpdatedRadiatedPower::excitation(BoutReal Te) {
//  double fEXC;	//<sigma*v> [m3/s]
//  double TT,Y;
//  
//  TT=Te;
//  
//  if (TT<1.0) TT=1.0;
//  Y=10.2/TT;
//  fEXC=49.0E-14/(0.28+Y)*exp(-Y)*sqrt(Y*(1.0+Y));
//  
//  return fEXC;
//}

// MANUALLY FIT BY MK TO MATCH STEFAN MIJIN THESIS E_iz (SOLKIT)
BoutReal UpdatedRadiatedPower::excitation_old(BoutReal Te) {
  double fEXC;	//<sigma*v> [m3/s]
  double TT,Y;
  
  TT=Te;
  
  if (TT<1.0) TT=1.0;
  Y=10.2/TT;
  fEXC=1E-13*5.27370587/(1.14166254+Y)*exp(-Y*1.24326264);
  
  return fEXC;
}


// Collision rate coefficient <sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::ionisation_old(BoutReal T) {
    double fION; // Rate coefficient
    double TT;

    if (T < 0.025) {
      T = 0.025; // 300K
    }
    
    TT = T;
	
// ORIGINAL SD1D RATE
    double ioncoeffs[9] = {-3.271397E1, 1.353656E1, -5.739329, 1.563155, \
			   -2.877056E-1, 3.482560e-2, -2.631976E-3, \
			   1.119544E-4, -2.039150E-6};
    
    double lograte = 0.0;
    for (int i=0;i<=8;i++)
      {
	lograte = lograte + ioncoeffs[i]*pow(log(TT),i);
      }

    fION = exp(lograte)*1.0E-6;

    return fION;
}

// <sigma*v> [m3/s]
// COMES FROM AMJUEL H.4 2.1.5 (SAWADA)
BoutReal UpdatedRadiatedPower::ionisation(BoutReal n, BoutReal T) {
  // double TT, RDNE, RTE, DNE, E, RN, RT, RNJ, RTI, suma, fION;
  // int i, j, i1, j1;
  
  double E, suma, fION;
  int i, j;

  if (n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;

  if (T < 0.025) {
    T = 0.025; // 300K
  }

  double MATA[9][9] = {
      {
          -3.248025330340E+01, 1.425332391510E+01, -6.632235026785E+00,
          2.059544135448E+00, -4.425370331410E-01, 6.309381861496E-02,
          -5.620091829261E-03, 2.812016578355E-04, -6.011143453374E-06,
      },
      {
          -5.440669186583E-02, -3.594347160760E-02, 9.255558353174E-02,
          -7.562462086943E-02, 2.882634019199E-02, -5.788686535780E-03,
          6.329105568040E-04, -3.564132950345E-05, 8.089651265488E-07,
      },
      {
          9.048888225109E-02, -2.014729121556E-02, -5.580210154625E-03,
          1.519595967433E-02, -7.285771485050E-03, 1.507382955250E-03,
          -1.527777697951E-04, 7.222726811078E-06, -1.186212683668E-07,
      },
      {
          -4.054078993576E-02, 1.039773615730E-02, -5.902218748238E-03,
          5.803498098354E-04, 4.643389885987E-04, -1.201550548662E-04,
          8.270124691336E-06, 1.433018694347E-07, -2.381080756307E-08,
      },
      {
          8.976513750477E-03, -1.771792153042E-03, 1.295609806553E-03,
          -3.527285012725E-04, 1.145700685235E-06, 6.574487543511E-06,
          3.224101773605E-08, -1.097431215601E-07, 6.271173694534E-09,
      },
      {
          -1.060334011186E-03, 1.237467264294E-04, -1.056721622588E-04,
          3.201533740322E-05, 8.493662724988E-07, -9.678782818849E-07,
          4.377402649057E-08, 7.789031791949E-09, -5.483010244930E-10,
      },
      {
          6.846238436472E-05, -3.130184159149E-06, 4.646310029498E-06,
          -1.835196889733E-06, -1.001032516512E-08, 5.176265845225E-08,
          -2.622921686955E-09, -4.197728680251E-10, 3.064611702159E-11,
      },
      {
          -2.242955329604E-06, -3.051994601527E-08, -1.479612391848E-07,
          9.474014343303E-08, -1.476839184318E-08, 1.291551676860E-09,
          -2.259663431436E-10, 3.032260338723E-11, -1.355903284487E-12,
      },
      {
          2.890437688072E-08, 1.888148175469E-09, 2.852251258320E-09,
          -2.342505583774E-09, 6.047700368169E-10, -9.685157340473E-11,
          1.161438990709E-11, -8.911076930014E-13, 2.935080031599E-14,
      },
  };

  // RDNE = n;
  // RTE = T;

  // DNE = RDNE;
  // TT = RTE;
  // E = DNE * 1.0E-14;

  // RN = log(E);
  // RT = log(TT);

  // suma = 0.0;
  // for (i = 1; i <= 9; i++) {
    // i1 = i - 1;
    // for (j = 1; j <= 9; j++) {
      // j1 = j - 1;
      // RNJ = pow(RN, j1);
      // if ((RN == 0.0) && (j1 == 0))
        // RNJ = 1.0;
      // RTI = pow(RT, i1);
      // if ((RT == 0.0) && (i1 == 0))
        // RTI = 1.0;
      // suma = suma + MATA[j - 1][i - 1] * RNJ * RTI;
    // }
  // }
  
  E = n * 1.0E-14;
  E = log(E);
  T = log(T);
  
  suma = 0.0;
  
  for (i = 0; i <=8; i++) {
	  for (j = 0; j <= 8; j++) {
		  suma = suma + MATA[j][i] * pow(E,j) * pow(T,i);
	  }
  }

  fION = exp(suma)*1.0E-6;

  return fION;
}

// Collision rate coefficient <sigma*v> [m3/s]
// COMES FROM AMJUEL H.4 2.1.5 (SAWADA) E-index 0 (no density dependence, i.e. coronal approximation)
BoutReal UpdatedRadiatedPower::ionisation_coronal(BoutReal T) {
    double fION; // Rate coefficient
    double TT;

    if (T < 0.025) {
      T = 0.025; // 300K
    }
    
    TT = T;
	
    double ioncoeffs[9] = {-3.248025330340E+01, 1.425332391510E+01, -6.632235026785E+00, \
			   1.425332391510E+01, -6.632235026785E+00, 2.059544135448E+00, \
			   -6.632235026785E+00, 2.059544135448E+00, -4.425370331410E-01};
    
    double lograte = 0.0;
    for (int i=0;i<=8;i++)
      {
	lograte = lograte + ioncoeffs[i]*pow(log(TT),i);
      }

    fION = exp(lograte)*1.0E-6;

    return fION;
}

// COMES FROM AMJUEL H.10 2.1.5 (SAWADA)
// This is an energy weighted rate (m-3s-1eV) for energy loss due to multistep ionisation in a 9 coefficient 2D polynomial fit
// It includes energy loss due to ionisation (i.e. 13.6eV) within it, so for SD1D's definition we need to separate this out later
// this is done in sd1d.cxx
BoutReal UpdatedRadiatedPower::excitation(BoutReal n, BoutReal T) {
  double E, suma, fEXC;
  int i, j;

  if (n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;

  if (T < 0.025) {
    T = 0.025; // 300K
  }

  double MATA[9][9] = {
      {
          -2.497580168306E+01, 1.004448839974E+01, -4.867952931298E+00,
          1.689422238067E+00, -4.103532320100E-01, 6.469718387357E-02,
          -6.215861314764E-03, 3.289809895460E-04, -7.335808238917E-06,
      },
      {
          1.081653961822E-03, -3.189474633369E-03, -5.852267850690E-03,
          7.744372210287E-03, -3.622291213236E-03, 8.268567898126E-04,
          -9.836595524255E-05, 5.845697922558E-06, -1.367574486885E-07,
      },
      {
          -7.358936044605E-04, 2.510128351932E-03, 2.867458651322E-03,
          -3.087364236497E-03, 1.327415215304E-03, -2.830939623802E-04,
          3.017296919092E-05, -1.479323780613E-06, 2.423236476442E-08,
      },
      {
          4.122398646951E-04, -7.707040988954E-04, -8.328668093987E-04,
          4.707676288420E-04, -1.424078519508E-04, 2.411848024960E-05,
          -1.474253805845E-06, -4.633029022577E-08, 5.733871119707E-09,
      },
      {
          -1.408153300988E-04, 1.031309578578E-04, 2.056134355492E-04,
          -5.508611815406E-05, 3.307339563081E-06, 5.707984861100E-07,
          -2.397868837417E-07, 3.337390374041E-08, -1.512777532459E-09,
      },
      {
          2.469730836220E-05, -3.716939423005E-06, -3.301570807523E-05,
          7.305867762241E-06, 5.256679519499E-09, -1.016945693300E-07,
          1.518743025531E-08, -1.770252084837E-09, 8.733801272834E-11,
      },
      {
          -2.212823709798E-06, -4.249704742353E-07, 2.831739755462E-06,
          -6.000115718138E-07, 7.597020291557E-10, 3.517154874443E-09,
          4.149084521319E-10, -5.289806153651E-11, 7.196798841269E-13,
      },
      {
          9.648139704737E-08, 4.164960852522E-08, -1.164969298033E-07,
          2.045211951761E-08, 1.799505288362E-09, -4.453195673947E-10,
          -6.803200444549E-12, 3.864394776250E-12, -1.441033650378E-13,
      },
      {
          -1.611904413846E-09, -9.893423877739E-10, 1.785440278790E-09,
          -1.790312871690E-10, -9.280890205774E-11, 2.002478264932E-11,
          -1.151855939531E-12, -8.694978774411E-15, 1.734769090475E-15,
      },
  };

  // RDNE = n;
  // RTE = T;

  // DNE = RDNE;
  // TT = RTE;
  // E = DNE * 1.0E-14;

  // RN = log(E);
  // RT = log(TT);

  // suma = 0.0;
  // for (i = 1; i <= 9; i++) {
    // i1 = i - 1;
    // for (j = 1; j <= 9; j++) {
      // j1 = j - 1;
      // RNJ = pow(RN, j1);
      // if ((RN == 0.0) && (j1 == 0))
        // RNJ = 1.0;
      // RTI = pow(RT, i1);
      // if ((RT == 0.0) && (i1 == 0))
        // RTI = 1.0;
      // suma = suma + MATA[j - 1][i - 1] * RNJ * RTI;
    // }
  // }

  // fEXC = exp(suma)*1.0E-6;
  
  E = n * 1.0E-14;
  E = log(E);
  T = log(T);
  
  suma = 0.0;
  
  for (i = 0; i <=8; i++) {
	  for (j = 0; j <= 8; j++) {
		  suma = suma + MATA[j][i] * pow(E,j) * pow(T,i);
	  }
  }

  fEXC = exp(suma)*1.0E-6;
  return fEXC;
}

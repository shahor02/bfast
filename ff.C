#include <math.h>
#include <fstream>
#include <sstream>


const int kNSolRRanges = 5;
const int kNSolZRanges = 2;
const int kNQuadrants = 4;
const float kSolR2Max0 = 80.*80.;
const float kSolR2Max1 = 250.*250.;
const float kSolR2Max2 = 400.*400.;
const float kSolR2Max3 = 423.*423.;
const float kSolR2Max4 = 500.*500.;

const float kSolZMax = 260.;

int GetQuadrant(float x,float y) {
  return y>0 ? (x>0 ? 0:1) : (x>0 ? 3:2);
}

struct SolParam {
  float mParBxyz[3][20];
};

typedef SolParam SolParam;

SolParam solPar[kNSolRRanges][kNSolZRanges][kNQuadrants];


void ff(double xyz[3])
{
  float x(xyz[0]),y(xyz[1]),z(xyz[2]);
  int izSeg = -1;
  if (z<kSolZMax) {
    if (izSeg>-kSolZMax) izSeg = z<0.f ? 0:1; // solenoid params
    else { // need to check dipole params

    }
  }
  
  float xx = x*x, yy = y*y, rr = xx+yy;
  
  if (rr < kSolR2Max0) {

  }
  else if (rr < kSolR2Max1) {

  }
  else if (rr < kSolR2Max2) {

  }
  else if (rr < kSolR2Max3) {

  }
  else if (rr < kSolR2Max4) {

  }
  else { // redefine to max radius
    float scl2Bond = sqrtf(kSolR2Max4/rr);
    x *= scl2Bond;
    y *= scl2Bond;
  }
  
}

double pl(const double* cf, double x,double y, double z)
{

  /** calculate polynomial
   *   cf[0] + cf[1]*x + cf[2]*y + cf[3]*z + cf[4]*xx + cf[5]*xy + cf[6]*xz + cf[7]*yy + cf[8]*yz + cf[9]*zz +
   *   cf[10]*xxx + cf[11]*xxy + cf[12]*xxz + cf[13]*xyy + cf[14]*xyz + cf[15]*xzz + cf[16]*yyy + cf[17]*yyz + cf[18]*yzz + cf[19]*zzz
  **/

    double val = cf[0] +
    x*(cf[1] + x*(cf[4] + x*cf[10] + y*cf[11] + z*cf[12]) + y*(cf[5]+z*cf[14]) ) +
    y*(cf[2] + y*(cf[7] + x*cf[13] + y*cf[16] + z*cf[17]) + z*(cf[8]) ) +
    z*(cf[3] + z*(cf[9] + x*cf[15] + y*cf[18] + z*cf[19]) + x*(cf[6]) );
    
  return val;
}



bool loadData(const char* inpFName)
{
  std::ifstream in(inpFName,ios::in);
  std::string line;
  int valI,component = -1, nParams = 0, header[4] = {-1,-1,-1,-1}; // iR, iZ, iQuadrant, nVal
  SolParam* curParam = 0; //std::nullptr;
  
  while (std::getline(in, line)) {
    if (line.empty() || line[0]=='#') continue; // empy or comment
    std::stringstream ss(line);
    int cnt=0;
    
    if (component<0) {
      while (cnt<4 && (ss>>header[cnt++]));
      if (cnt!=4) {
	printf("Wrong header: %s\n",line.c_str());
	return false;
      }
      curParam = &solPar[header[0]][header[1]][header[2]];
    }
    else {
      while (cnt<header[3] && (ss>>curParam->mParBxyz[component][cnt++]));
      if (cnt!=header[3]) {
	printf("Wrong data (npar=%d) for param %d %d %d %d: %s\n",cnt,header[0],header[1],header[2],header[3],line.c_str());
	return false;	
      }
    }    
    component++;
    if (component>2) {
      component = -1; // next header expected
      nParams++;
    }
  }
  //
  printf("loaded %d params\n",nParams);
}

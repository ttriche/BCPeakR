#include	<iomanip>
#include	<fstream>
#include	<cmath>
#include	<iostream>
#include	"./TNT/tnt.h"
#include	"./MyFun/myFun.h"

using namespace std;
using namespace TNT;


extern "C" {

class cppoisson {

  public:
    Vector<double>      estPara;
    //cppoisson::cppoisson(Matrix<double>,double,double,double,int,int);
    //cppoisson::~cppoisson(){};
    cppoisson(Matrix<double>,double,double,double,int,int);
    ~cppoisson(){};

    void	BcmixSmooth();

  private:
    int			N, K, M, emIterMax;
    double		p, alpha, beta, invbeta, log00, sumvec;
    Matrix<double>      obs;
    Matrix<int>         fFilter, bFilter;
    Matrix<double>      falpha, fIbeta, fQ, fLog, balpha, bIbeta, bQ, bLog;

    void    SubFun1(double&, double&, double&, double, double, double, double); 
    void    BcmixForward();
    void    BcmixBackward();

}; /* end of -- class cppoisson */


/* hyperprior: Gamma(alpha,beta), corresponding to my SS paper:
	alpha --- gamma, 	beta --- lam
*/
cppoisson::cppoisson(Matrix<double> obs1, double p1,
	double alpha1, double beta1, int K1, int M1)
{
  K=K1;	M=M1;  N=obs1.num_rows()-1;	p=p1;  alpha=alpha1;	beta=beta1;
  invbeta=1.0/beta;	log00=LogGamma(alpha)-alpha*log(invbeta);
  obs.newsize(N+1,4);	obs=obs1;
  fFilter.newsize(N+1, K+1);	bFilter.newsize(N+1, K+1);
  falpha.newsize(N+1, K+1);	balpha.newsize(N+1, K+1);
  fIbeta.newsize(N+1, K+1);	bIbeta.newsize(N+1, K+1);
  fQ.newsize(N+1, K+1);		bQ.newsize(N+1, K+1);
  fLog.newsize(N+1, K+1);	bLog.newsize(N+1, K+1);
  estPara.newsize(N+1);
  BcmixForward();               BcmixBackward();
}


void cppoisson::BcmixSmooth()
{
  double	total, logij, tmpmax;
  Vector<double>	tmpA(K+1);
  Matrix<double>	tmpAlpha(K+1,K+1),tmpIBeta(K+1,K+1),tmpB(K+1,K+1);
  estPara=0.0;
  for (int t=1; t<N; t++) {	total=0.0;	tmpmax=-1.0e+300;
    for (int i=1;i<=min(t,K);i++) {	tmpA[i]=log(p)+log(fQ[t][i]);
      if (tmpmax<tmpA[i])	tmpmax=tmpA[i];
    }
    for (int i=1;i<=min(t,K);i++) for (int j=1;j<=min(N-t,K);j++) {
      tmpAlpha[i][j]=falpha[t][i]+balpha[t+1][j]-alpha;
      tmpIBeta[i][j]=fIbeta[t][i]+bIbeta[t+1][j]-invbeta;
      logij=LogGamma(tmpAlpha[i][j])-tmpAlpha[i][j]*log(tmpIBeta[i][j]);
      tmpB[i][j]= log(1-p) + log(fQ[t][i]) + log(bQ[t+1][j]) +
	logij + log00 -fLog[t][i] - bLog[t+1][j];
      if (tmpmax<tmpB[i][j])	tmpmax=tmpB[i][j];
    }
    for (int i=1; i<=min(t,K); i++) {
      total += tmpA[i] = exp(tmpA[i]-tmpmax);
      estPara[t] += tmpA[i]*falpha[t][i]/fIbeta[t][i];
    }
    for (int i=1;i<=min(t,K);i++) for (int j=1;j<=min(N-t,K);j++) {
      total += tmpB[i][j] = exp(tmpB[i][j]-tmpmax);
      estPara[t] += tmpB[i][j]*tmpAlpha[i][j]/tmpIBeta[i][j];
    }
    estPara[t] /= total;
  }
  for (int i=1;i<=K;i++)
    estPara[N]+=fQ[N][i]*falpha[N][i]/fIbeta[N][i];
}


/****************** Beginning of prIbetaate functions **********************/

void cppoisson::BcmixBackward()
{
  int		s;
  Vector<double>        bQstar(K+2), numTmp(K+2);
  Vector<int>           tmpBackFil(K+2);
  Matrix<double>        tmpPost(K+2,4);
  bQ[N][1]=1.0;		bFilter[N][1] = N;
  SubFun1(balpha[N][1],bIbeta[N][1],bLog[N][1],alpha,invbeta,obs[N][2],obs[N][3]);
  for (int t=N-1; t>=N-K+1; t--) {	s=N+1-t;
    SubFun1(tmpPost[s][1],tmpPost[s][2],tmpPost[s][3],alpha,invbeta,obs[t][2],obs[t][3]);
    numTmp[s] = tmpPost[s][3]-log00;
    for (int i=1; i<s; i++) {
      SubFun1(tmpPost[i][1],tmpPost[i][2],tmpPost[i][3],
        balpha[t+1][i], bIbeta[t+1][i], obs[t][2],obs[t][3]);
      numTmp[i] = tmpPost[i][3]- bLog[t+1][i];
    }
    double      numTmpMax=numTmp[1];
    for (int i=2;i<=s;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1;i<=s;i++) numTmp[i] -= numTmpMax;
    sumvec=bQstar[s]=p*exp(numTmp[s]);
    for (int i=1;i<s;i++) sumvec+=bQstar[i]=(1-p)*bQ[t+1][i]*exp(numTmp[i]);
    for (int i=1; i<=s; i++)  {
      balpha[t][i] = tmpPost[i][1];	bIbeta[t][i] = tmpPost[i][2];
      bLog[t][i] = tmpPost[i][3];     bQ[t][i]=bQstar[i]/sumvec;
      bFilter[t][i] = N+1-i;
    }
  }
  for (int t=N-K; t>=1; t--) {		s=N+1-t;
    SubFun1(tmpPost[K+1][1],tmpPost[K+1][2],tmpPost[K+1][3],
      alpha,invbeta, obs[t][2],obs[t][3]);
    numTmp[K+1]=tmpPost[K+1][3]-log00;		tmpBackFil[K+1]=t;
    for (int i=1; i<K+1; i++) {
      SubFun1(tmpPost[i][1], tmpPost[i][2], tmpPost[i][3],
        balpha[t+1][i], bIbeta[t+1][i], obs[t][2],obs[t][3]);
      numTmp[i] = tmpPost[i][3]- bLog[t+1][i];
      tmpBackFil[i] = bFilter[t+1][i];
    }
    double      numTmpMax=numTmp[1];
    for (int i=1;i<=K+1;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1;i<=K+1;i++) numTmp[i] -= numTmpMax;
    sumvec=bQstar[K+1]=p*exp(numTmp[K+1]);
    for (int i=1;i<K+1;i++) sumvec+=bQstar[i]=(1-p)*bQ[t+1][i]*exp(numTmp[i]);
    int		minind = 1;
    for (int i=1;i<=K+1-M;i++) if (bQstar[minind]>bQstar[i]) minind=i;
    double delta=1.0/(sumvec-bQstar[minind]);
    for (int i=1; i<minind; i++) {
      balpha[t][i]=tmpPost[i][1];	bIbeta[t][i]=tmpPost[i][2];
      bLog[t][i]=tmpPost[i][3];       bFilter[t][i]=tmpBackFil[i];
      bQ[t][i]=bQstar[i]*delta;
    }
    for (int i=minind+1; i<=K+1; i++) {
      balpha[t][i-1]=tmpPost[i][1];	bIbeta[t][i-1]=tmpPost[i][2];
      bLog[t][i-1]=tmpPost[i][3];   bFilter[t][i-1]=tmpBackFil[i];
      bQ[t][i-1]=bQstar[i]*delta;
    }
  }
}

void cppoisson::BcmixForward()
{
  Vector<double>	fQstar(K+2), numTmp(K+2);
  Vector<int>		tmpForFil(K+2);
  Matrix<double>	tmpPost(K+2,4);
  fQ[1][1]=1.0;		fFilter[1][1]=1;
  SubFun1(falpha[1][1],fIbeta[1][1],fLog[1][1],alpha,invbeta,obs[1][2],obs[1][3]);
  for (int t=2; t<=K; t++) {
    SubFun1(tmpPost[t][1],tmpPost[t][2],tmpPost[t][3],alpha,invbeta,obs[t][2],obs[t][3]);
    numTmp[t]=tmpPost[t][3]-log00;
    for (int i=1; i<t; i++) {
      SubFun1(tmpPost[i][1],tmpPost[i][2],tmpPost[i][3],falpha[t-1][i],
	fIbeta[t-1][i],obs[t][2],obs[t][3]);
      numTmp[i] = tmpPost[i][3]-fLog[t-1][i];
    }
    double      numTmpMax=numTmp[1];
    for (int i=2;i<=t;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1;i<=t;i++) numTmp[i] -= numTmpMax;
    sumvec=fQstar[t]=p*exp(numTmp[t]);
    for (int i=1;i<t;i++) sumvec+=fQstar[i]=(1-p)*fQ[t-1][i]*exp(numTmp[i]);
    for (int i=1; i<=t; i++)  {
      falpha[t][i] = tmpPost[i][1];	fIbeta[t][i] = tmpPost[i][2];
      fLog[t][i] = tmpPost[i][3];	fQ[t][i] = fQstar[i]/sumvec;
      fFilter[t][i] = i;
    }
  }

  for (int t=K+1; t<N+1; t++) {
    SubFun1(tmpPost[K+1][1],tmpPost[K+1][2],tmpPost[K+1][3],
      alpha, invbeta,obs[t][2],obs[t][3]);
    numTmp[K+1]=tmpPost[K+1][3]-log00;		tmpForFil[K+1]=t;
    for (int i=1; i<K+1; i++) {
      SubFun1(tmpPost[i][1],tmpPost[i][2],tmpPost[i][3],falpha[t-1][i],
	fIbeta[t-1][i], obs[t][2],obs[t][3]);
      numTmp[i]=tmpPost[i][3]-fLog[t-1][i];
      tmpForFil[i] = fFilter[t-1][i];
    }
    double	numTmpMax=numTmp[1];
    for (int i=2;i<=K+1;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1; i<=K+1; i++)  numTmp[i] -= numTmpMax;
    sumvec=fQstar[K+1]=p*exp(numTmp[K+1]);
    for (int i=1;i<K+1;i++) sumvec+=fQstar[i]=(1-p)*fQ[t-1][i]*exp(numTmp[i]);
    int		minind = 1;
    for (int i=1;i<=K+1-M;i++) if (fQstar[minind]>fQstar[i]) minind=i;
    double delta=1.0/(sumvec-fQstar[minind]);
    for (int i=1; i<minind; i++) {
      falpha[t][i]=tmpPost[i][1];	fIbeta[t][i]=tmpPost[i][2];
      fLog[t][i]=tmpPost[i][3];   fFilter[t][i]=tmpForFil[i];
      fQ[t][i]=fQstar[i]*delta;
    }
    for (int i=minind+1; i<=K+1; i++) {
      falpha[t][i-1]=tmpPost[i][1];	fIbeta[t][i-1]=tmpPost[i][2];
      fLog[t][i-1]=tmpPost[i][3];	fFilter[t][i-1]=tmpForFil[i];
      fQ[t][i-1]=fQstar[i]*delta;
    }
  }
}

void    cppoisson::SubFun1(double &newalpha, double &newinvbeta,
		double &newlog, double oldalpha, double oldinvbeta, 
		double count,double length)
{
  newalpha = oldalpha + (double)(length*count);
  newinvbeta = oldinvbeta + (double)(length);
  newlog = LogGamma(newalpha)-newalpha*log(newinvbeta);
//  cout << "---" << newalpha << "-" << newinvbeta << "-" << LogGamma(newalpha)
//       << "-" << log((double)newinvbeta) << "-" << log((double)65) << endl;
//  printf("%lf %lf\n", newinvbeta, log(newinvbeta));
}

/***************** End of prIbetaate functions ***************************/


}

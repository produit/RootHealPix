#ifndef T_HEAL_FFT_H
#define T_HEAL_FFT_H

#include <cmath>
#include <complex>
#include <vector>

#include "TMath.h"

namespace THealFFT {
  /*
// from bluestein.h of HEALPix C++
int prime_factor_sum (int n);
void bluestein_i (int n, double **tstorage);
void bluestein (int n, double *data, double *tstorage, int isign);

// from fftpack.h of HEALPix C++
void cfftf(int N, double complex_data[], double wrk[]);
void cfftb(int N, double complex_data[], double wrk[]);
void cffti(int N, double wrk[]);
void rfftf(int N, double data[], double wrk[]);             
void rfftb(int N, double data[], double wrk[]);
void rffti(int N, double wrk[]);
  */
// from ls_fft.h of HEALPix C++
typedef struct
  {
  double *work;
  int length;
  int bluestein;
  } complex_plan_i;

typedef complex_plan_i * complex_plan;
/*
complex_plan make_complex_plan (int length);
void kill_complex_plan (complex_plan plan);
void complex_plan_forward (complex_plan plan, double *data);
void complex_plan_backward (complex_plan plan, double *data);
*/
typedef struct
  {
  double *work;
  int length;
  int bluestein;
  } real_plan_i;

typedef real_plan_i * real_plan;
/*
real_plan make_real_plan (int length);
void kill_real_plan (real_plan plan);
void real_plan_forward_fftpack (real_plan plan, double *data);
void real_plan_backward_fftpack (real_plan plan, double *data);
void real_plan_forward_fftw (real_plan plan, double *data);
void real_plan_backward_fftw (real_plan plan, double *data);
void real_plan_forward_c (real_plan plan, double *data);
void real_plan_backward_c (real_plan plan, double *data);
  */

// from fftpack_support.h of HEALPix C++
// slightly modified.

void complex_plan_forward (complex_plan plan, double *data);
void complex_plan_backward (complex_plan plan, double *data);
void real_plan_forward_fftpack (real_plan plan, double *data);
void real_plan_backward_fftpack (real_plan plan, double *data);
void real_plan_forward_c (real_plan plan, double *data);
void real_plan_backward_c (real_plan plan, double *data);

class cfft
  {
  private:
    int n;
    complex_plan plan;

  public:
    cfft () : n(-1), plan(0) {}
    cfft (int size_);
    ~cfft ();
    void Set (int size_);

    int size() const
      { return n; }
    void forward (double *data)
      { complex_plan_forward(plan,data); }
    void backward (double *data)
      { complex_plan_backward(plan,data); }
    void forward (std::vector<std::complex<double> >&data)
      { forward(&(data[0].real())); }
    void backward (std::vector<std::complex<double> >&data)
      { backward(&(data[0].real())); }
  };

class rfft
  {
  private:
    int n;
    real_plan plan;

  public:
    rfft () : n(-1), plan(0) {}
    rfft (int size_);
    ~rfft ();
    void Set (int size_);

    int size() const
      { return n; }

    void forward_fftpack (double *data)
      { real_plan_forward_fftpack(plan,data); }
    void backward_fftpack (double *data)
      { real_plan_backward_fftpack(plan,data); }
    void forward_fftpack (std::vector<double> &data)
      { forward_fftpack(&(data[0])); }
    void backward_fftpack (std::vector<double> &data)
      { backward_fftpack(&(data[0])); }
    void forward_c (std::vector<std::complex<double> >&data)
      { real_plan_forward_c(plan,&(data[0].real())); }
    void backward_c (std::vector<std::complex<double> >&data)
      { real_plan_backward_c(plan,&(data[0].real())); }
  };

// from ylmgen.h of HEALPix C++
class Ylmgen
  {
  private:
    double fsmall, fbig, eps, cth_crit;
    int lmax, mmax, m_last, m_crit;
    std::vector<double> cf;
    std::vector<double> recfac;
    std::vector<double> mfac;

    enum { large_exponent2 = 90, minscale=-4 };

    void recalc_recfac (int m)
      {
      using namespace std;

      if (m_last==m) return;

      int m2 = m*m;
      double f_old=1;
      for (int l=m; l<recfac.size()/2; ++l)
        {
        int l2 = (l+1)*(l+1);
        recfac[2*l]   = sqrt(double(4*l2 - 1) / (l2-m2));
        recfac[2*l+1] = recfac[2*l]/f_old;
        f_old = recfac[2*l];
        }

      m_last=m;
      }

  public:
    /*! Creates a generator which will calculate Y_lm(theta,phi=0)
        up to \a l=l_max and \a m=m_max. It may regard Y_lm whose absolute
        magnitude is smaller than \a epsilon as zero. */
    Ylmgen (int l_max, int m_max, double epsilon=1e-30)
      : eps(epsilon), cth_crit(2.), lmax(l_max), mmax(m_max), m_last(-1),
      m_crit(mmax+1), cf(-minscale+11), recfac((lmax+1)*2), mfac(mmax+1)
      {
      using namespace std;

      fsmall = ldexp(1.,-large_exponent2);
      fbig   = ldexp(1., large_exponent2);
      for (int m=0; m<cf.size(); ++m)
        cf[m] = ldexp(1.,(m+minscale)*large_exponent2);

      mfac[0] = 1;
      for (int m=1; m<mfac.size(); ++m)
        mfac[m] = mfac[m-1]*sqrt((2*m+1.)/(2*m));
      for (int m=0; m<mfac.size(); ++m)
        mfac[m] = (1/TMath::Log(2))*log((1/TMath::Sqrt(4*TMath::Pi()))*mfac[m]);
      }

    /*! For a colatitude given by \a cth and \a sth (representing cos(theta)
        and sin(theta)) and a multipole moment \a m, calculate the
        Y_lm(theta,phi=0) for \a m<=l<=lmax and return in it \a result[l].
        On exit, \a firstl is the \a l index of the first Y_lm with an
        absolute magnitude larger than \a epsilon. If \a firstl>lmax, all
        absolute values are smaller than \a epsilon.
        \a result[l] is undefined for all \a l<firstl. */
    void get_Ylm (double cth, double sth, int m, std::vector<double> &result,
      int &firstl)
      {
      using namespace std;

      if(m <= mmax){
	//planck_assert (m<=mmax, "get_Ylm: m larger than mmax");
	return;
      } // if

      if (((m>=m_crit)&&(abs(cth)>=cth_crit)) || ((m>0)&&(sth==0)))
        { firstl=lmax+1; return; }

      recalc_recfac(m);
      result.resize(lmax+1);

      double logval = mfac[m];
      if (m>0) logval += m*(1/TMath::Log(2))*log(sth);
      int scale = int (logval/large_exponent2)-minscale;
      double corfac = (scale<0) ? 0. : cf[scale];
      double lam_1 = 0;
      double lam_2 = exp(TMath::Log(2)*(logval-(scale+minscale)*large_exponent2));
      if (m&1) lam_2 = -lam_2;

      int l=m;
      while (true)
        {
        if (abs(lam_2*corfac)>eps) break;
        if (++l>lmax) break;
        double lam_0 = cth*lam_2*recfac[2*(l-1)] - lam_1*recfac[2*(l-1)+1];
        if (abs(lam_0*corfac)>eps) { lam_1=lam_2; lam_2=lam_0; break; }
        if (++l>lmax) break;
        lam_1 = cth*lam_0*recfac[2*(l-1)] - lam_2*recfac[2*(l-1)+1];
        if (abs(lam_1*corfac)>eps) { lam_2=lam_1; lam_1=lam_0; break; }
        if (++l>lmax) break;
        lam_2 = cth*lam_1*recfac[2*(l-1)] - lam_0*recfac[2*(l-1)+1];

        while (abs(lam_2)>fbig)
          {
          lam_1 *= fsmall;
          lam_2 *= fsmall;
          ++scale;
          corfac = (scale<0) ? 0. : cf[scale];
          }
        }

      firstl=l;
      if (l>lmax)
        { m_crit=m; cth_crit=abs(cth); return; }

      lam_1*=corfac;
      lam_2*=corfac;

      for (;l<lmax-2;l+=3)
        {
        result[l]=lam_2;
        double lam_0 = cth*lam_2*recfac[2*l] - lam_1*recfac[2*l+1];
        result[l+1] = lam_0;
        lam_1 = cth*lam_0*recfac[2*(l+1)] - lam_2*recfac[2*(l+1)+1];
        result[l+2] = lam_1;
        lam_2 = cth*lam_1*recfac[2*(l+2)] - lam_0*recfac[2*(l+2)+1];
        }
      while (true)
        {
        result[l]=lam_2;
        if (++l>lmax) break;
        double lam_0 = cth*lam_2*recfac[2*(l-1)] - lam_1*recfac[2*(l-1)+1];
        result[l] = lam_0;
        if (++l>lmax) break;
        lam_1 = cth*lam_0*recfac[2*(l-1)] - lam_2*recfac[2*(l-1)+1];
        result[l] = lam_1;
        if (++l>lmax) break;
        lam_2 = cth*lam_1*recfac[2*(l-1)] - lam_0*recfac[2*(l-1)+1];
        }
      }

  };
}

#endif // T_HEAL_FFT_H

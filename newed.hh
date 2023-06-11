//
// newed.hh
// Developer : Jingyi Zhou
// Based on ##
// History:
//   Feb 2023, 
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NEWEP_h
#define NEWEP_h 1

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

#include "Math/Interpolator.h"

#include <gsl/gsl_integration.h>

#include <iostream>
#include <cmath>

#define InterpolPoints 1001
#define InterpolType ROOT::Math::Interpolation::kCSPLINE

#define Abs   TMath::Abs
#define Exp   TMath::Exp
#define Log   TMath::Log
#define DiLog TMath::DiLog
#define Sqrt  TMath::Sqrt
#define Sin   TMath::Sin
#define Cos   TMath::Cos
#define Tan   TMath::Tan
#define ASin  TMath::ASin
#define ACos  TMath::ACos
#define ATan  TMath::ATan
#define ATan2 TMath::ATan2

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const double pi = 3.14159265358979323846;
const double pi2 = pi * pi;
const double deg = pi / 180.0;
const double m = 0.5109989461e-3; // GeV
const double m2 = m * m;
const double m4 = m2 * m2;
const double Md = 1875.612928e-3; 
const double M = 1875.612928e-3; 
//const double M = 938.272046e-3; // GeV
const double M2 = M * M;
const double M4 = M2 * M2;
const double mmu = 105.6583745e-3; // GeV
const double mtau = 1776.82e-3; // GeV
const double alpha = 0.72973525664e-2;
const double alpha_pi = alpha / pi;
const double alpha2 = alpha * alpha;
const double alpha3 = alpha2 * alpha;
const double mu = 2.792782;
const double e = Sqrt(4.0 * pi *alpha);
const double mkb = 389.379; // GeV^{-2} to mkbarn conversion
const double GeVtofm = 25.68189504;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ed.dat";
const char *ifilename = "elastic_ed.info";

int Use_TPE;

double Ei_e;

double theta_min, theta_max;

double Ef_e, theta_e, phi_e;
double Ef_d, theta_d, phi_d;

TLorentzVector vi_e, vi_d;
TLorentzVector vf_e, vf_d;

double omega;

double E_g_min, E_g_cut;
double v_min, v_cut;

double E_g, theta_g, phi_g;
TLorentzVector v_g;

double theta[InterpolPoints];
double xs_elastic_sin[InterpolPoints];
double xs_born_sin[InterpolPoints];

double mode;

TRandom *PseRan;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS(double s, double q2);
double NonRadXS(double s, double q2);

//double BornXS_Sin(double theta);
//double ElasticXS_Sin(double theta);
double BornXS_Sin_ed(double theta);
double ElasticXS_Sin_ed(double theta);

ROOT::Math::Interpolator Interpolator_ElasticXS_Sin_ed(InterpolPoints, InterpolType);
//ROOT::Math::Interpolator Interpolator_BornXS_Sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator Interpolator_BornXS_Sin_ed(InterpolPoints, InterpolType);
ROOT::Math::Interpolator TPE_Feshbach(32, InterpolType);
ROOT::Math::Interpolator TPE_Oleksandr(32, InterpolType);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline double Pow2(double arg) // arg^2
{
    return TMath::Power(arg, 2);
}

inline double Pow3(double arg) // arg^3
{
    return TMath::Power(arg, 3);
}

inline double Pow4(double arg) // arg^4
{
    return TMath::Power(arg, 4);
}

inline double Pow5(double arg) // arg^5
{
    return TMath::Power(arg, 5);
}

inline double Pow6(double arg) // arg^6
{
    return TMath::Power(arg, 6);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GCGMGQ_Model1(double Q2, double &gc, double &gm, double &gq)
{//abbott_2000_1_Model
    double q2 = Q2*GeVtofm;
    double gm0= 1.714;
    double qm0 = 7.37;
    double gq0 = 25.83;
    double qq0 = 8.1;
    double gc0 = 1.;
    double qc0 = 4.21;
    double ac[6] = {0, 6.740e-1, 2.246e-2, 9.806e-3, -2.709e-4, 3.793e-6};
    double am[6] = {0, 5.804e-1, 8.701e-2, -3.624e-3, 3.448e-4, -2.818e-6};
    double aq[6] = {0, 8.796e-1, -5.656e-2, 1.933e-2, -6.734e-4, 9.438e-6};
    double  sum1 = 0.;
    double  sum2 = 0.;
    double  sum3 = 0.;
    for (int n=1; n<=5; n++){
        sum1 = sum1 + ac[n]*pow(q2,n);
    }
    for (int m=1; m<=5; m++){
        sum2 = sum2 + am[m]*pow(q2,m);
    }
    for (int l=1; l<=5; l++){
        sum3 = sum3 + aq[l]*pow(q2,l);
    }
    gc = gc0 * (1. - q2 / (qc0*qc0)) / (1. + sum1);
    gm = gm0 * (1. - q2 / (qm0*qm0)) / (1. + sum2);
    gq = gq0 * (1. - q2 / (qq0*qq0)) / (1. + sum3);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GCGMGQ_Model2(double Q2, double &gc, double &gm, double &gq)
{//abbott_2000_2_Model
    double gev_to_inv_fm = 5.067730758534282;
    double q2=Q2*GeVtofm;
    double eta = Q2 / (4. * pow(Md*gev_to_inv_fm,2));
    double delta = pow((0.89852 * gev_to_inv_fm),2);
    double gq2 = 1. / pow((1. + q2 / (4. * delta)),2);
    double a[4] = {1.57057, 12.23792, -42.04576, 27.92014};
    double alpha[4] = {1.52501, 8.75139, 15.97777, 23.20415};
    double b[4] = {0.07043, 0.14443, -0.27343, 0.05856};
    double beta[4] = {43.67795, 30.05435, 16.43075, 2.80716};
    double c[4] = {-0.16577, 0.27557, -0.05382, -0.05598};
    double gamma[4] = {1.87055, 14.95683, 28.04312, 41.12940};
    double g0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int j=0; j<=3; j++){
       g0 = g0 + (a[j]/(alpha[j]+q2));
    }
    for (int k=0; k<=3; k++){
       sum1  = sum1 + (b[k] / (beta[k] + q2));
    }
    Double_t g1 =sqrt(q2)*sum1;
    for (int l=0; l<=3; l++){
       sum2  = sum2 + (c[l] / (gamma[l] + q2));
    }
    double g2 =q2*sum2;
    double C = gq2*gq2 / (2. * eta + 1.);
    double sqrt_2_eta = sqrt(2 * eta);
    gc = C * ((1. - 2. / 3. * eta) * g0 + 8. / 3. * sqrt_2_eta * g1 + 2. / 3. * (2. * eta - 1.) * g2);
    gm = C * (2. * g0 + 2. * (2. * eta - 1.) / sqrt_2_eta * g1 - 2. * g2);
    gq = C * (-1. * g0 + 2. / sqrt_2_eta * g1 - (1. + 1. / eta) * g2);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GCGMGQ_Model3(double qq, double &gc, double &gm, double &gq)
{//Parker's Model
    double q2 = fabs(qq);
    double gc0 = 1.;
    double qc0 = 4.21;
    double gm0 = 1.714;
    double qm0 = 7.37;
    double gq0 = 25.83;
    double qq0 = 8.1;
    double ac[5] = {0.02708986,0.54510802,0.02708986,0.02708142,0.02708986};
    double am[5] = {0.03189419,0.03189573,0.4433256,0.03189394,0.0318955};
    double aq[5] = {0.04249515,0.04249517,0.0013469,0.04250929,0.51344307};
    double sum1 = 1.0;
    double sum2 = 1.0;
    double sum3 = 1.0;
    for (int m=0; m<5; m++){
        sum1 = sum1 * (1. + ac[m]*q2);
    }
    for (int m=0; m<5; m++){
        sum2 = sum2 * (1. + am[m]*q2);
    }
    for (int m=0; m<5; m++){
        sum3 = sum3 * (1. + aq[m]*q2);
    }
    gc = gc0 * (1. - fabs(qq) / (qc0*qc0))/sum1;
    gm = gm0 * (1. - fabs(qq) / (qm0*qm0))/sum2;
    gq = gq0 * (1. - fabs(qq) / (qq0*qq0))/sum3;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GCGMGQ_Model4(double qq, double &gc, double &gm, double &gq)
{//JZ's refit of SOG
    double q2 = fabs(qq);
    double x = sqrt(q2);
    double gc0 = 1.0;
    double gm0 = 1.714;
    double gq0 = 25.83;
    double gamma = 0.8/Sqrt(3./2.);
    double gamma2 = gamma*gamma;
    double A1[12] = {0.0716166, 236.436, 1454.11, 59.1315, -1052.67, -17.7894, -8.49156, 995.005, -1323.4, -365.277, 23.5912, 0.276615};
    double A2[12] = {-0.0792563, -17.8047, 1.01881, -22.0732, 3.33656, -147.885, 2.72721, -5.09261, 28.6759, 0.280952, 158.688,-0.792919};
    double A3[12] = {0.358712, 3.22108, 0.445367, -9.95379, -7.35359, -247.766, 14.7679, -302.738, 548.153, 1.60897, -0.26087, 0.51755};
    double R1[12] = {5.19665, 3.75539, 3.863, 1.56157, 1.49323, 2.83786, 2.46828, 1.4896, 3.83144, 3.91101, 2.70553, 0.4};
    double R2[12] = {7.00828, 3.16036, 0.533948, 3.72463, 1.95537, 3.34131, 2.33181, 2.08879, 3.6881, 0.0441331, 3.31351, 0.40};
    double R3[12] = {9.51058, 2.77812, 1.36345, 3.73107, 3.41977, 2.2224, 3.59642, 2.27404, 2.24954, 3.98724, 0.57017,0.40};
    double sumA1 = 0;
    double sumA2 = 0;
    double sumA3 = 0;
    for (int k = 0; k<11; k++){
        sumA1 = sumA1 + A1[k];
    }
    for (int k = 0; k<11; k++){
        sumA2 = sumA2 + A2[k];
    }
    for (int k = 0; k<11; k++){
        sumA3 = sumA3 + A3[k];
    }
    A1[11] = 1. - sumA1;
    A2[11] = 1. - sumA2;
    A3[11] = 1. - sumA3;
    double sum1 = 0.;
    double sum2 = 0.;
    double sum3 = 0.;
    for (int k = 0; k<12; k++){
        sum1 = sum1 + A1[k]/(1.+2.*R1[k]*R1[k]/gamma2)*(cos(x*R1[k])+(2.*R1[k]*R1[k]/gamma2)*sin(x*R1[k])/(x*R1[k]));
    }
    for (int k = 0; k<12; k++){
        sum2 = sum2 + A2[k]/(1.+2.*R2[k]*R2[k]/gamma2)*(cos(x*R2[k])+(2.*R2[k]*R2[k]/gamma2)*sin(x*R2[k])/(x*R2[k]));
    }
    for (int k = 0; k<12; k++){
        sum3 = sum3 + A3[k]/(1.+2.*R3[k]*R3[k]/gamma2)*(cos(x*R3[k])+(2.*R3[k]*R3[k]/gamma2)*sin(x*R3[k])/(x*R3[k]));
    }
    gc = gc0*Exp(-q2*gamma2/4.)*sum1;
    gm = gm0*Exp(-q2*gamma2/4.)*sum2;
    gq = gq0*Exp(-q2*gamma2/4.)*sum3;
}

void GetGCGMGQ(double qq, double &gc, double &gm, double &gq)
{
    if (mode == 0){
        GCGMGQ_Model1(qq, gc, gm, gq);
    }
    if (mode == 1){
        GCGMGQ_Model1(qq, gc, gm, gq);
    }
    if (mode == 2){
        GCGMGQ_Model2(qq, gc, gm, gq);
    }
    if (mode == 3){
        GCGMGQ_Model3(qq, gc, gm, gq);
    }
    if (mode == 4){
        GCGMGQ_Model4(qq, gc, gm, gq);
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double A(double QQ) // Structure function A
{
  double tau = QQ/Pow2(2.*M);
  double gc, gm, gq;
  GetGCGMGQ(QQ, gc, gm, gq);
  return Pow2(gc) + 8.*Pow2(tau*gq)/9. + 2.*tau*Pow2(gm)/3.;//Eq.(8)
}
double B(double QQ) // Structure function B
{
  double tau = QQ/Pow2(2.*M);
  double gc, gm, gq;
  GetGCGMGQ(QQ, gc, gm, gq);
  return 4.*tau*(1. + tau)*Pow2(gm)/3.;//Eq.(8)
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)//GeV
{
    return ((Ei_e + M) * (M * Ei_e + m2) + Sqrt(M2 - Pow2(m * Sin(theta))) * (Pow2(Ei_e) - m2) * Cos(theta)) / (Pow2(Ei_e + M) - (Pow2(Ei_e) - m2) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalSQ2(double theta, double &s, double &q2)
{
    Ef_e = ElasticEnergy(theta);

    vf_e.SetPxPyPzE(Sqrt(Pow2(Ef_e) - m2) * Sin(theta), 0.0, Sqrt(Pow2(Ef_e) - m2) * Cos(theta), Ef_e);
    vf_d = vi_e + vi_d - vf_e;

    s = 2.0 * vi_e * vi_d;
    q2 = -(vi_e - vf_e) * (vi_e - vf_e);//Q2 in the paper, unit: GeV2
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double Theta_B1(double q2)
{
    return q2 - 2.*m2;
}
double Theta_B2(double s, double q2)
{
    double x = s - q2; 

    return (s*x-M2*q2)/(2.*M2);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double BornXS_dQ(double theta)//dsigmaBorn/dQ2
{
    double s, q2; 
    CalSQ2(theta, s, q2);

    double lambda_s = s * s - 4.0 * m2 * M2;

    double W1 = 2.*M2*B(q2);
    double W2 = 4.*M2*A(q2);

    return 2. * pi * alpha2 * ( Theta_B1(q2) * W1 + Theta_B2(s, q2) * W2 )/( lambda_s * Pow2(q2) );
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double BornXS_Sin_ed(double theta)//dsigmaBorn/dtheta,already the integrand
{
    double s, q2;
    CalSQ2(theta, s, q2);
    
    double Ef_e = ElasticEnergy(theta);

    double MottXS = Pow2(2.*alpha*Ef_e/q2)*(Ef_e/Ei_e)*Pow2(Cos(theta/2.)); //Eq.(9)

    double BornXSsin = 2.*pi*MottXS*(A(q2)+B(q2)*Pow2(Tan(theta/2.)))*Sin(theta);//Eq.(7)

    return BornXSsin;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double dsigma_AMM_dqQ2(double theta)//anomalous magnetic moments
{
    double s, q2;
    CalSQ2(theta, s, q2);

    double x = s - q2;

    double W1 = 2.*M2*B(q2);
    double W2 = 4.*M2*A(q2);

    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq.(26)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2)); // eq.(43)
    double lambda_s = s * s - 4.0 * m2 * M2; // eq.(43)
    double sigma_AMM = alpha3 * m2 * l_m * (12.0 * M2 * W1 - (q2 + 4.0 * M2) * W2) / (2.0 * M2 * q2 * lambda_s);//Eq.(53)
    return sigma_AMM;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double delta_vac(double theta)//vacuum polarization of lepton and hardon
{
    double s, q2;
    CalSQ2(theta, s, q2);
    // eq. (48) lepton vacuum polarization
    double delta_vac_l = 0.0;
    double lepton_mass[3] = {m, mmu, mtau};

    for (auto &vac_m : lepton_mass) {
        double vac_m2 = vac_m * vac_m;
        double vac_slambda_m = Sqrt(q2 * q2 + 4.0 * vac_m2 * q2);
        double vac_l_m = 1.0 / vac_slambda_m * Log((vac_slambda_m + q2) / (vac_slambda_m - q2));//Eq.(49)
        delta_vac_l += 2.0 / 3.0 * (q2 + 2.0 * vac_m2) * vac_l_m - 10.0 / 9.0 + 8.0 / 3.0 * vac_m2 / q2 * (1.0 - 2.0 * vac_m2 * vac_l_m);
    }
    // eq. (51) hardon vacuum polarization
    // fit parameters A,B,C from Table 1 of arXiv:2210.03785 [hep-ph]
    double A = -1.345e-9;
    double B = -2.302e-3;
    double C = 4.091;
    double delta_vac_h = -2.*pi*(A+B*Log(1.+C*q2))/alpha;

    //return delta_vac_l;//just for comparision with ep
    return delta_vac_l+delta_vac_h;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double delta_VR(double theta)//sum of infrared divergence terms and leptonic vertex correction
{
    double s, q2;
    CalSQ2(theta, s, q2);
    double x = s - q2;
    double Sp = s + x;//Eq.(35)

    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq. (26)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2)); // eq. (43)

    double lambda_s = s * s - 4.0 * m2 * M2; // eq. (21)
    double slambda_s = Sqrt(lambda_s);
    double l_s = 1.0 / slambda_s * Log((s + slambda_s) / (s - slambda_s)); // eq. (43)

    double lambda_x = x * x - 4.0 * m2 * M2; // eq. (43)
    double slambda_x = Sqrt(lambda_x);
    double l_x = 1.0 / slambda_x * Log((x + slambda_x) / (x - slambda_x)); // eq. (43)

    double rho = (q2+slambda_m)*(Sp-slambda_m)/slambda_m;

    double S_phi = ((q2+2.*m2)/slambda_m)*(0.25*lambda_x*l_x*l_x - 0.25*lambda_s*l_s*l_s + DiLog(1-(x+slambda_x)*rho/(8.*m2*M2)) + DiLog(1-(rho/(2.*(x+slambda_x)))) - DiLog(1-(q2*(s+slambda_s)*rho)/(2.*M2*Pow2(q2+slambda_m))) - DiLog(1-(2.*m2*q2*rho)/(Pow2(q2+slambda_m)*(s+slambda_s))));//Eq.(44)

    double q2_m = q2 + 2.0 * m2;
    double q2mlm = q2_m * l_m - 1.0;

    //double v_limit = 0.99 * (slambda_s*slambda_m-q2*(s+2.*m2))/(2.*m2);//0.99 of the upper limit for v at fixed q2;//Eq.(23)
    //double v_cut_max = (v_limit > v_cut) ? v_cut : v_limit;
    
    return 2.0 * q2mlm * Log(v_min / m / M) + 0.5 * (s * l_s + x * l_x) + S_phi + (1.5 * q2 + 4.0 * m2) * l_m - 2.0 - q2_m / slambda_m * (0.5 * lambda_m * l_m * l_m + 2.0 * DiLog(2.0 * slambda_m / (q2 + slambda_m)) - pi2 / 2.0); // eq. (47)
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double delta_inf(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);
    
    double x = s - q2;

    double q2_m = q2 + 2.0 * m2;
    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq. (26)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2));

    double lambda_s = s * s - 4.0 * m2 * M2; // eq. (21)
    double v_limit = 0.99 * 2.0 * q2 * (lambda_s - q2 * (s + m2 + M2)) / (q2 * (s + 2.0 * m2) + Sqrt(q2 * lambda_s * (q2 + 4.0 * m2)));
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;

    return (q2_m * l_m - 1.) * Log(v_max * v_max / s / x);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double delta_add(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);
    
    double x = s - q2;

    double q2_m = q2 + 2.0 * m2;
    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq. (26)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2));
    double q2mlm = q2_m * l_m - 1.0;

    double lambda_s = s * s - 4.0 * m2 * M2; // eq. (21)
    double v_limit = 0.99 * 2.0 * q2 * (lambda_s - q2 * (s + m2 + M2)) / (q2 * (s + 2.0 * m2) + Sqrt(q2 * lambda_s * (q2 + 4.0 * m2)));
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;

    return -2.0 * alpha_pi * q2mlm * Log(v_max / v_min);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double Jacobian_theta(double s, double q2)
{
    double x = s - q2;
    double Sp = s + x;//Eq.(35)

    double lambda_s = s * s - 4.0 * m2 * M2; // eq. (21)
    double slambda_s = Sqrt(lambda_s);
    double l_s = 1.0 / slambda_s * Log((s + slambda_s) / (s - slambda_s)); // eq. (43)

    double lambda_x = x * x - 4.0 * m2 * M2; // eq. (43)
    double slambda_x = Sqrt(lambda_x);
    double l_x = 1.0 / slambda_x * Log((x + slambda_x) / (x - slambda_x)); // eq. (43)


    return -(slambda_s*Pow3(slambda_x))/(2.*M2*(s*x-2.*m2*(q2+2.*M2)));//eq. (61) 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Finite part of the infrared free bremsstrahlung cross-section integrant of tau
double sig_rad_1_dQ(double s, double q2, double v, double tau)
{
    double x = s - q2 ;//def in Sasha
    double Sp = 2*s - q2;

    double lambda_s = s * s - 4.0 * m2 * M2;
    double lambda_q = Pow2(v + q2) + 4.0 * M2 * q2;//==lambda_y
   
    double C1 = 4. * m2 * ( q2 + tau * (q2+v) - Pow2(tau) * M2)+ Pow2( q2 + tau * s );
    double C2 = 4. * m2 * ( q2 + tau * (q2+v) - Pow2(tau) * M2)+ Pow2( q2 + tau * (v-x) );
    double B1 = tau * ( s * (q2+v) + 2. * M2 * q2) + q2 * (Sp-v);
    double B2 = tau * ( (x-v) * (q2+v) - 2. * M2 * q2) + q2*(Sp-v); 
    double sc1 = Sqrt(C1);
    double sc2 = Sqrt(C2);
 
    double Fd = (1./sc2 - 1./sc1)/tau;
    double F1p = 1./sc1 + 1./sc2;
    double F2p = m2 * (B2/sc2/C2 + B1/sc1/C1);
    double F2m = m2 * (B2/sc2/C2 - B1/sc1/C1);
    double FIR = F2p - (q2 + 2.*m2) * Fd;
    double F = 1.0 / Sqrt(lambda_q);

    double theta_11 = 4.0 * (q2 - 2.0 * m2) * FIR;
    double theta_12 = 4.0 * FIR * tau;
    double theta_13 = - 4.0 * F - 2.0 * Pow2(tau) * Fd;
    double theta_21 = 2. * ( s*x - M2*q2) * FIR/M2;
    double theta_22 = (2. * (q2 - 2.*tau*M2 - 2.*(1.+tau)*s) * FIR + Sp * (q2*F1p + 2.*F2m - tau*Sp*Fd))/(2.*M2);
    double theta_23 = ((tau * (2.*tau*M2 - q2) + 4.*m2) * Fd - Sp*F1p + 2.*(1.+tau)*(tau*Sp*Fd + x*F1p + FIR - F2m))/(2.*M2) + 2.*F;
    double theta_24 = -tau * (1.+tau) * (F1p + (tau+2.)*Fd)/(2.*M2);   

    double R = v/(1.+tau);
    double Q2tilde = q2 + R*tau;
    double gc, gm, gq;
    GetGCGMGQ(Q2tilde, gc, gm, gq);
    double tau_t = Q2tilde / 4.0 / M2;
    double A =  Pow2(gc) + 8. * Pow2(tau_t * gq)/9. + 2. * tau_t * Pow2(gm)/3.;//Eq.(8)
    double B =  4. * tau_t * (1. + tau_t) * Pow2(gm)/3.;//Eq.(8)
    double W1 = 2. * M2 * B;
    double W2 = 4. * M2 * A;

    double factor1 = W1 * ( theta_11 / R + theta_12 + theta_13 * R );
    double factor2 = W2 * ( theta_21 / R + theta_22 + theta_23 * R + theta_24 * Pow2(R) );

    double integrand = ( factor1 + factor2 ) / ((1+tau)*Pow2(Q2tilde));
    double fac = - alpha3 / ( 2. * lambda_s);
    //std::cout<<s<<" "<<q2<<" "<<v<<" "<<tau<<" "<<fac<<" "<<integrand<<" "<<fac*integrand<<std::endl;
    return integrand;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double func1_dQ(double x, void *params)
{
    double *p = (double *)params;

    double F1 = sig_rad_1_dQ(p[0], p[1], p[2], x);
    return F1;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Finite part of the infrared free bremsstrahlung cross-section integrated over tau
double sig_rad_2_dQ2(double s, double q2, double v)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    double result, error;
    double par[3] = {s, q2, v};

    gsl_function F;
    F.function = &func1_dQ;
    F.params = &par;

    double x = s - q2;
    double lambda_q = Pow2(v + q2) + 4.0 * M2 * q2;

    double tau_max = (v + q2 + Sqrt(lambda_q)) / 2.0 / M2;
    double tau_min = (v + q2 - Sqrt(lambda_q)) / 2.0 / M2;


    double pts[5] = {tau_min, 0.9*(tau_min+tau_max)/2., (tau_min+tau_max)/2., 1.1*(tau_min+tau_max)/2.,tau_max};
    //gsl_integration_qagp(&F, pts, 5, 0, 1e-6, 10000, w, &result, &error);//2.2GeV
    gsl_integration_qags(&F, tau_min, tau_max, 0, 1e-6, 10000, w, &result, &error);//1.1GeV

    gsl_integration_workspace_free(w);

    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq. (26)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2)); // eq.

    double J0 = 2. * (( q2 + 2. * m2 ) * l_m - 1.); 

    double W1 = 2.*M2*B(q2);
    double W2 = 4.*M2*A(q2);
 
    double sum1 = 4. * J0 * ( Theta_B1(q2) * W1 + Theta_B2(s,q2) * W2 ) / ( v * Pow2(q2) );

    return sum1 + result;


}

double func2(double x, void *params)
{
    double *p = (double *)params;

    double F2 = sig_rad_2_dQ2(p[0], p[1], x);
    return F2;
}

//double sig_rad_test(double s, double q2)
//{
//    double v, sig_rad2;
//    for (int i=0; i<1000; i++){
//        v = 0.001*i/1000;
//        sig_rad2 = sig_rad_2_dQ2(s,q2,v);
//        //std::cout<<s<<" "<<q2<<" "<<v<<" "<<sig_rad2<<std::endl;
//    }
//    //std::cout<<s<<" "<<q2<<" "<<v<<" "<<sig_rad2<<std::endl;
//    return 0.;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Finite part of the infrared free bremsstrahlung cross-section integrated over tau and v
double sig_rad_3_dQ2(double s, double q2)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

    double result, error;
    double par[2] = {s, q2};

    gsl_function F;
    F.function = &func2;
    F.params = &par;

    gsl_integration_qags(&F, 0, v_min, 0, 1e-6, 10000, w, &result, &error);

    gsl_integration_workspace_free(w);

    double lambda_s = s * s - 4.0 * m2 * M2;
    
    return - alpha3 * result / ( 2. * lambda_s);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Finite part of the infrared free bremsstrahlung cross-section integrant of tau
double sig_rad_1_dtheta(double theta, double s, double q2, double v, double tau)
{
    double x = s - q2 ;//def in Sasha
    double Sp = 2*s - q2;

    double lambda_s = s * s - 4.0 * m2 * M2;
    double lambda_q = Pow2(v + q2) + 4.0 * M2 * q2;//==lambda_y
    double R = v/(1.+tau);
    double Q2tilde = q2 + R*tau;
   
    double D = M2 * (lambda_s + v * (v - 2.*s)) - m2 * (lambda_s * Pow2(Sin(theta)) + 4.* v * M2);
    double Q2_R = ( (s + 2. * M2)*( lambda_s - v * s) - lambda_s * (s - v)*Pow2(Cos(theta)) - 2. * M * Sqrt(lambda_s) * Sqrt(D) * Cos(theta))/(Pow2( s + 2. * M2) - lambda_s * Pow2(Cos(theta)));

    q2 = Q2_R;

    double C1 = 4. * m2 * ( q2 + tau * (q2+v) - Pow2(tau) * M2)+ Pow2( q2 + tau * s );
    double C2 = 4. * m2 * ( q2 + tau * (q2+v) - Pow2(tau) * M2)+ Pow2( q2 + tau * (v-x) );
    double B1 = tau * ( s * (q2+v) + 2. * M2 * q2) + q2 * (Sp-v);
    double B2 = tau * ( (x-v) * (q2+v) - 2. * M2 * q2) + q2*(Sp-v); 
    double sc1 = Sqrt(C1);
    double sc2 = Sqrt(C2);
 
    double Fd = (1./sc2 - 1./sc1)/tau;
    double F1p = 1./sc1 + 1./sc2;
    double F2p = m2 * (B2/sc2/C2 + B1/sc1/C1);
    double F2m = m2 * (B2/sc2/C2 - B1/sc1/C1);
    double FIR = F2p - (q2 + 2.*m2) * Fd;
    double F = 1.0 / Sqrt(lambda_q);

    double theta_11 = 4.0 * (q2 - 2.0 * m2) * FIR;
    double theta_12 = 4.0 * FIR * tau;
    double theta_13 = - 4.0 * F - 2.0 * Pow2(tau) * Fd;
    double theta_21 = 2. * ( s*x - M2*q2) * FIR/M2;
    double theta_22 = (2. * (q2 - 2.*tau*M2 - 2.*(1.+tau)*s) * FIR + Sp * (q2*F1p + 2.*F2m - tau*Sp*Fd))/(2.*M2);
    double theta_23 = ((tau * (2.*tau*M2 - q2) + 4.*m2) * Fd - Sp*F1p + 2.*(1.+tau)*(tau*Sp*Fd + x*F1p + FIR - F2m))/(2.*M2) + 2.*F;
    double theta_24 = -tau * (1.+tau) * (F1p + (tau+2.)*Fd)/(2.*M2);   

    double gc, gm, gq;
    GetGCGMGQ(Q2tilde, gc, gm, gq);
    double tau_t = Q2tilde / 4.0 / M2;
    double A =  Pow2(gc) + 8. * Pow2(tau_t * gq)/9. + 2. * tau_t * Pow2(gm)/3.;//Eq.(8)
    double B =  4. * tau_t * (1. + tau_t) * Pow2(gm)/3.;//Eq.(8)
    double W1 = 2. * M2 * B;
    double W2 = 4. * M2 * A;

    double factor1 = W1 * ( theta_11 / R + theta_12 + theta_13 * R );
    double factor2 = W2 * ( theta_21 / R + theta_22 + theta_23 * R + theta_24 * Pow2(R) );

    double integrand = ( factor1 + factor2 ) / ((1+tau)*Pow2(Q2tilde));

    return integrand;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double func1_dtheta(double x, void *params)
{
    double *p = (double *)params;

    double F1 = sig_rad_1_dtheta(p[0], p[1], p[2], p[3],x);
    return F1;
}

// Finite part of the infrared free bremsstrahlung cross-section integrated over tau
double sig_rad_2_dtheta(double theta, double v)
{
    double s, q2;
    CalSQ2(theta, s, q2);
 

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    double result, error;
    double par[4] = {theta, s, q2, v};

    gsl_function F;
    F.function = &func1_dtheta;
    F.params = &par;

    double x = s - q2;
    double lambda_q = Pow2(v + q2) + 4.0 * M2 * q2;

    double tau_max = (v + q2 + Sqrt(lambda_q)) / 2.0 / M2;
    double tau_min = (v + q2 - Sqrt(lambda_q)) / 2.0 / M2;

    double pts[5] = {tau_min, 0.9*(tau_min+tau_max)/2., (tau_min+tau_max)/2., 1.1*(tau_min+tau_max)/2.,tau_max};
    gsl_integration_qagp(&F, pts, 5, 0, 1e-7, 10000, w, &result, &error);//1.1GeV
    //gsl_integration_qags(&F, tau_min, tau_max, 0, 1e-5, 1000, w, &result, &error);
    //gsl_integration_qags(&F, tau_min, tau_max, 0, 1e-6, 10000, w, &result, &error);

    gsl_integration_workspace_free(w);


    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq. (26)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2)); // eq.

    double J0 = 2. * (( q2 + 2. * m2 ) * l_m - 1.); 

    double W1 = 2.*M2*B(q2);
    double W2 = 4.*M2*A(q2);
 
    double sum1 = 4. * Jacobian_theta(s,q2) * J0 * ( Theta_B1(q2) * W1 + Theta_B2(s,q2) * W2 ) / ( v * Pow2(q2) );

    double lambda_s = s * s - 4.0 * m2 * M2;
    double D = M2 * (lambda_s + v * (v - 2.*s)) - m2 * (lambda_s * Pow2(Sin(theta)) + 4.* v * M2);
    double Q2_R = ( (s + 2. * M2)*( lambda_s - v * s) - lambda_s * (s - v)*Pow2(Cos(theta)) - 2. * M * Sqrt(lambda_s) * Sqrt(D) * Cos(theta))/(Pow2( s + 2. * M2) - lambda_s * Pow2(Cos(theta)));

    double fac1 = (lambda_s - v * s - Q2_R * (s + 2. * M2))/( Pow2( s + 2. * M2) - lambda_s * Pow2(Cos(theta)) );
    double fac2 = ( s + 2. * M2)/(Cos(theta)) + M * Sqrt( lambda_s / D) * (s - v + 2.*m2); 
    double J_theta = - fac1 * fac2;
    
    return sum1 + J_theta * result;
}

double func3(double x, void *params)
{
    double *p = (double *)params;

    double F2 = sig_rad_2_dtheta(p[0], x);
    return F2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Finite part of the infrared free bremsstrahlung cross-section integrated over tau and v
double sig_rad_3_dtheta(double theta)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    double result, error;
    double par[1] = {theta};

    gsl_function F;
    F.function = &func3;
    F.params = &par;

    gsl_integration_qags(&F, 0, v_min, 0, 1e-6, 10000, w, &result, &error);

    gsl_integration_workspace_free(w);

    double s, q2;
    CalSQ2(theta, s, q2);
 
    double lambda_s = s * s - 4.0 * m2 * M2;

    return Sin(theta) * alpha3 * result / ( 2. * lambda_s);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double ElasticXS_Sin_ed(double theta)//dtheta
{
    double s, q2;
    CalSQ2(theta, s, q2);

    //double sig_0 = BornXS_Sin_ed(theta);

    double jacob_sin = - Jacobian_theta(s,q2) * Sin(theta);
    double sig_0 = BornXS_dQ(theta) * jacob_sin;


    double sig_Fs = sig_rad_3_dtheta(theta);
    //double sig_Fs = sig_rad_3_dQ2(s,q2) * jacob_sin;
    double sigma_AMM = dsigma_AMM_dqQ2(theta) * jacob_sin;

    double Efe = ElasticEnergy(theta);
    double result = (1.0 + alpha_pi * (delta_VR(theta) + delta_vac(theta) - delta_inf(theta))) * sig_0 * Exp(alpha_pi * delta_inf(theta)) + sigma_AMM  + sig_Fs;
    //std::cout<<theta/deg<<" "<<result<<" "<<sig_0<<" "<<((1.0 + alpha_pi * (delta_VR(theta) + delta_vac(theta) - delta_inf(theta))) * sig_0 * Exp(alpha_pi * delta_inf(theta)) )/sig_0<<" "<<sigma_AMM/sig_0<<" "<<sig_Fs/sig_0<<std::endl;
    //std::cout<<theta/deg<<" "<<Efe<<" "<<result<<" "<<sig_0<<std::endl;

    return result;
}




class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_e = theta_min + arg[0] * (theta_max - theta_min);

        return (theta_max - theta_min) * Interpolator_ElasticXS_Sin_ed.Eval(theta_e);
    }
};

double matrix_r(double s, double q2, double v, double tau, double phi)
{
    double x = s - q2;
    double Sp = s + x;//Eq.(35)
    double lambda_q = Pow2(v + q2) + 4.0 * M2 * q2;

    double tau_max = (v + q2 + Sqrt(lambda_q)) / 2.0 / M2;
    double tau_min = (v + q2 - Sqrt(lambda_q)) / 2.0 / M2;

    if (tau <= tau_min || tau >= tau_max) return 0.0;

    double lambda_z = (tau - tau_min) * (tau_max - tau) * (s * (s - q2 - v) * q2 - M2 * Pow2(q2) - m2 * lambda_q);
    double slambda_z = (lambda_z > 0.0) ? Sqrt(lambda_z) : 0.0;

    double z1 = (q2 * (s + (s - q2 - v)) + tau * (s * (q2 + v) + 2.0 * M2 * q2) - 2.0 * M * slambda_z * Cos(phi)) / lambda_q;
    double z2 = (q2 * (s + (s - q2 - v)) + tau * ((s - q2 - v) * (q2 + v) - 2.0 * M2 * q2) - 2.0 * M * slambda_z * Cos(phi)) / lambda_q;

    double F = 0.5 / pi / Sqrt(lambda_q);
    double Fd = F / z1 / z2;
    double F1p = F / z1 + F / z2;
    double F2p = F * m2 * (1.0 / Pow2(z2) + 1.0 / Pow2(z1));
    double F2m = F * m2 * (1.0 / Pow2(z2) - 1.0 / Pow2(z1));
    double FIR = F2p - (q2 + 2.0 * m2) * Fd;

    double theta_11 = 4.0 * (q2 - 2.0 * m2) * FIR;
    double theta_12 = 4.0 * FIR * tau;
    double theta_13 = - 4.0 * F - 2.0 * Pow2(tau) * Fd;
    double theta_21 = 2. * ( s*x - M2*q2) * FIR/M2;
    double theta_22 = (2. * (q2 - 2.*tau*M2 - 2.*(1.+tau)*s) * FIR + Sp * (q2*F1p + 2.*F2m - tau*Sp*Fd))/(2.*M2);
    double theta_23 = ((tau * (2.*tau*M2 - q2) + 4.*m2) * Fd - Sp*F1p + 2.*(1.+tau)*(tau*Sp*Fd + x*F1p + FIR - F2m))/(2.*M2) + 2.*F;
    double theta_24 = -tau * (1.+tau) * (F1p + (tau+2.)*Fd)/(2.*M2);   

    //double theta_11 = 4.0 * (q2 - 2.0 * m2) * F_IR;
    //double theta_12 = 4.0 * F_IR * tau;
    //double theta_13 = -4.0 * F - 2.0 * Pow2(tau) * F_d;
    //double theta_21 = 2.0 * (s * x - M2 * q2) * F_IR / M2;
    //double theta_22 = (2.0 * (s + x) * F_2m + (Pow2(s) - Pow2(x)) * F_1p + 2.0 * (s - x - 2.0 * M2 * tau) * F_IR - tau * Pow2(s + x) * F_d) / 2.0 / M2;
    //double theta_23 = (4.0 * M2 * F + (4.0 * m2 + 2.0 * M2 * Pow2(tau) - (s - x) * tau) * F_d - (s + x) * F_1p) / 2.0 / M2;

    double factor = 1.0 + tau;
    double r = v / factor;
    double t = q2 + r * tau;
    double tau_t = t / 4.0 / M2;

    double gc, gm, gq;
    GetGCGMGQ(t, gc, gm, gq);
    double A =  Pow2(gc) + 8. * Pow2(tau_t * gq)/9. + 2. * tau_t * Pow2(gm)/3.;//Eq.(8)
    double B =  4. * tau_t * (1. + tau_t) * Pow2(gm)/3.;//Eq.(8)
    double W1 = 2. * M2 * B;
    double W2 = 4. * M2 * A;

    double theta_1 = theta_11 / r + theta_12 + theta_13 * r;
    double theta_2 = theta_21 / r + theta_22 + theta_23 * r + theta_24 * r * r;

    double result = Pow6(e) / Pow2(t) * (-4.0 * pi) * Sqrt(lambda_q) * (theta_1 * W1 + theta_2 * W2) / r; // matrix, Eq(14),(15)

    return result;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BremsIntegrand: public TFoamIntegrand
{
public:
    BremsIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_e = theta_min + arg[0] * (theta_max - theta_min);
        E_g = E_g_min + arg[1] * (E_g_cut - E_g_min);
        theta_g = arg[2] * pi;
        phi_g = arg[3] * 2.0 * pi;

        if (E_g > M * (Ei_e - m) / (M + Ei_e - Sqrt(Pow2(Ei_e) - m2) * Cos(theta_g))) return 0.0;

        double A = Sqrt(Pow2(Ei_e) - m2) * Cos(theta_e) - E_g * (Cos(theta_e) * Cos(theta_g) + Sin(theta_e) * Sin(theta_g) * Cos(phi_g));
        double B = Ei_e + M - E_g;
        double C = E_g * (Ei_e + M - Sqrt(Pow2(Ei_e) - m2) * Cos(theta_g)) - M * Ei_e - m2;

        if (m2 * (Pow2(A) - Pow2(B)) + Pow2(C) < 0.0) return 0.0;

        if (Abs(Pow2(A) - Pow2(B)) < 1.0e-12) return 0.0;

        Ef_e = (B * C - A * Sqrt(m2 * (Pow2(A) - Pow2(B)) + Pow2(C))) / (Pow2(A) - Pow2(B));

        if (Abs(A * Sqrt(Pow2(Ef_e) - m2) - B * Ef_e - C) > 1.0e-9) return 0.0;

        if (Ef_e < m || Ef_e > Ei_e - E_g) return 0.0;

        v_g.SetPxPyPzE(E_g * Sin(theta_g) * Cos(phi_g), E_g * Sin(theta_g) * Sin(phi_g), E_g * Cos(theta_g), E_g);
        vf_e.SetPxPyPzE(Sqrt(Pow2(Ef_e) - m2) * Sin(theta_e), 0.0, Sqrt(Pow2(Ef_e) - m2) * Cos(theta_e), Ef_e);
        vf_d = vi_e + vi_d - vf_e - v_g;

        TVector3 k1 = vi_e.Vect(), k2 = vf_e.Vect();
        TVector3 q = (vi_e - vf_e).Vect(), k = v_g.Vect();
        TVector3 k1k2 = k1.Cross(k2);
        TVector3 qk = q.Cross(k);

        double s = 2.0 * Ei_e * M;
        double q2 = -(vi_e - vf_e) * (vi_e - vf_e);
        double v = (vi_e + vi_d - vf_e) * (vi_e + vi_d - vf_e) - M2;
        double tau = v_g * (vi_e - vf_e) / (v_g * vi_d);
        double phik = qk.Angle(k1k2);

        double matrix = matrix_r(s, q2, v, tau, phik);

        double xs = matrix * (1.0 / (1024.0 * Pow5(pi) * Sqrt(Pow2(Ei_e) - m2) * M)) * E_g * ((Pow2(Ef_e) - m2) / Abs(A * Ef_e - B * Sqrt(Pow2(Ef_e) - m2)));

        return (theta_max - theta_min) * (E_g_cut - E_g_min) * pi * (2.0 * pi) * (2.0 * pi) * xs * Sin(theta_e) * Sin(theta_g) * mkb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

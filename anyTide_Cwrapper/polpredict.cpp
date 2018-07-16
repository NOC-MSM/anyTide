//
//  cppTest.cpp
//  POL Tide Predictor
//
//  Created by Nick Thorne on 06/12/2012.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

extern "C"
{

    //Test function
    int aFunc(int n)
    {

        return n*3;
    }
        
}


#include <stdio.h>
#include <cstring>
#include <string.h>
#include <iostream>
#include <stdlib.h>     /* strtod */
#include <cmath>  
#include <Python.h>
#define MAXCONST 120


double dmod(double a, double b)
{
    while (a>=b)
        a=a-b;
    while (a<0.0)
        a=a+b;
    return a;
}

int IsLeapYear(int yy)
{
    int x=0;
    if ((yy % 4) == 0) x=1;
    if ((yy % 100)==0) x=0;
    if ((yy % 400)==0) x=1;
    return x;
}

int ddmm2jul(int dd, int mm, int yy)
{
    int ndys[]={0,31,59,90,120,151,181,212,243,273,304,334,365};
    
    if (IsLeapYear(yy))
        for (int i=2; i<13; i++) ndys[i]+=1;
    
    int cnt=ndys[mm-1]+dd;
    return cnt;
}

void Sphen(int year, int day, double* s, double* p, double* h, double* en, double* p1)
{
    static double delt[] = {0.00, -5.04, -3.90, -2.87, -0.58,  0.71,  1.80,
        3.08,  4.63,  5.86,  7.21,  8.58, 10.50, 12.10, 12.49, 14.41, 15.59,
        15.81, 17.52, 19.01, 18.39, 19.55, 20.36, 21.01, 21.81, 21.76, 22.35,
        22.68, 22.94, 22.93, 22.69, 22.94, 23.20, 23.31, 23.63, 23.47, 23.68,
        23.62, 23.53, 23.59, 23.99, 23.80, 24.20, 24.99, 24.97, 25.72, 26.21,
        26.37, 26.89, 27.68, 28.13, 28.94, 29.42, 29.66, 30.29, 30.96, 31.09,
        31.59, 31.52, 31.92, 32.45, 32.91, 33.39, 33.80, 34.23, 34.73, 35.40,
        36.14, 36.99, 37.87, 38.75, 39.70, 40.70, 41.68, 42.82, 43.96, 45.00,
        45.98, 47.00, 48.03, 49.10, 50.10, 50.97, 51.81, 52.57};
    
    double cycle,t,td;
    double delta,deltat;
    long  yr,ilc,icent;
    long  it,ild,iday;
    long  iyd,ipos;
    
    cycle=360.000;
    ilc=0;
    icent=year / 100;
    yr=year-icent*100;
    t=icent-20.0;
    
    /*  For the following equations, time origin is fixed at 00hr of 1 Jan, 2000  */
    it=icent-20;
    if (it<0)
        iday=(it/4)-it;
    else
        iday=((it+3)/4)-it;
    
    /*  T is in julian century. Correction in gregorian calendar where only century
     year divisible by 4 is a leap year  */
    td=0.0;
    if (yr!=0) {
        iyd=365*yr;
        ild=(yr-1)/4;
        if ((icent-(icent/4)*4)==0)
            ilc=1;
        td=iyd+ild+ilc;
    }
    td=td+iday+day-1.5;
    t=t+(td/36525.0);
    deltat=0.0;
    ipos=year-1899;
    if (ipos>0) {
        if (ipos<83)
            delta=(delt[ipos+1]+delt[ipos])/2.0;
        else
            delta=(65.0-50.5)/20.0*(year-1980)+50.5;
        deltat=delta*1.0E-6;
    }
    
    *s = 218.3165 + 481267.8813*t - 0.0016*t*t + 152.0*deltat;
    *h = 280.4661 + 36000.7698*t  + 0.0003*t*t + 11.0*deltat;
    *p =  83.3535 + 4069.0139*t   - 0.0103*t*t + deltat;
    *en= 125.0445 - 1934.1363*t   + 0.0021*t*t - deltat;
    *p1= 282.9384 + 1.7195*t      + 0.0005*t*t;
    
    *s=dmod(*s,cycle);
    *h=dmod(*h,cycle);
    *p=dmod(*p,cycle);
    *en=dmod(*en,cycle);
    *p1=dmod(*p1,cycle);
}


void Vset(double S, double P, double P1, double H, double* V)
{
    double H2,H3,H4,P2,S2,S3,S4;
    
    H2=2.0*H;
    H3=3.0*H;
    H4=4.0*H;
    S2=2.0*S;
    S3=3.0*S;
    S4=4.0*S;
    P2=P+P;
    
    /*	V's computed in degrees  */
    V[0]   = 0; //Added by Nick Thorne
    V[1]   = H;
    V[2]   = H2;
    V[3]   = S-P;
    V[4]   = S2-H2;
    V[5]   = S2;
    V[6]   = H-S4+P2+270.0;
    V[7]   = H3-S4+270.0;
    V[8]   = H-S3+P+270.0;
    V[9]   = H3-S3-P+270.0;
    V[10]  = H-S2+270.0;
    V[11]  = H3-S2+90.0;
    V[12]  = H-S+90.0;
    V[13]  = H3-S-P+90.0;
    V[14]  = P1-H2+270.0;
    V[15]  = 270.0-H;
    V[16]  = 180.0;
    V[17]  = H+90.0;
    V[18]  = H2-P1+90.0;
    V[19]  = H3+90.0;
    V[20]  = S-H+P+90.0;
    V[21]  = S+H-P+90.0;
    V[23]  = S2+H+90.0;
    V[26]  = H2-S4+P2;
    V[27]  = H4-S4;
    V[28]  = H2-S3+P;
    V[29]  = H4-S3-P;
    V[31]  = H2-S2;
    V[32]  = H4-S2;
    V[33]  = P-S+180.0;
    V[34]  = H2-S-P+180.0;
    V[35]  = P1-H;
    V[36]  = 0.0;
    V[37]  = H-P1+180.0;
    V[38]  = H2;
    
    V[22]  = -V[10];
    V[24]  = V[10]+V[8];
    V[25]  = V[31]+V[28];
    V[30]  = V[10]+V[15];
    V[39]  = V[31]-V[28];
    V[40]  = V[17]+V[21];
    V[41]  = -V[31];
    V[42]  = V[31]+V[10];
    V[43]  = H3-S3+180.0;
    V[44]  = V[10];
    V[45]  = V[31]+V[17];
    V[46]  = V[17];
    V[47]  = V[25];
    V[48]  = V[31]+V[31];
    V[49]  = V[28];
    V[50]  = V[31];
    V[51]  = V[31]+V[38];
    V[52]  = 0.0;
    V[53]  = V[38];
    V[54]  = V[48]+V[28];
    V[55]  = V[48]+V[31];
    V[56]  = V[47];
    V[57]  = V[48];
    V[58]  = V[48]+V[38];
    V[59]  = V[31];
    V[60]  = V[51];
    V[61]  = V[54];
    V[62]  = V[55]-V[38];
    V[63]  = V[55];
    V[64]  = V[25]+V[38];
    V[65]  = V[28]-V[38];
    V[66]  = -V[38];
    V[67]  = V[48]-V[28]-V[28];
    V[68]  = V[31]+V[8];
    V[69]  = V[48]-V[15];
    V[70]  = V[48]-V[8];
    V[71]  = V[62];
    V[72]  = V[55];
    V[73]  = V[48]-V[38];
    V[74]  = V[55]-V[17];
    V[75]  = V[31]+V[43];
    V[76]  = V[55]-V[10];
    V[77]  = V[54]+V[28];
    V[78]  = V[55]+V[28];
    
    V[89]  = V[48]+V[48];
    
    V[79]  = V[89]-V[38];
    V[80]  = V[89];
    V[81]  = V[54]-V[38];
    V[82]  = V[48]+V[29];
    V[83]  = V[62];
    V[84]  = V[89]-V[28];
    V[85]  = V[55]-V[28];
    V[86]  = V[51]+V[34];
    V[87]  = V[77];
    V[88]  = V[78];
    
    V[90]  = V[54];
    V[91]  = V[55];
    V[92]  = V[55]+V[38];
    V[93]  = V[64];
    V[94]  = V[48];
    V[95]  = V[58];
    V[96]  = V[89];
    V[97]  = V[55];
    V[98]  = V[89]+V[28];
    V[99]  = V[89]+V[31];
    V[100] = V[89];
    V[101] = V[31]+V[29];
    V[102] = V[73];
    
    /*  V[103],V[104] changed according to Dr. Cartright's notes of 15/11/1977  */
    
    V[103] = V[31]-H;
    V[104] = V[31]+H;
    V[105] = V[31]-V[29];
    V[106] = V[38]-V[31];
    V[107] = V[54];
    V[108] = V[101];
    V[109] = V[85];
    V[110] = V[34];
    V[111] = V[28]+V[35];
    V[112] = V[28]-V[35];
    V[113] = V[42]-270.0;
    V[114] = V[45]-90.0;
    V[115] = V[34];
    
    for (int i=1; i<=115; i++)  /* Put in range [0..360] */
        V[i]=dmod(V[i],360.0);
}


void UFset(double P, double CN, double* U, double* F)
{
    double RAD=6.28318530717959/360.0;
    double DEG=360.0/6.28318530717959;
    double PI =3.14159265358979323846;
    
    double NW,PW;
    double W1,W2,W3,W4,W5,W6;
    double A1,A2,A3,A4,X,Y;
    
    int i;
    
    PW = P*RAD;
    NW = CN*RAD;
    
    W1 = cos(NW);
    W2 = cos(2.0*NW);
    W3 = cos(3.0*NW);
    W4 = sin(NW);
    W5 = sin(2.0*NW);
    W6 = sin(3.0*NW);
    
    A1 = PW-NW;
    A2 = 2.0*PW;
    A3 = A2-NW;
    A4 = A2-2.0*NW;
    
    
    /* U's are computed in radians */
    
    U[3] = 0.0;
    F[3] = 1.0 - 0.1300*W1 + 0.0013*W2;
    
    U[5] = -0.4143*W4 + 0.0468*W5 - 0.0066*W6;
    F[5] = 1.0429 + 0.4135*W1 - 0.004*W2;
    
    U[10]= 0.1885*W4 - 0.0234*W5 + 0.0033*W6;
    F[10]= 1.0089 + 0.1871*W1 - 0.0147*W2 + 0.0014*W3;
    
    X = 2.0*cos(PW)+0.4*cos(A1);
    Y = sin(PW)+0.2*sin(A1);
    
    U[12] = atan (Y/X);
    if(X<0.0) U[12] = U[12]+PI;
    F[12] = sqrt(X*X+Y*Y);
    
    U[17] = -0.1546*W4 + 0.0119*W5 - 0.0012*W6;
    F[17] = 1.0060 + 0.1150*W1 - 0.0088*W2 + 0.0006*W3;
    
    U[21] = -0.2258*W4 + 0.0234*W5 - 0.0033*W6;
    F[21] = 1.0129 + 0.1676*W1 - 0.0170*W2 + 0.0016*W3;
    
    F[23] = 1.1027 + 0.6504*W1 + 0.0317*W2 - 0.0014*W3;
    U[23] = -0.6402*W4 + 0.0702*W5 - 0.0099*W6;
    
    U[31] = -0.0374*W4;
    F[31] = 1.0004 - 0.0373*W1 + 0.0002*W2;
    
    X = 1.0-0.2505*cos(A2)-0.1102*cos(A3)-0.0156*cos(A4)-0.037*W1;
    Y = -0.2505*sin(A2)-0.1102*sin(A3)-0.0156*sin(A4)-0.037*W4;
    
    U[34] = atan (Y/X);
    if (X<0.0) U[34] = U[34]+PI;
    F[34] = sqrt (X*X+Y*Y);
    
    U[38] = -0.3096*W4 + 0.0119*W5 - 0.0007*W6;
    F[38] = 1.0241 + 0.2863*W1 + 0.0083*W2 - 0.0015*W3;
   
    U[0]   = 0; //Added by Nick Thorne
    U[1]   = 0.0;
    U[2]   = 0.0;
    U[4]   = -U[31];
    U[6]   = U[10];
    U[7]   = U[10];
    U[8]   = U[10];
    U[9]   = U[10];
    U[11]  = U[31];
    U[13]  = U[21];
    U[14]  = 0.0;
    U[15]  = 0.0;
    U[16]  = 0.0;
    U[18]  = 0.0;
    U[19]  = 0.0;
    U[20]  = U[21];
    U[22]  = -U[10];
    U[24]  = 2.0*U[10];
    U[25]  = 2.0*U[31];
    U[26]  = U[31];
    U[27]  = U[31];
    U[28]  = U[31];
    U[29]  = U[31];
    U[30]  = U[10];
    U[32]  = U[31]+U[38];
    U[33]  = U[31];
    U[35]  = 0.0;
    U[36]  = 0.0;
    U[37]  = 0.0;
    U[39]  = 0.0;
    U[40]  = U[17]+U[21];
    U[41]  = U[4];
    U[42]  = U[31]+U[10];
    U[43]  = U[31]*1.5;
    U[44]  = U[10];
    U[45]  = U[31]+U[17];
    U[46]  = U[17];
    U[47]  = U[25];
    U[48]  = U[25];
    U[49]  = U[31];
    U[50]  = U[31];
    U[51]  = U[32];
    U[52]  = 0.0;
    U[53]  = U[38];
    U[54]  = U[25]+U[31];
    U[55]  = U[54];
    U[56]  = U[25];
    U[57]  = U[25];
    U[58]  = U[25]+U[38];
    U[59]  = U[31];
    U[60]  = U[32];
    U[61]  = 0.0;
    U[62]  = U[54]-U[38];
    U[63]  = U[54];
    U[64]  = U[58];
    U[65]  = U[31]-U[38];
    U[66]  = -U[38];
    U[67]  = 0.0;
    U[68]  = U[42];
    U[69]  = U[25];
    U[70]  = U[25]-U[10];
    U[71]  = U[54]-U[38];
    U[72]  = U[54];
    U[73]  = U[25]-U[38];
    U[74]  = U[54]-U[17];
    U[75]  = 2.5*U[31];
    U[76]  = U[54]-U[10];
    U[77]  = 2.0*U[25];
    U[78]  = U[77];
    U[79]  = U[77]-U[38];
    U[80]  = U[77];
    U[81]  = U[71];
    U[82]  = U[54];
    U[83]  = U[71];
    U[84]  = U[54];
    U[85]  = U[25];
    U[86]  = U[51]+U[34];
    U[87]  = U[77];
    U[88]  = U[77];
    U[89]  = U[77];
    U[90]  = U[54];
    U[91]  = U[54];
    U[92]  = U[54]+U[38];
    U[93]  = U[58];
    U[94]  = U[25];
    U[95]  = U[58];
    U[96]  = U[77];
    U[97]  = U[54];
    U[98]  = 5.0*U[31];
    U[99]  = U[98];
    U[100] = U[77];
    U[101] = U[25];
    U[102] = U[73];
    
    /* U[103],U[104] changed according to Dr.Cartwright's notes 15/11/1977 */
    
    U[103] = 0.0;
    U[104] = 0.0;
    U[105] = 0.0;
    U[106] = -U[65];
    U[107] = U[54];
    U[108] = U[25];
    U[109] = 0.0;
    U[110] = U[31];
    U[111] = 0.0;
    U[112] = 0.0;
    U[113] = U[42];
    U[114] = U[45];
    U[115] = U[31];
    
    /* Convert into degrees */
    
    for (i=1; i<=115; i++)
        U[i]=dmod(U[i]*DEG,360.0);
    F[0]   = 0.0; //added by Nick Thorne
    F[1]   = 1.0;
    F[2]   = 1.0;
    F[4]   = F[31];
    F[6]   = F[10];
    F[7]   = F[10];
    F[8]   = F[10];
    F[9]   = F[10];
    
    F[11]  = F[31];
    F[13]  = F[21];
    F[14]  = 1.0;
    F[15]  = 1.0;
    F[16]  = 1.0;
    F[18]  = 1.0;
    F[19]  = 1.0;
    F[20]  = F[21];
    F[22]  = F[10];
    F[24]  = F[10]*F[10];
    F[25]  = F[31]*F[31];
    F[26]  = F[31];
    F[27]  = F[31];
    F[28]  = F[31];
    F[29]  = F[31];
    F[30]  = F[10];
    F[32]  = F[31]*F[38];
    F[33]  = F[31];
    F[35]  = 1.0;
    F[36]  = 1.0;
    F[37]  = 1.0;
    F[39]  = F[25];
    F[40]  = F[17]*F[21];
    F[41]  = F[31];
    F[42]  = F[31]*F[10];
    F[43]  = pow(F[31],1.5);
    F[44]  = F[10];
    F[45]  = F[31]*F[17];
    F[46]  = F[17];
    F[47]  = F[25];
    F[48]  = F[25];
    F[49]  = F[31];
    F[50]  = F[31];
    F[51]  = F[32];
    F[52]  = 1.0;
    F[53]  = F[38];
    F[54]  = F[25]*F[31];
    F[55]  = F[54];
    F[56]  = F[25];
    F[57]  = F[25];
    F[58]  = F[25]*F[38];
    F[59]  = F[31];
    F[60]  = F[32];
    F[61]  = 1.0;
    F[62]  = F[54]*F[38];
    F[63]  = F[54];
    F[64]  = F[58];
    F[65]  = F[32];
    F[66]  = F[38];
    F[67]  = F[25]*F[25];
    F[68]  = F[42];
    F[69]  = F[25];
    F[70]  = F[25]*F[10];
    F[71]  = F[54]*F[38];
    F[72]  = F[54];
    F[73]  = F[58];
    F[74]  = F[54]*F[17];
    F[75]  = 0.5*(F[25]+F[54]);
    F[76]  = F[54]*F[10];
    F[77]  = F[67];
    F[78]  = F[67];
    F[79]  = F[67]*F[38];
    F[80]  = F[67];
    F[81]  = F[71];
    F[82]  = F[54];
    F[83]  = F[71];
    F[84]  = F[54]*F[25];
    F[85]  = F[67];
    F[86]  = F[51]*F[34];
    F[87]  = F[67];
    F[88]  = F[67];
    F[89]  = F[67];
    F[90]  = F[54];
    F[91]  = F[54];
    F[92]  = F[71];
    F[93]  = F[58];
    F[94]  = F[25];
    F[95]  = F[58];
    F[96]  = F[67];
    F[97]  = F[54];
    F[98]  = F[84];
    F[99]  = F[84];
    F[100] = F[67];
    F[101] = F[25];
    F[102] = F[58];
    
    /*	F[103],F[104] changed according to Dr.Cartwright's notes 15/11/1977   */
    F[103] = 1.0;
    F[104] = 1.0;
    F[105] = 1.0;
    F[106] = F[32];
    F[107] = F[54];
    F[108] = F[25];
    F[109] = 1.0;
    F[110] = F[54];
    F[111] = 1.0;
    F[112] = 1.0;
    F[113] = F[42];
    F[114] = F[45];
    F[115] = F[54];
}


void SigmaSet(double* SI)
{
    SI[0]   = 0.0; //added by Nick Thorne
    SI[1]   = 0.410686E-01;
    SI[2]   = 0.821373E-01;
    SI[3]   = 0.5443747E+00;
    SI[4]   = 0.10158958E+01;
    SI[5]   = 0.10980331E+01;
    SI[6]   = 0.128542862E+02;
    SI[7]   = 0.129271398E+02;
    SI[8]   = 0.133986609E+02;
    SI[9]   = 0.134715145E+02;
    SI[10]  = 0.139430356E+02;
    SI[11]  = 0.140251729E+02;
    SI[12]  = 0.144920521E+02;
    SI[13]  = 0.145695476E+02;
    SI[14]  = 0.149178647E+02;
    SI[15]  = 0.149589314E+02;
    SI[16]  = 0.150000000E+02;
    SI[17]  = 0.150410686E+02;
    SI[18]  = 0.150821353E+02;
    SI[19]  = 0.151232059E+02;
    SI[20]  = 0.155125897E+02;
    SI[21]  = 0.155854433E+02;
    SI[22]  = 0.160569644E+02;
    SI[23]  = 0.161391017E+02;
    SI[24]  = 0.273416965E+02;
    SI[25]  = 0.274238337E+02;
    SI[26]  = 0.278953548E+02;
    SI[27]  = 0.279682084E+02;
    SI[28]  = 0.284397295E+02;
    SI[29]  = 0.285125831E+02;
    SI[30]  = 0.289019669E+02;
    SI[31]  = 0.289841042E+02;
    SI[32]  = 0.290662415E+02;
    SI[33]  = 0.294556253E+02;
    SI[34]  = 0.295284789E+02;
    SI[35]  = 0.299589333E+02;
    SI[36]  = 0.300000000E+02;
    SI[37]  = 0.300410667E+02;
    SI[38]  = 0.300821373E+02;
    SI[39]  = 0.305443747E+02;
    SI[40]  = 0.306265120E+02;
    SI[41]  = 0.310158958E+02;
    SI[42]  = 0.429271398E+02;
    SI[43]  = 0.434761563E+02;
    SI[44]  = 0.439430356E+02;
    SI[45]  = 0.440251729E+02;
    SI[46]  = 0.450410686E+02;
    SI[47]  = 0.574238337E+02;
    SI[48]  = 0.579682084E+02;
    SI[49]  = 0.584397295E+02;
    SI[50]  = 0.589841042E+02;
    SI[51]  = 0.590662415E+02;
    SI[52]  = 0.600000000E+02;
    SI[53]  = 0.600821373E+02;
    SI[54]  = 0.864079380E+02;
    SI[55]  = 0.869523127E+02;
    SI[56]  = 0.874238337E+02;
    SI[57]  = 0.879682084E+02;
    SI[58]  = 0.880503457E+02;
    SI[59]  = 0.889841042E+02;
    SI[60]  = 0.890662415E+02;
    SI[61]  = 0.264079379E+02;
    SI[62]  = 0.268701754E+02;
    SI[63]  = 0.269523127E+02;
    SI[64]  = 0.275059710E+02;
    SI[65]  = 0.283575922E+02;
    SI[66]  = 0.299178627E+02;
    SI[67]  = 0.310887494E+02;
    SI[68]  = 0.423827651E+02;
    SI[69]  = 0.430092770E+02;
    SI[70]  = 0.445695475E+02;
    SI[71]  = 0.568701754E+02;
    SI[72]  = 0.569523127E+02;
    SI[73]  = 0.578860711E+02;
    SI[74]  = 0.719112441E+02;
    SI[75]  = 0.724602605E+02;
    SI[76]  = 0.730092771E+02;
    SI[77]  = 0.848476674E+02;
    SI[78]  = 0.853920422E+02;
    SI[79]  = 0.858542795E+02;
    SI[80]  = 0.859364168E+02;
    SI[81]  = 0.863258006E+02;
    SI[82]  = 0.864807915E+02;
    SI[83]  = 0.868701754E+02;
    SI[84]  = 0.874966873E+02;
    SI[85]  = 0.885125832E+02;
    SI[86]  = 0.885947204E+02;
    SI[87]  = 0.1148476674E+03;
    SI[88]  = 0.1153920422E+03;
    SI[89]  = 0.1159364168E+03;
    SI[90]  = 0.1164079379E+03;
    SI[91]  = 0.1169523127E+03;
    SI[92]  = 0.1170344500E+03;
    SI[93]  = 0.1175059710E+03;
    SI[94]  = 0.1179682084E+03;
    SI[95]  = 0.1180503457E+03;
    SI[96]  = 0.1459364168E+03;
    SI[97]  = 0.1469523127E+03;
    SI[98]  = 0.1743761463E+03;
    SI[99]  = 0.1749205210E+03;
    SI[100] = 0.1759364168E+03;
    SI[101] = 0.274966873E+02;
    SI[102] = 0.278860711E+02;
    SI[103] = 0.289430356E+02;
    SI[104] = 0.290251728E+02;
    SI[105] = 0.304715211E+02;
    SI[106] = 0.310980331E+02;
    SI[107] = 0.564079379E+02;
    SI[108] = 0.574966873E+02;
    SI[109] = 0.585125830E+02;
    SI[110] = 0.595284789E+02;
    SI[111] = 0.283986609E+02;
    SI[112] = 0.284807981E+02;
    SI[113] = 0.729271398E+02;
    SI[114] = 0.740251728E+02;
    SI[115] = 0.295284789E+02;
    
}

extern "C"
{

int POLPredict( double * resultArrayPtr, double duration, int nvals, double * H,double * G,int * K, double ref_height, int Base_dd, int Base_mm, int Base_yy, float startTime, float *maxRangePtr, int constit_count)  //constit_count = number of constituents (15 for currents 40 for tides)
{
    //printf("Tidal Prediction Code\n=====================\n\n");
    //int Base_dd, Base_mm, Base_yy;
    
    // We assume if 40 consitutents then is a tide calc, else currents.
    int isCurrents = 0;
    // Currents Doodson numbers are always the same for currents
    
    int KDoodsonCurrents[MAXCONST] = {8, 10,  15, 16, 17, 26, 27, 28, 29, 31, 34, 35, 36, 38, 48};
    
  //  if (constit_count !=40)
  //  {
        isCurrents= 1;
        //Now redirect the K pointer to our new Doodson array
      //  Ki = KDoodsonCurrents;
  //  }
   
    double t;
    double RAD = 6.28318530717959/360.0;
    //double H[MAXCONST]={0};
    //double G[MAXCONST]={0};
    //int K[MAXCONST]={0};
    int ndc;
    int i;
    //double z0;
    double V[MAXCONST],U[MAXCONST];
    double F[MAXCONST],Sig[MAXCONST];
    double s,p,hh,en,p1;
    double Z; 
    

    int n=0;
    
   // K = KDoodsonCurrents;


    ndc=constit_count;  /* ndc is the number of Doodson Constituents - sample data is 40 HCs from Christchurch */

    
    if (1 == 0)
    {
        printf("Assigned results2 :\n");
        for (n=0;n<constit_count;n++)
        {
            printf("> H[%d]=%f; G[%d]=%f; K[%d]=%d\n",n,H[n],n,G[n],n,K[n]);
        }
        printf("\n");
    }
    
    /* Compute the basic V,U,F and Sigma tables */
    int jul=ddmm2jul(Base_dd,Base_mm,Base_yy);   /* jul=Julian DayNum 1=1-Jan ; 365/366=31-Dec */
    Sphen(Base_yy,jul,&s,&p,&hh,&en,&p1);
    SigmaSet(Sig);
    Vset(s,p,p1,hh,V);
    UFset(p,en,U,F);


    /* Need to write a function that computes t, number of hours required. Time t is from the base time.
     for example if Base date is 1/1/2000,
     t=12.5  would be 1/1/2000 12:30pm
     t=30.25 would be 2/1/2000 6:15am  */
    
    /* Convert into radians */
    for (i=0; i<ndc; i++) G[i]*=RAD;
    for (i=1; i<MAXCONST; i++) {
        V[i]*=RAD;
        U[i]*=RAD;
        Sig[i]*=RAD;
    }

    /* Alter H and G to include the V,U and f values (shifts to a common time base) 
     This will increase calculation speed later  */
    for (i=0; i<ndc; i++) {
	//printf("_________________\n");

	//printf("H[%d] = %lf \n",i,H[i]);
	//printf("G[%d] = %lf \n",i,G[i]);
	//printf("......\n");

        H[i]*=F[K[i]];
        G[i]=-G[i]+V[K[i]]+U[K[i]];
         //printf("\nk[%d] = %d\n",i,K[i]);
         //printf("V[%d]=%lf\n",K[i],V[K[i]]);
         //printf("U[%d] =%lf\n", K[i],U[K[i]]);
         //printf("G[%d] = %lf \n",i,G[i]);
         //printf("V[0]=%lf\n",V[0]);
	 //printf("Sig[%d] =%lf\n", i,Sig[i]);
    }
 	
    /* TIME SERIES LOOP */
    int resultIndex=0;
    //nvals = 193; //number of points
    float incrementVal = duration/(nvals); // -1 because we have 1 extra point to plot. E.g. 
    // 2 hours has 3 points, 0, 1, 2

    //return 0;

    for (t=0.0; t<duration; t+=incrementVal) {  //loops 193 for .125, 385 for .0625
        
        /* Loop to sum up each constituent */
        Z=ref_height;   

        for (i=0; i<ndc; i++) {

            Z+=(H[i]*cos(Sig[K[i]]*(t+startTime) + G[i] ));
	    //printf("Z=%6.2f\n",Z);
        }

        //printf("t=%6.5f  Z=%6.5f\n",(t+ startTime),Z);
	//printf("%6.5f,%6.5f\n",(t+ startTime),Z);


	//resultArrayPtr[resultIndex] =  Z;
        *(resultArrayPtr+resultIndex) =  Z; //test pattern = ((int)(t/4)) +1;
        resultIndex++;
	

    }

    //printf("Number of points = %d\n",resultIndex);
    //printf("resultArrayPtr size: %lu", sizeof(*resultArrayPtr));

    return 0;
}

}

//-------------------------------------------------------------------------
int split_HGK (char *str, double * H,double * G,int * K)
{
  char* pEnd;
  char * pch;
  int counter = 0;
  int position = 0;
	
  pch = strtok (str,",");
  while (pch != NULL)
  {
    //printf ("%s\n",pch);
    double tempD = strtod (pch, &pEnd);

    if(counter == 0) { H[position] = tempD; }
    if(counter == 1) { G[position] = tempD; }	
    if(counter == 2) { K[position] = static_cast<int>(tempD); }	

    if(counter == 2) {
	counter = 0;
	position++;
    } else {	    				
    	counter++;	
    }


    pch = strtok (NULL, ",");
  }

  return 0;
}



static PyObject* predictor(PyObject* self, PyObject* args)
{
	
	int arr_size = 40;

	double H[arr_size];
	double G[arr_size];
	int K[arr_size];	

	const char *lat, *lng, *date, *hgk;
	
	PyArg_ParseTuple(args, "ssss", &date, &lat, &lng, &hgk);

	split_HGK(const_cast<char *>(hgk), H, G, K);

	//generate date and time
	std::string date_full = std::string(date);
	std::string year_str = date_full.substr(0,4);
	std::string month_str = date_full.substr (5,2);
	std::string day_str = date_full.substr (8,2);
	int year = atoi(year_str.c_str());
	int month = atoi(month_str.c_str());
	int day = atoi(day_str.c_str());
	int Base_dd = day; // from script comments
	int Base_mm = month; // from script comments
	int Base_yy = year; // from script comments
	float startTime = 0; //this could be timestamp from base datetime ?

	double resultArrayPtr[100];   // this is var where final results are stored
	double duration = 24;         // 24, from script comments
	int nvals = 96;               // ?  //nvals = 193, from script comments
	double ref_height = 0; // we have this value in json data (API)
	float maxRangePtr[arr_size];  // 0 
	int constit_count = 40;  // 40 -> we know this from json data (API)

	/*
	std::cout << "Date: " << date_full;
	std::cout << "Data: " << hgk;

	for (int i = 0; i < arr_size; ++i)
 	{ 
		printf("%f %f %i\n", H[i], G[i], K[i]);
	}
	*/

	int z35 = POLPredict(resultArrayPtr, duration, nvals, H, G, K, ref_height, Base_dd, Base_mm, Base_yy, startTime, maxRangePtr, constit_count);
   	/*--------------------------------------------------*/
	/*std::cout << "\n---------------------------\n";
	for(int j=0; j<24; j++) {
		std::cout << resultArrayPtr[j];
		std::cout << "\n";
	}*/



	PyObject* tuple = PyTuple_New(nvals); 

	for (int i = 0; i < nvals; i++) { 
	double rVal = resultArrayPtr[i]; 
		PyTuple_SetItem(tuple, i, Py_BuildValue("d", rVal)); 
	} 

	return tuple; 

	

    //return Py_BuildValue("[d]", resultArrayPtr);
}


static PyMethodDef PPpredictor_funcs[] = {
    {"predictor", (PyCFunction)predictor, METH_VARARGS},	
    {NULL}
};

extern "C" void initpolpredict(void)
{
    Py_InitModule3("polpredict", PPpredictor_funcs,
                   "Extension Pol Predictor module");
}


// build: c++ -o polp POLPredict.cpp
/*
int main(int argc, char * argv[]) {
	
	int arr_size = 40;
			      
	double resultArrayPtr[200];   // this is var where final results are stored
	double duration = 24;         // 24, from script comments
	int nvals = 96;               // ?  //nvals = 193, from script comments

	double ref_height = 0; // we have this value in json data (API)

    	double H[arr_size];
    	double G[arr_size];
    	int K[arr_size];


	char* date1 = argv[1];
	char* lat = argv[3];
	char* lng = argv[4];

	std::string date_full = std::string(date1);
	std::string year_str = date_full.substr(0,4);
	std::string month_str = date_full.substr (5,2);
	std::string day_str = date_full.substr (8,2);

	int year = atoi(year_str.c_str());
	int month = atoi(month_str.c_str());
	int day = atoi(day_str.c_str());


	int Base_dd = day; // from script comments
	int Base_mm = month; // from script comments
	int Base_yy = year; // from script comments
	float startTime = 0; //this could be timestamp from base datetime ?



	float maxRangePtr[arr_size];  // 0 
	int constit_count = 40;  // 40 -> we know this from json data (API)


	int c1 = 0;

	for (int i = 5; i < 45; ++i)
 	{ 
				
		H[i-5] = strtod(argv[c1+5], NULL);
		G[i-5] = strtod(argv[c1+1+5], NULL);
		K[i-5] = strtod(argv[c1+2+5], NULL);

		//printf("%s %f\n", argv[c1+5], H[i-5]);

		c1 = c1 + 3;

	}


	int z35 = POLPredict(resultArrayPtr, duration, nvals, H, G, K, ref_height, Base_dd, Base_mm, Base_yy, startTime, maxRangePtr, constit_count);
	
	//printf("resultArrayPtr size: %lu", sizeof(resultArrayPtr));	
	//printf("%lf", *resultArrayPtr);
        return 0;
}
*/

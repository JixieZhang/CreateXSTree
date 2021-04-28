//awk -v zz=-10.0 'BEGIN {for(j=1;j<=8;j++) { for(i=1;i<=10;i++) {printf("%6.2f,",zz);zz=zz+0.5;};printf("\n    ");}}'
//awk -v zz=-25.0 'BEGIN {for(j=1;j<=8;j++) { for(i=1;i<=10;i++) {printf("%6.2f,",zz);zz=zz+1.0;};printf("\n    ");}}'
//awk -v zz=-0.08 'BEGIN {for(j=1;j<=9;j++) { for(i=1;i<=10;i++) {printf("%6.3f,",zz);zz=zz+0.002;};printf("\n    ");}}'
//awk -v zz=-0.05 'BEGIN {for(j=1;j<=8;j++) { for(i=1;i<=10;i++) {printf("%6.3f,",zz);zz=zz+0.002;};printf("\n    ");}}'

#ifndef _ACCEPTANCE_EXC_BIN_
#define _ACCEPTANCE_EXC_BIN_

namespace ACCEPTANCE {

  static const int kYtarBinNum=64;
  static const int kDeltaBinNum=60;
  static const int kThetaBinNum=80;
  static const int kPhiBinNum=50;

  //-12 < SHMS < 18, -14 < HMS < 14
  static double kYtarBin[]={
    -14.00,-13.50,-13.00,-12.50,-12.00,-11.50,-11.00,-10.50,-10.00, -9.50,
     -9.00, -8.50, -8.00, -7.50, -7.00, -6.50, -6.00, -5.50, -5.00, -4.50,
     -4.00, -3.50, -3.00, -2.50, -2.00, -1.50, -1.00, -0.50,  0.00,  0.50,
      1.00,  1.50,  2.00,  2.50,  3.00,  3.50,  4.00,  4.50,  5.00,  5.50,
      6.00,  6.50,  7.00,  7.50,  8.00,  8.50,  9.00,  9.50, 10.00, 10.50,
     11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00, 14.50, 15.00, 15.50,
     16.00, 16.50, 17.00, 17.50, 18.00
  };

  static double kDeltaBin[]={
    -25.00,-24.00,-23.00,-22.00,-21.00,-20.00,-19.00,-18.00,-17.00,-16.00,
    -15.00,-14.00,-13.00,-12.00,-11.00,-10.00, -9.00, -8.00, -7.00, -6.00,
     -5.00, -4.00, -3.00, -2.00, -1.00,  0.00,  1.00,  2.00,  3.00,  4.00,
      5.00,  6.00,  7.00,  8.00,  9.00, 10.00, 11.00, 12.00, 13.00, 14.00,
     15.00, 16.00, 17.00, 18.00, 19.00, 20.00, 21.00, 22.00, 23.00, 24.00,
     25.00, 26.00, 27.00, 28.00, 29.00, 30.00, 31.00, 32.00, 33.00, 34.00,
     35.00
  };

  static  double kThetaBin[]={
    -0.080,-0.078,-0.076,-0.074,-0.072,-0.070,-0.068,-0.066,-0.064,-0.062,
    -0.060,-0.058,-0.056,-0.054,-0.052,-0.050,-0.048,-0.046,-0.044,-0.042,
    -0.040,-0.038,-0.036,-0.034,-0.032,-0.030,-0.028,-0.026,-0.024,-0.022,
    -0.020,-0.018,-0.016,-0.014,-0.012,-0.010,-0.008,-0.006,-0.004,-0.002,
     0.000, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018,
     0.020, 0.022, 0.024, 0.026, 0.028, 0.030, 0.032, 0.034, 0.036, 0.038,
     0.040, 0.042, 0.044, 0.046, 0.048, 0.050, 0.052, 0.054, 0.056, 0.058,
     0.060, 0.062, 0.064, 0.066, 0.068, 0.070, 0.072, 0.074, 0.076, 0.078,
     0.080
  };

  static  double kPhiBin[]={
    -0.050,-0.048,-0.046,-0.044,-0.042,-0.040,-0.038,-0.036,-0.034,-0.032,
    -0.030,-0.028,-0.026,-0.024,-0.022,-0.020,-0.018,-0.016,-0.014,-0.012,
    -0.010,-0.008,-0.006,-0.004,-0.002, 0.000, 0.002, 0.004, 0.006, 0.008,
     0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024, 0.026, 0.028,
     0.030, 0.032, 0.034, 0.036, 0.038, 0.040, 0.042, 0.044, 0.046, 0.048,
     0.050
  };

};

#endif //_ACCEPTANCE_EXC_BIN_


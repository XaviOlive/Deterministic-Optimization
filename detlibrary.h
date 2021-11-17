#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

//---MACROS----
#define INITIALx_1 -1.5             //Initial position indicated of x_1 (seed)
#define INITIALx_2 -1               //Initial position indicated of x_2 (seed)
#define EPSILON 0.000001            //Epsilon factor stablished in the program (realted with accuracy desired)
#define RR 3                        //Replacement factor of lambda for the L-M method.

double solve_Rosen_function(double position[2]) {               //Calculation of the value of the function
    double x,y;
    x=position[0]; y = position[1];
    return 100.0*pow(y - x*x,2) + pow(1 - x,2);
}

void Gradient_calc(double position[2], double Gradient[2]) {    //Calculation of the gradient of the function
    double x,y;
    x=position[0]; y = position[1];
    Gradient[0] = -400*x*(y-x*x) - 2*(1-x);
    Gradient[1] = 200*(y-x*x);
}

void Hessian_calc(double position[2], double Hessian[2][2]) {   //Calculation of the Hessian of the function
    double x,y;
    x=position[0]; y = position[1];
    Hessian[0][0] = 1200*x*x - 400*y + 2; Hessian[0][1] = -400*x;
    Hessian[1][0] = -400*x; Hessian[1][1] = 200;
}

void d_CGM_calc(double d[2], double g[2], double betta) {       //Calculation of the direction in which the next step will be done for the Conjugate gradient method
    d[0]= -g[0]+betta*d[0];
    d[1]= -g[1]+betta*d[1];
}

double alpha_calc(double d[2], double g[2], double H[2][2]) {   //Calculation of alpha
    double alpha;
    alpha=(-g[0]*d[0]-g[1]*d[1]) / ((d[0]*H[0][0]+d[1]*H[1][0])*d[0] + (d[0]*H[0][1]+d[1]*H[1][1])*d[1]);
    return alpha;
}

void step_CGM(double d[2], double alpha, double p[2]){          //Calculation of the new step for Conjugate gradient method
    p[0]=p[0]+alpha*d[0];
    p[1]=p[1]+alpha*d[1];
}

void p_LM_calc2(double q[2], double H[2][2], double g[2], double lambda) { //h calculation for Levenberg–Marquardt Method (equivalent to the direction d)
    double M[2][2];
    double det;
    det=(H[0][0]+lambda)*(H[1][1]+lambda)-M[0][1]*M[1][0];
    M[0][0]=(H[1][1]+lambda)/det; M[0][1]=-H[1][0]/det;
    M[1][0]=-H[0][1]/det; M[1][1]=(H[0][0]+lambda)/det;
    q[0]= -(M[0][0]*g[0]+M[0][1]*g[1]);
    q[1]= -(M[1][0]*g[0]+M[1][1]*g[1]);
}

void step_LM(double d[2], double p[2], double p_new[2]){        //Calculation of the new step for Levenberg–Marquardt Method
    p_new[0]=p[0]+d[0];
    p_new[1]=p[1]+d[1];
}

double rho_calculation(double p[2], double p_new[2], double d[2], double g[2], double lambda){ //Calculation of rho
    double rho;
    rho=(solve_Rosen_function(p)-solve_Rosen_function(p_new)) /( lambda*(d[0]*d[0]+d[1]*d[1])/2-(d[0]*g[0]+d[1]*g[1]));
    return rho;
}

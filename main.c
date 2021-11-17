//  main.c         Deterministic          by Xavi Olivé Fernández on 03/01/2021.

#include "detlibrary.h"


int main(int argc, char* argv[]){

    //DECLARATION OF VARIABLES USED:
    FILE * Cord; //File used to save the steps done
    int i=0;
    double x[2], xx[2], direction[2];                   //vectors positions and direction
    double grad[2], Hes[2][2];                          //Calcul elements (first and second derivatives)
    double betta, alpha, rho, lambda;                   //Parameters
    double residu, oldresidu=1, tol;                    //Tolerance and residus parameters
    clock_t start, end;                                 //Time releated parameters
   
    x[0] = INITIALx_1; x[1] = INITIALx_2;             //Initial position
    
    //INITIAL CONDITIONS (for the two methods)
    Gradient_calc(x, grad);
    Hessian_calc(x, Hes);
    direction[0]=-grad[0]; direction[1]=-grad[1];
    residu=grad[0]*grad[0]+grad[1]*grad[1];
    tol = residu * EPSILON * EPSILON ;                     //Tolerance is got considerin epsilon but also the residu for the intital (i=0) case.
    lambda=Hes[0][0];                               //Lambda is got as the maximum term of the Hessian (as assumed tradionally)
    
    Cord = fopen("Coordenades.txt", "w");             //File creation where the steps will be stored
    fprintf (Cord, "%.4f\t%.4f\n",x[0],x[1]);
    
    start = clock();                                  //Clock start counter
    
    if (*argv[1] == 'g') {                            //METHODE 1: CONJUGATE GRADIENT METHOD
        printf("\t CONJUGATE GRADIENT METHOD SELECTED\n\t .\n\t .\n\t .\n");
        
        while (i < 10000 && residu > tol) {
            
            alpha=alpha_calc(direction, grad, Hes);
            step_CGM(direction, alpha, x);            //Calculation of the first step (and rpevioulsy alpha_k)
            
            Gradient_calc(x, grad);                   //Update of the gradient and the Hessian (as they have dependence of (x and y)
            Hessian_calc(x, Hes);                     //They will be used in the previous steps
            
            oldresidu=residu;
            residu=grad[0]*grad[0]+grad[1]*grad[1];   //Calculation of the residus
            betta=residu/oldresidu;                   //Calculation of the parametre beta_k with the residus (as the residus are combinations of the gradient needed to calculate beta)
            d_CGM_calc(direction, grad, betta);       //Update of the new direction
            
            fprintf (Cord, "%.4f\t%.4f\n",x[0],x[1]); //printing the results got (steps in our file)
            i++;
        }
    }
    else if (*argv[1] == 'l') {                       //METHODE 2: LEVENBERG-MARQUARDT METHOD
        printf("\t LEVENBERG-MARQUARDT METHOD SELECTED\n\t .\n\t .\n\t .\n");
        
        while (i < 10000 && residu > tol) {
            
            p_LM_calc2(direction, Hes, grad, lambda);
            step_LM(direction, x, xx);               //Calculation of the new direction and the new step (as we need the old values to update rho, we save them).
            rho=rho_calculation(x, xx, direction, grad, lambda); //Verification of rho which will be use to decide whicc action to do next
            
            if (rho > 0) {                          //Step accetable
                x[0]=xx[0]; x[1]=xx[1];             //We actualiza the new step
                Gradient_calc(x, grad);
                Hessian_calc(x, Hes);
                residu=grad[0]*grad[0]+grad[1]*grad[1];
                lambda=lambda*fmax(1/RR,1-pow((2*rho-1),3)); //We actualize the gradient, the hessian, the residu and lambda (which is reduced at least for 1/RR)
            }
            else {                                  //Step not accepted. Lambda then is increased in order to get a better direction for the same position we had.
                lambda=RR*lambda;
            }
            fprintf (Cord, "%.4f\t%.4f\n",x[0],x[1]); //printing the results got (steps in our file)
            i++;
        }
    }
    else {
        printf("\t NO METHOD SELECTED\n\n\n");
        exit(0);                                    //When no correct argument selected we leave the program
    }
                 
    end = clock();                                  //Clock end counter
    
    printf("\t CONVERGENCE COMPLETED\n\n");
    //RESULTS PRINTED ON SCREEN
    printf("FINAL RESULT:\n -position got: (%.9f, %.9f)\n -function value: %.9f\n -number of iterations: %d\n -time needed: %.4fms\n", x[0], x[1], solve_Rosen_function(x), i, ((double) end-start)*1000/CLOCKS_PER_SEC);
    fclose (Cord);
    return 0;
}


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define R_W 0.33
#define L_F 3
#define N 4
#define halfWidth 0.5
double dJ[4] = {0,0,0,0}; // gradient of cost function
double d2J[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; // hessian of cost function 
double A[4][4];  // quadatric term
double B[4]; // linear term
double C; // const term
double dBarrier[4]; // barrier gradient 
double W_v[2][2] = {{70,0},{0,100}}; // should always be [W1 0; 0 W2]


// motor and diff limits
const double frontLimit = 15*13.5; // front motor
const double rearLimit = 130*3.23; // rear motor
const double minDiff = 1;
const double bias_ratio = 4; // ratio of more traction wheel to less traction wheel torque

// barrier function
const double t = 0.000001; // strength of barrier

// dotProduct
double dotProduct(double * a, double * b, int size) {
    double result = 0;

    for(int i = 0; i < size;i ++) {
        result += a[i]*b[i];
    }
    return result;
}
// compute A and B and C
void initialize(double steering_angle, double normal_force[4],double total_torque, double yaw_moment_request) {
    // steering angle to wheel angle conversion
    double leftWheelAngle = steering_angle;
    double rightWheelAngle = steering_angle;
    double B_w[4];
    double B_v[2][4] = {{1,1,1,1},{(L_F*sin(leftWheelAngle)-halfWidth*cos(leftWheelAngle))/R_W,(L_F*sin(rightWheelAngle)+halfWidth*cos(rightWheelAngle))/R_W,-halfWidth/R_W,halfWidth/R_W}};
    
    for(int i = 0;i < 4; i++) {
        B_w[i] = 1/normal_force[i]/R_W;
    }
    double W_vB_v[2][4];
    double W_vc[2];
    W_vc[0] = W_v[0][0]*total_torque;
    W_vc[1] = W_v[1][1]*yaw_moment_request;

    for(int k = 0; k < 2; k++) { 
        for (int j = 0; j < 4; j++) {
            W_vB_v[k][j] = W_v[k][k]*B_v[k][j];
       }
    }
    // B
    for(int j = 0; j < 4;j++) {
        B[j] = W_vc[0]*W_vB_v[0][j]+W_vc[1]*W_vB_v[1][j];
    }
    // store A
    // initialize d2J
    for(int i = 0; i < 4; i++) {
        for(int j = i; j < 4; j++) {
            A[i][j] = W_vB_v[0][i]*W_vB_v[0][j]+W_vB_v[1][i]*W_vB_v[1][j]+B_w[i]*B_w[j];
            // adding identity matrix
            if(i == j) {
                A[i][j] += 1;
                d2J[i][j] = 2*A[i][j];
            }
            // symmetric matrix
            else{
                A[j][i] = A[i][j];
                d2J[i][j] = 2*A[i][j];
                d2J[j][i] = 2*A[i][j];
            }
        }
    }
    C = W_vc[0]*W_vc[0]+W_vc[1]*W_vc[1];
    return;
}

void dJCalc(double * u, double * torqueLimit) {
    double dBarrier[4];
    double f1 = (u[0]-torqueLimit[0]);
    double f2 = (-u[0]+torqueLimit[1]);
    double f3 = (u[1]-torqueLimit[2]);
    double f4 = (-u[1]+torqueLimit[3]);
    double f5 = (u[2]-torqueLimit[4]);
    double f6 = -u[2]+minDiff;
    double f7 = (u[3]-torqueLimit[5]);
    double f8 = -u[3]+minDiff;
    double f9 = (u[2]+u[3]-rearLimit);
    double f10 = (u[2]-bias_ratio*u[3]);
    double f11 = (-bias_ratio*u[2]+u[3]);


    dBarrier[0] = 1/(-f1)-1/(-f2);
    dBarrier[1] = 1/(-f3)-1/(-f4);
    dBarrier[2] = 1/(-f5)-1/(-f6)+1/(-f9)+(1/(-f10)-bias_ratio/(-f11));
    dBarrier[3] = 1/(-f7)-1/(-f8)+1/(-f9)-(bias_ratio/(-f10)+1/(-f11));
    
    for(int i = 0; i < 4; i++) {
        dJ[i] = 2*(A[i][0]*u[0]+A[i][1]*u[1]+A[i][2]*u[2]+A[i][3]*u[3]- B[i]) + dBarrier[i]/t;
    }
}

void dJd2JCalc(double * u, double * torqueLimit) {
    double dBarrier[4];
    double f1 = (u[0]-torqueLimit[0]);
    double f2 = (-u[0]+torqueLimit[1]);
    double f3 = (u[1]-torqueLimit[2]);
    double f4 = (-u[1]+torqueLimit[3]);
    double f5 = (u[2]-torqueLimit[4]);
    double f6 = -u[2]+minDiff;
    double f7 = (u[3]-torqueLimit[5]);
    double f8 = -u[3]+minDiff;
    double f9 = (u[2]+u[3]-rearLimit);
    double f10 = (u[2]-bias_ratio*u[3]);
    double f11 = (-bias_ratio*u[2]+u[3]);


    dBarrier[0] = 1/(-f1)-1/(-f2);
    dBarrier[1] = 1/(-f3)-1/(-f4);
    dBarrier[2] = 1/(-f5)-1/(-f6)+1/(-f9)+(1/(-f10)-bias_ratio/(-f11));
    dBarrier[3] = 1/(-f7)-1/(-f8)+1/(-f9)-(bias_ratio/(-f10)+1/(-f11));
    
    for(int i = 0; i < 4; i++) {
        dJ[i] = 2*(A[i][0]*u[0]+A[i][1]*u[1]+A[i][2]*u[2]+A[i][3]*u[3]- B[i]) + dBarrier[i]/t;
    }
    // matrix is symmetric so we only have to compute one side
    d2J[0][0] = 2*A[0][0]+1/t*(1/(f1*f1)+1/(f2*f2));
    d2J[1][1] = 2*A[1][1]+1/t*(1/(f3*f3)+1/(f4*f4));
    d2J[2][2] = 2*A[2][2]+1/t*(1/(f5*f5)+1/(f6*f6)+1/(f9*f9)+1/(f10*f10)+bias_ratio*bias_ratio/(f11*f11));
    d2J[3][3] = 2*A[3][3]+1/t*(1/(f7*f7)+1/(f8*f8)+1/(f9*f9)+bias_ratio*bias_ratio/(f10*f10)+1/(f11*f11));
    d2J[2][3] = 2*A[2][3]+1/t*(1/(f9*f9)-bias_ratio/(f10*f10)-bias_ratio/(f11*f11));
    d2J[3][2] = d2J[2][3];
}

double JCalc(double * u,double *torqueLimit) {
    double temp[4];
    for(int i = 0; i < 4; i++) {
        temp[i] = dotProduct(A[i],u,4); 
    }
    double f1 = (u[0]-torqueLimit[0]);
    double f2 = (-u[0]+torqueLimit[1]);
    double f3 = (u[1]-torqueLimit[2]);
    double f4 = (-u[1]+torqueLimit[3]);
    double f5 = (u[2]-torqueLimit[4]);
    double f6 = -u[2]+minDiff;
    double f7 = (u[3]-torqueLimit[5]);
    double f8 = -u[3]+minDiff;
    double f9 = (u[2]+u[3]-rearLimit);
    double f10 = (u[2]-bias_ratio*u[3]);
    double f11 = (-bias_ratio*u[2]+u[3]);

    double sum =dotProduct(u,temp,4)-2*dotProduct(B,u,4)+C;
    sum -= 1/t*(log(-f1)+log(-f2)+log(-f3)+log(-f4)+log(-f5)+log(-f6)+log(-f7)+log(-f8)+log(-f9)+log(-f11));
    return sum;
}



// checking if traction limits are beyond motor limts and correcting
void torqueLimitCorrection(double * torqueLimit) {
    const double bias_ratio = 4;
    double frontLeftConstraint,frontRightConstraint; 
    if(torqueLimit[0]>frontLimit){
        torqueLimit[0] = frontLimit;
    }
    if(torqueLimit[1]<-frontLimit) {
        torqueLimit[1] = -frontLimit;
    }

    if(torqueLimit[2]>frontLimit){
        torqueLimit[2] = frontLimit;
    }
    if(torqueLimit[3]<-frontLimit) {
        torqueLimit[3] = -frontLimit;
    }
}

// Function to solve A*x = b using Gaussian elimination with partial pivoting
short solve_system(double A[N][N], double b[N], double x[N]) {
    int i, j, k, max_row;
    double max_val, temp, factor;
    double aug[N][N+1];  // Augmented matrix [A|b]
    
    // Create augmented matrix
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][N] = b[i];
    }
    
    // Gaussian elimination with partial pivoting
    for (i = 0; i < N; i++) {
        // Find pivot row
        max_row = i;
        max_val = fabs(aug[i][i]);
        for (k = i + 1; k < N; k++) {
            if (fabs(aug[k][i]) > max_val) {
                max_val = fabs(aug[k][i]);
                max_row = k;
            }
        }
        
        // Swap rows if necessary
        if (max_row != i) {
            for (j = i; j <= N; j++) {
                temp = aug[i][j];
                aug[i][j] = aug[max_row][j];
                aug[max_row][j] = temp;
            }
        }
        
        // Check for singular matrix
        if (fabs(aug[i][i]) < 1e-10) {
            printf("Matrix is singular or nearly singular\n");
            return 0;
        }
        
        // Eliminate below
        for (k = i + 1; k < N; k++) {
            factor = aug[k][i] / aug[i][i];
            for (j = i; j <= N; j++) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }
    
    // Back substitution
    for (i = N - 1; i >= 0; i--) {
        x[i] = aug[i][N];
        for (j = i + 1; j < N; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    return 1;
}

void initialTorqueCorrection(double * torque, double * torqueLimits) {
    const int indArray[4] = {0,2,4,5};
    const double epsilon = 10;
    for(int i = 0; i < 4; i++) {
        if(torque[i]>torqueLimits[indArray[i]]) {
            torque[i] = torqueLimits[indArray[i]]-epsilon;
        }
    }
    
    if(torque[0]<torqueLimits[1]) {
        torque[0] = torqueLimits[1]+epsilon;
    }
    
    if(torque[1]<torqueLimits[3]) {
        torque[1]=torqueLimits[3]+epsilon;
    }

    if(torque[2]<minDiff) {
        torque[2] = minDiff+epsilon;
    }
    if(torque[3]<minDiff) {
        torque[3] = minDiff+epsilon;
    }

    if(-4*torque[2]+torque[3]>-1) {
        torque[2] = torque[3]/4+epsilon;
    }
    else if(torque[2]-4*torque[3]) {
        torque[3] = torque[2]/4+epsilon;
    }
}

void computeYawMoment(double * torque,double steering_angle) {
    double total_torque = torque[0]+torque[1]+torque[2]+torque[3];
    double yaw_Moment = 0;
    double leftWheelAngle = steering_angle;
    double rightWheelAngle = steering_angle;
    double B_v[4] = {(L_F*sin(leftWheelAngle)-halfWidth*cos(leftWheelAngle))/R_W,(L_F*sin(rightWheelAngle)+halfWidth*cos(rightWheelAngle))/R_W,-halfWidth/R_W,halfWidth/R_W};
    for(int i = 0; i <4; i++) {
        yaw_Moment += B_v[i]*torque[i];
    }
    printf("Total Torque:%f,Yaw Moment:%f\n",total_torque,yaw_Moment);
}



int main(int argc, char*argv[]) {
    double steering_angle,yawMomentRequest,total_torque;
    if(argc == 1) {
        steering_angle = 3.14/12;  // angle of steering wheel
        yawMomentRequest = 1000;
        total_torque = 547;
    }
    else if(argc == 4) {
        for (int i = 1; i < argc; i++) {
            char *endptr;
            double val = strtod(argv[i], &endptr);
            int unsafe = 0;
            if (*endptr != '\0') {
                // Not a valid double
                unsafe = 1;
                printf("Invalid number: %s\n", argv[i]);
            } 
            else {
                if(i == 2) {
                    yawMomentRequest = val;
                }
                else if(i == 1) {
                    steering_angle = val/180*3.141519;
                }
                else if(i == 3) {
                    total_torque = val;
                }
                printf("Argument %d as double: %f\n", i, val);
            }
            if(unsafe) {
                return 1;
            }
            
        }
    }
    else {
        printf("Usage: %s  (Steering Angle Degrees) (Yaw Moment Request Nm) (Total Torque Nm)\n",argv[0]);
        return 1;
    }


    clock_t start,end;
    // initalization of constants
    double elapsedTime;  
    
    start = clock();

    
    double normalForce[4] = {350*0.25,350*0.25,350*0.25,350*0.25};
    
    double torqueLimits[] = {200,-150,200,-150,150,150};
    initialize(steering_angle,normalForce,total_torque,yawMomentRequest);
    
    // torqueLimits
    unsigned int iteration = 0;
    torqueLimitCorrection(torqueLimits);
    double correctionTorqueF = yawMomentRequest*R_W/4.0/halfWidth;
    double torque[4] = {total_torque/4-correctionTorqueF,total_torque/4+correctionTorqueF,total_torque/4-correctionTorqueF,total_torque/4+correctionTorqueF};
    initialTorqueCorrection(torque,torqueLimits);
    double torque_initial[4] = {torque[0],torque[1],torque[2],torque[3]};
    // Gradient Descent 
    double step = 0.0005; // Step size
    double beta = 0.5;
    double alpha = 0.1;
    double nextTorque[4];
    double initalCost = JCalc(torque,torqueLimits);
    double JCost = initalCost;
    dJCalc(torque,torqueLimits);
    double currentCost = initalCost;
    double newCost;
    int broke = 0;
    int averageIter = 0;
    // Gradient Method
    
    iteration = 0;
    // while(iteration < 10000) {
    //     dJCalc(torque,torqueLimits);
    //     double gradDotDir = dotProduct(dJ,dJ,4);
    //     if(gradDotDir < 30000) {
    //         break;
    //     }
    //     step = 0.000000000001;
    //     for(short j = 0; j < 4; j++) {
    //         nextTorque[j] = torque[j]-step*dJ[j];
    //     }
    //     newCost = JCalc(nextTorque, torqueLimits);
        
    //     // Backtracking
    //     while(isnan(newCost) || newCost > currentCost + alpha * step * gradDotDir ) {
    //         step *= beta;
    //         for(short j = 0; j < 4; j++) {
    //             nextTorque[j] = torque[j] - step * dJ[j];
    //         }
    //         newCost = JCalc(nextTorque, torqueLimits);
    //         if(step<1e-100) {
    //             broke = 1;
    //             printf("BROKEN\n");
    //             break;
    //         }
    //         averageIter++;
    //     }
    //     if(broke) {
    //         break;
    //     }
    //     if(-4*nextTorque[2]+nextTorque[3]>-1) {
    //             printf("BROKEN\n");
    //             break;
    //     }
    //     currentCost = newCost;
        
    //     for(short j = 0; j < 4;j ++) {
    //         torque[j] = nextTorque[j];
    //     }

    //     iteration += 1;
    // }
    // averageIter /= iteration;
    
    // Newton Method
    double newtonStep[4];
    double lambda_squared;
    double tolerance =10000;
    iteration = 0;
    newCost = initalCost;
    dJd2JCalc(torque,torqueLimits);
    while(iteration < 12000 && solve_system(d2J,dJ,newtonStep)) {
        lambda_squared = dotProduct(dJ,newtonStep,4);
        if(lambda_squared/2 < tolerance) {
            break;
        }
        step = 0.0001;
        currentCost = newCost;
        
        for(int i = 0; i < 4; i++) {
            nextTorque[i] = torque[i] - step*newtonStep[i];
        }
        
        double newCost = JCalc(nextTorque, torqueLimits);
        // int backIter = 1;
        int broken = 0;
        while( isnan(newCost) || (newCost > currentCost+alpha*step*lambda_squared)){
            step = beta*step;
            for(int i = 0; i < 4; i++) {
                nextTorque[i] = torque[i] - step*newtonStep[i];
            }
            newCost = JCalc(nextTorque,torqueLimits);
            if(step<1e-100) {
                printf("BROKEN\n");
                broken = 1;
                break;
            }
            if(-4*nextTorque[3]+nextTorque[2]>-30) {
                printf("BROKEN\n");
                broken = 1;
                break;
            }
            // backIter++;
        }
        if(broken ) {
            break;
        }
        
        /* if(minBackIter > backIter) {
           minBackIter = backIter; 
        }
        else if(maxBackIter < backIter) {
            maxBackIter = backIter;
        }
        averageIter += backIter; */
        for(int i = 0; i < 4; i++) {
            torque[i] = nextTorque[i];
        }
        
        dJd2JCalc(torque,torqueLimits);
        iteration += 1;
    }
    int nIteration = iteration;
    
    iteration = 0;
    while(iteration < 7000) {
        dJCalc(torque,torqueLimits);
        double gradDotDir = dotProduct(dJ,dJ,4);
        if(gradDotDir < 500) {
            break;
        }
        step = 0.00000001;
        for(short j = 0; j < 4; j++) {
            nextTorque[j] = torque[j]-step*dJ[j];
        }
        newCost = JCalc(nextTorque, torqueLimits);
        
        // Backtracking
        while(isnan(newCost) || newCost > currentCost + alpha * step * gradDotDir ) {
            step *= beta;
            for(short j = 0; j < 4; j++) {
                nextTorque[j] = torque[j] - step * dJ[j];
            }
            newCost = JCalc(nextTorque, torqueLimits);
            if(step<1e-100) {
                broke = 1;
                printf("BROKENG\n");
                break;
            } 
            averageIter++;
        }
        if(-4*nextTorque[3]+nextTorque[2]>-1) {
                printf("BROKENG\n");
                broke = 1;
                break;
        }
        if(broke) {
            break;
        }
        currentCost = newCost;
        for(short j = 0; j < 4;j ++) {
            torque[j] = nextTorque[j];
        }

        iteration += 1;
    }
    end = clock();
    if(iteration != 0){ averageIter /= iteration;}
    
    elapsedTime = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Steering Angle %f,Yaw Request:%f, Total Torque:%f\n",steering_angle,yawMomentRequest,total_torque);
    printf("Torques: %f,%f,%f,%f\n",torque_initial[0],torque_initial[1],torque_initial[2],torque_initial[3]);
    printf("dJ: %f,%f,%f,%f\n",dJ[0],dJ[1],dJ[2],dJ[3]);
    printf("Total Elapsed time: %lf seconds\n\n", elapsedTime);
    printf("Initial Cost: %f\n",JCost);
    printf("Newton Iterations:%d\n",nIteration);
    printf("Grad Iterations:%d\n",iteration);
    printf("AverageIter:%d\n",averageIter);
    printf("Final Cost: %f\n",newCost);
    printf("Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);
    computeYawMoment(torque,steering_angle);
    return 0;
}
#include "Optimization.h"
#include "math.h"
#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
// global variables
double dJ[4] = {0,0,0,0}; // gradient of cost function
double d2J[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; // hessian of cost function 
double A[4][4];  // quadatric term
double B[4]; // linear term
double C; // const term
double dBarrier[4]; // barrier gradient 
double Weights_v[2][2] = {{50,0},{0,50}}; // should always be [W1 0; 0 W2]

// motor and diff limits
const double frontLimit = 15*13.5; // front motor
const double rearLimit = 130*3.23; // rear motor
const double minDiff = 1;
const double bias_ratio = 4; // ratio of more traction wheel to less traction wheel torque
const double lambda_max = 0.1;
// barrierfunction
const double t = 0.0000001; // strength of barrier


double dotProduct(double * a, double * b, int size) {
    double result = 0;

    for(int i = 0; i < size;i ++) {
        result += a[i]*b[i];
    }
    return result;
}


void initialize(double slip_ratio[4],double total_torque, double yaw_moment_request,struct tMyModel *mp) {
    double leftWheelAngle = Car.Tire[0].rot;
    double rightWheelAngle = Car.Tire[1].rot;
    double R_W = mp->wheelRadius;
    double L_F = mp->L_F;    
    double halfWidth = mp->halfWidth;

    double B_w[4];
    double B_v[2][4] = {{1,1,1,1},{(L_F*sin(leftWheelAngle)-halfWidth*cos(leftWheelAngle))/R_W,(L_F*sin(rightWheelAngle)+halfWidth*cos(rightWheelAngle))/R_W,-halfWidth/R_W,halfWidth/R_W}};
    
    for(int i = 0;i < 4; i++) {
        B_w[i] = 10/slip_ratio_max*slip_ratio[i];
    }
    double Weights_vB_v[2][4];
    double Weights_vc[2];
    Weights_vc[0] = Weights_v[0][0]*total_torque;
    Weights_vc[1] = Weights_v[1][1]*yaw_moment_request;

    for(int k = 0; k < 2; k++) { 
        for (int j = 0; j < 4; j++) {
            Weights_vB_v[k][j] = Weights_v[k][k]*B_v[k][j];
       }
    }
    // B
    for(int j = 0; j < 4;j++) {
        B[j] = Weights_vc[0]*Weights_vB_v[0][j]+Weights_vc[1]*Weights_vB_v[1][j];
    }
    // store A
    // initialize d2J
    for(int i = 0; i < 4; i++) {
        for(int j = i; j < 4; j++) {
            A[i][j] = Weights_vB_v[0][i]*Weights_vB_v[0][j]+Weights_vB_v[1][i]*Weights_vB_v[1][j]+B_w[i]*B_w[j];
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
    C = Weights_vc[0]*Weights_vc[0]+Weights_vc[1]*Weights_vc[1];
    return;
}

void dJCalc(double u[4], double torqueLimit[6]){ 
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
    dBarrier[2] = 1/(-f5)-1/(-f6)+1/(-f9)+1/(-f10)-bias_ratio/(-f11);
    dBarrier[3] = 1/(-f7)-1/(-f8)+1/(-f9)-bias_ratio/(-f10)+1/(-f11);
    // dBarrier[2] = 1/(-f5)+1/(-f9)+(1/(-f10)-bias_ratio/(-f11));
    // dBarrier[3] = 1/(-f7)+1/(-f9)-bias_ratio/(-f10)+1/(-f11);
    
    for(int i = 0; i < 4; i++) {
        dJ[i] = 2*(A[i][0]*u[0]+A[i][1]*u[1]+A[i][2]*u[2]+A[i][3]*u[3]- B[i]) + dBarrier[i]/t;
    }
}

void dJd2JCalc(double u[4], double torqueLimit[6]) {
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
    dBarrier[2] = 1/(-f5)-1/(-f6)+1/(-f9)+1/(-f10)-bias_ratio/(-f11);
    dBarrier[3] = 1/(-f7)-1/(-f8)+1/(-f9)-bias_ratio/(-f10)+1/(-f11);
    
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

double JCalc(double u[4],double torqueLimit[6]) {
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
    sum -= 1/t*(log(-f1)+log(-f2)+log(-f3)+log(-f4)+log(-f5)+log(-f6)+log(-f7)+log(-f8)+log(-f9)+log(-f10)+log(-f11));
    return sum;
}

// checking if traction limits are beyond motor limts and correcting
void torqueLimitCorrection(double torqueLimit[6]) {
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
short solve_system(double A[4][4], double b[4], double x[4]) {
    const int N = 4;
    int i, j, k, max_row;
    double max_val, temp, factor;
    double aug[4][4+1];  // Augmented matrix [A|b]
    
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

void initialTorqueCorrection(double torque[4], double torqueLimits[6]) {
    const int indArray[4] = {0,2,4,5};
    const double epsilon = 1;
    for(int i = 0; i < 4; i++) {
        if(torque[i]-torqueLimits[indArray[i]]>-epsilon) {
            torque[i] = torqueLimits[indArray[i]]-epsilon;
        }
    }
    
    if(torque[0]-torqueLimits[1]<epsilon) {
        torque[0] = torqueLimits[1]+epsilon;
    }
    
    if(torque[1]-torqueLimits[3]<epsilon) {
        torque[1]=torqueLimits[3]+epsilon;
    }

    if(torque[2]-minDiff<epsilon) {
        torque[2] = minDiff+epsilon;
    }
    if(torque[3]-minDiff<epsilon) {
        torque[3] = minDiff+epsilon;
    }

    if(-4*torque[2]+torque[3]>-epsilon) {
        torque[2] = torque[3]/4+epsilon;
    }
    else if(torque[2]-4*torque[3]>-epsilon) {
        torque[3] = torque[2]/4+epsilon;
    }
}

short correctionTorques(double torque[4],double steering_angle,double total_torque, double yawMomentRequest, struct tMyModel *mp) {
    double R_W = mp->wheelRadius;
    double L_F = mp->L_F;    
    double halfWidth = mp->halfWidth;
    
    double actual_total_torque = torque[0]+torque[1]+torque[2]+torque[3];
    double torque_error = actual_total_torque-total_torque; // difference between requested total torque and actual optimal torque
    double yaw_Moment = 0;
    double yaw_Moment_error = 0;
    double leftWheelAngle = Car.Tire[0].rot;
    double rightWheelAngle = Car.Tire[1].rot;
    double B_v[4] = {(L_F*sin(leftWheelAngle)-halfWidth*cos(leftWheelAngle))/R_W,(L_F*sin(rightWheelAngle)+halfWidth*cos(rightWheelAngle))/R_W,-halfWidth/R_W,halfWidth/R_W};
    
    if(torque_error > 0) {
        for(int i = 0; i < 4;i++) {
            torque[i] *= actual_total_torque/total_torque;
        }
    }
    
    for(int i = 0; i <4; i++) {
        yaw_Moment += B_v[i]*torque[i];
    }
    yaw_Moment_error = yaw_Moment-yawMomentRequest;

    if(yaw_Moment_error*yaw_Moment_error < 9000) {
        return 1;
    }
    else{
        return 0;
    }
}

// returns successfull run of torque optimization
short optimalTorque(double torque[4], struct tMyModel *mp, double yawMomentRequest, double total_torque, double slipRatio[4], double tractionLimits[6]) {
    double steering_angle = mp->steerAngle;
    double torqueLimits[6] = {tractionLimits[0],tractionLimits[1],tractionLimits[2],tractionLimits[3],tractionLimits[4],tractionLimits[5]}; // wheel torque on each tires
    initialize(slipRatio,total_torque,yawMomentRequest,mp);
    torqueLimitCorrection(torqueLimits);
    initialTorqueCorrection(torque,torqueLimits);

    double nextTorque[4]; // next iteration of torque 
    double initialCost = JCalc(torque,torqueLimits); // cost of initial torque
    if(isnan(initialCost)) {
        return 0;
    }
    double currentCost = initialCost;
    double iteration = 0;
    double newtonStep[4];
    double step = 0.0001; // step size
    double beta = 0.05;
    double alpha = 0.2; 
    double lambda_squared = 0;
    double tolerance = 1000;
    double newCost;
    int totalIter = 0; // total backiterations
    int maxTotalIter = 50000;

    // Newton Method
    dJd2JCalc(torque,torqueLimits);
    
    while(iteration < 12000 && solve_system(d2J,dJ,newtonStep)) {
        lambda_squared = dotProduct(dJ,newtonStep,4);
        if(lambda_squared/2 < tolerance) {
            break;
        }
        step = 1;
        currentCost = newCost;
        
        for(int i = 0; i < 4; i++) {
            nextTorque[i] = torque[i] - step*newtonStep[i];
        }
        
        double newCost = JCalc(nextTorque, torqueLimits);
        while( isnan(newCost) || (newCost > currentCost+alpha*step*lambda_squared)){
            step = beta*step;
            for(int i = 0; i < 4; i++) {
                nextTorque[i] = torque[i] - step*newtonStep[i];
            }
            newCost = JCalc(nextTorque,torqueLimits);
            if(step<1e-100) {
                return 0;
            }
            totalIter++;
        }

        for(int i = 0; i < 4; i++) {
            torque[i] = nextTorque[i];
        }
        
        if(totalIter>maxTotalIter){
            break;
        }
        
        dJd2JCalc(torque,torqueLimits);
        iteration += 1;
    }

    return correctionTorques(torque,steering_angle,total_torque,yawMomentRequest,mp); 
}
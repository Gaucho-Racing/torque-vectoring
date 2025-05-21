#include<stdio.h>
#include<math.h>
#include<windows.h>
#define R_W 1
#define L_F 1
#define L_S 1
#define N 4

double W_v[2] = {1000, 1000};
double B_w[4];
double B_v[2][4] = {{1,1,1,1},{1,1,-L_S/R_W,L_S/R_W}};
double dJ[4] = {0,0,0,0}; // gradient of cost function
double d2J[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; // hessian of cost function 
double A[4][4];  // quadatric term
double B[4]; // linear term
double c; // const term
double dBarrier[4]; // barrier gradient 
const double t = 16;
double newtonStep[4];
const double correctionFactor = 1;
// convention: 0 -> FL,1 -> FR, 2 -> RL, 3 -> RR
// convex surface can be represented as u^TAu+2Bu+c, c can be ignored
// barrier terms allow the solution to remain in the allowable region

// Initialize
void initialize(double steering_angle, double normal_force[4]) {
    B_v[1][0] = (L_F*sin(steering_angle)-L_S*cos(steering_angle))/R_W;
    B_v[1][1] = (L_F*sin(steering_angle)+L_S*cos(steering_angle))/R_W;
    for(short i = 0; i < 4;i ++) {
        B_w[i] = 1/normal_force[i]/R_W;
    }
    // printf("Intialized\n");
    return;
}

double dotProduct(double * a, double * b, int size) {
    double result = 0;

    for(int i = 0; i < size;i ++) {
        result += a[i]*b[i];
    }
    

    return result;
}

// compute A and B
void constCalc( double total_torque, double yaw_moment_request) {
    double W_vB_v[4];
    double W_vc = W_v[0]*total_torque+W_v[1]*yaw_moment_request;

    for (int j = 0; j < 4; j++) {
        W_vB_v[j] = 0;
        for (int k = 0; k < 2; k++) {
            W_vB_v[j] += W_v[k] * B_v[k][j];
        }
        // W_vB_v[j] *= 2;
        B[j] = W_vB_v[j] * W_vc;
    }
    // printf("(");
    // for(int i = 0; i < 4;i++) {
    //     printf("%f ",W_vB_v[i]);
    // }
    // printf(")");

    // store A
    // initialize d2J
    for(int i = 0; i < 4; i++) {
        for(int j = i; j < 4; j++) {
            A[i][j] = W_vB_v[i]*W_vB_v[j]+B_w[i]*B_w[j];
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
    c = W_vc*W_vc;
    return;
}

// 4 torque input, 3 motor where RR and RL are added together
// torque limit: [UB_FL,LB_FL,UB_FR,LB_FR,UB_RL,UB_RR ]
void dbarrierCalc(double * u , double * torqueLimit){
    const double frontLimit = 15*13.5;
    const double rearLimit = 230*3.23;
    const double bias_ratio = 4;
    const double minDiff = 10;
    double frontLeftConstraint,frontRightConstraint; 
    if(torqueLimit[0]<frontLimit){
        frontLeftConstraint = torqueLimit[0];
    }
    else {
        frontLeftConstraint = frontLimit;
    }
    if(torqueLimit[0]<frontLimit){
       frontRightConstraint = torqueLimit[2];
    }
    else {
       frontRightConstraint = frontLimit;
    }
    
    double f1 = (u[0]-frontLeftConstraint);
    double f2 = (-u[0]+torqueLimit[1]);
    double f3 = (u[1]-frontRightConstraint);
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
    dBarrier[2] = 1/(-f5)-correctionFactor/(-f6)+1/(-f9)+1/(-f10)-bias_ratio/(-f11);
    dBarrier[3] = 1/(-f7)-correctionFactor/(-f8)+1/(-f9)-bias_ratio/(-f10)+1/(-f11);
}


// gradient of J, cost function including barrier terms
void dJCalc(double * torque,double * limits) {
    dbarrierCalc(torque,limits);
    for(int i = 0; i < 4; i++) {
        dJ[i] = 2*(A[i][0]*torque[0]+A[i][1]*torque[1]+A[i][2]*torque[2]+A[i][3]*torque[3]- B[i]) + dBarrier[i]/t;
        // dJ[i] = 2*(dotProduct(A[i],torque,4)-B[i])+dBarrier[i]/t;
    }
    return;
}

// limits is FL/FR/R motor torque limits
// torque limit: [UB_FL,LB_FL,UB_FR,LB_FR,UB_RL,UB_RR ]
void d2JCalc(double * u,double * torqueLimit) {
    const double frontLimit = 15*13.5;
    const double rearLimit = 230*3.23;
    const double bias_ratio = 4;
    const double minDiff = 10;  
    double frontLeftConstraint,frontRightConstraint; 
    if(torqueLimit[0]<frontLimit){
        frontLeftConstraint = torqueLimit[0];
    }
    else {
        frontLeftConstraint = frontLimit;
    }
    if(torqueLimit[0]<frontLimit){
       frontRightConstraint = torqueLimit[2];
    }
    else {
       frontRightConstraint = frontLimit;
    }
    double f1 = (u[0]-frontLeftConstraint);
    double f2 = (-u[0]+torqueLimit[1]);
    double f3 = (u[1]-frontRightConstraint);
    double f4 = (-u[1]+torqueLimit[3]);
    double f5 = (u[2]-torqueLimit[4]);
    double f6 = -u[2]+minDiff;
    double f7 = (u[3]-torqueLimit[5]);
    double f8 = -u[3]+minDiff;
    double f9 = (u[2]+u[3]-rearLimit);
    double f10 = (u[2]-4*u[3]);
    double f11 = (-4*u[2]+u[3]);
    
    // matrix is symmetric so we only have to compute one side
    d2J[0][0] = 2*A[0][0]+1/t*(1/(f1*f1)+1/(f2*f2));
    d2J[1][1] = 2*A[1][1]+1/t*(1/(f3*f3)+1/(f4*f4));
    d2J[2][2] = 2*A[2][2]+1/t*(1/(f5*f5)+correctionFactor/(f6*f6)+1/(f9*f9)+1/(f10*f10)+16/(f11*f11));
    d2J[3][3] = 2*A[3][3]+1/t*(1/(f7*f7)+correctionFactor/(f8*f8)+1/(f9*f9)+16/(f10*f10)+1/(f11*f11));
    d2J[2][3] = 2*A[2][3]+1/t*(1/(f9*f9)-4/(f10*f10)-4/(f11*f11));
    d2J[3][2] = d2J[2][3];
    
}

double JCalc(double * u,double *torqueLimit) {
    double temp[4];
    for(int i = 0; i < 4; i++) {
        temp[i] = dotProduct(A[i],u,4); 
    }
    const double frontLimit = 15*13.5;
    const double rearLimit = 230*3.23;
    const double bias_ratio = 4;
    const double minDiff = 10;  
    double frontLeftConstraint,frontRightConstraint; 
    if(torqueLimit[0]<frontLimit){
        frontLeftConstraint = torqueLimit[0];
    }
    else {
        frontLeftConstraint = frontLimit;
    }
    if(torqueLimit[0]<frontLimit){
       frontRightConstraint = torqueLimit[2];
    }
    else {
       frontRightConstraint = frontLimit;
    }
    double f1 = (u[0]-frontLeftConstraint);
    double f2 = (-u[0]+torqueLimit[1]);
    double f3 = (u[1]-frontRightConstraint);
    double f4 = (-u[1]+torqueLimit[3]);
    double f5 = (u[2]-torqueLimit[4]);
    double f6 = -u[2]+minDiff;
    double f7 = (u[3]-torqueLimit[5]);
    double f8 = -u[3]+minDiff;
    double f9 = (u[2]+u[3]-rearLimit);
    double f10 = (u[2]-4*u[3]);
    double f11 = (-4*u[2]+u[3]);
    double sum =dotProduct(u,temp,4)-2*dotProduct(B,u,4)+c;
    sum -= 1/t*(log10(-f1)+log10(-f2)+log10(-f3)+log10(-f4)+log10(-f5)+correctionFactor*(-f6)+log10(-f7)+correctionFactor*(-f8)+log10(-f9)+log10(-f11));
    return sum;
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


int main() {
    LARGE_INTEGER frequency;
    LARGE_INTEGER start, end;
    double elapsedTime;

    // Get the frequency of the high-resolution performance counter
    QueryPerformanceFrequency(&frequency);

    // Get the starting time
    QueryPerformanceCounter(&start);
    double normal[4] = {350/4,350/4,350/4,350/4};
    double steering_angle = 3.14/6;
    initialize(steering_angle,normal);
    double total_torque = 847;
    constCalc(total_torque,3);
    double longVel = 30;
    double wheelBase = 3;
    double K_U = 0.025; 
	double yaw_ref = longVel/(wheelBase*(1+K_U*longVel*longVel))*steering_angle;
    double correctionTorqueF = 50*yaw_ref*R_W/2.0/1.6;


    double torqueLimits[] = {100,-100,100,-100,100,100};
    double torque[4] = {total_torque/4,total_torque/4,total_torque/4,total_torque/4};
    double torque_initial[4]= {total_torque/4,total_torque/4,total_torque/4,total_torque/4}; 
    double rate = 0.0000001;
    int r = 0;
    printf("Initial Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);
    double initialJ = JCalc(torque_initial,torqueLimits);
    dJCalc(torque,torqueLimits);
    // printf("Magn: %f\n",dotProduct(dJ,dJ,4));
    printf("dJ: %f,%f,%f,%f\n",dJ[0],dJ[1],dJ[2],dJ[3]);
    double step = 1.0;  // Start with a larger step
    double beta = 0.5;  // Same as in Newton's method
    double alpha = 0.0001;  // Line search parameter
    double currentCost = JCalc(torque, torqueLimits);
    double tempTorque[4];

    dJCalc(torque, torqueLimits);
    while(r < 10000 && dotProduct(dJ, dJ, 4) > 30000) {
        // Try step
        for(short j = 0; j < 4; j++) {
            tempTorque[j] = torque[j] - step * dJ[j];
        }
        
        // Check if cost decreases sufficiently
        double newCost = JCalc(tempTorque, torqueLimits);
        double gradDotDir = -dotProduct(dJ, dJ, 4);  // Dot product of gradient and -gradient
        
        // Backtracking
        while(newCost > currentCost + alpha * step * gradDotDir || isnan(newCost)) {
            step *= beta;
            for(short j = 0; j < 4; j++) {
                tempTorque[j] = torque[j] - step * dJ[j];
            }
            newCost = JCalc(tempTorque, torqueLimits);
            if(step < 1e-10) break;  // Prevent infinite loop
        }
        
        // Update torques
        for(short j = 0; j < 4; j++) {
            torque[j] = tempTorque[j];
        }
        
        currentCost = newCost;
        dJCalc(torque, torqueLimits);
        r += 1;
    }
    // while(r < 1000 && dotProduct(dJ,dJ,4)>30000) {
    //     dJCalc(torque,torqueLimits);
    //     for(short j = 0; j < 4; j++) {
    //         torque[j] -= rate*dJ[j];
    //     }
    //     r += 1;
    // }
    // QueryPerformanceCounter(&end);

    // Calculate the elapsed time in seconds
    // elapsedTime = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    // printf("Medium Elapsed time: %lf seconds\n\n", elapsedTime);
    // printf("Intial Cost: %f\n",JCalc(torque_initial,torqueLimits));
    // printf("Magn: %f\n",dotProduct(dJ,dJ,4));
    // printf("Magnitude: %f\n",dotProduct(dJ,dJ,4));
    // printf("Intermediate Cost: %f\nIterations:%i\n",JCalc(torque,torqueLimits),r);
    // printf("Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);

    // printf("Newtons's Method Started\n");
    
    d2JCalc(torque,torqueLimits);
    dJCalc(torque,torqueLimits);
    double tolerance = 10000;
    double temp[4];
    double lambda_squared = 0;
    alpha = 0.01;
    beta = 0.5;
    step = 0.000001;
    double initialCost=JCalc(torque,torqueLimits);
    int iteration = 0;
    double tempCost = initialCost;
    while(solve_system(d2J,dJ,newtonStep) && iteration < 5000) {
        lambda_squared = dotProduct(dJ,newtonStep,4);
        if(lambda_squared/2 < tolerance) {
            break;
        }
        step = 0.0001;
        initialCost = tempCost;
        for(int i = 0; i < 4; i++) {
            temp[i] = torque[i] - step*newtonStep[i];
        }
        double tempCost = JCalc(temp,torqueLimits);
        while( (tempCost > initialCost+alpha*step*lambda_squared) || isnan(tempCost)){
            step = beta*step;
            for(int i = 0; i < 4; i++) {
                temp[i] = torque[i] - step*newtonStep[i];
            }
            tempCost = JCalc(temp,torqueLimits);
            if(temp[2]<10) {
                printf("BROKEN");
            }
        }
        for(int i = 0; i < 4; i++) {
            torque[i] = temp[i];
        }
        
        iteration += 1;
    }
    
    // printf("Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);
    
    QueryPerformanceCounter(&end);
    elapsedTime = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
    printf("Total Elapsed time: %lf seconds\n\n", elapsedTime);
    printf("Lambda^2: %f\n", lambda_squared);
    printf("Intial Cost:%f\n",initialJ);
    printf("Final Cost: %f\nIterations:%i\n",JCalc(torque,torqueLimits),iteration);
    printf("Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);


    return 0;
}

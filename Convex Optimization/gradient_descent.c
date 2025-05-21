#include<stdio.h>
#include<math.h>
#include<windows.h>
#define R_W 1
#define L_F 0.7
#define L_S 0.3

double W_v[2] = {1000, 1000};
double B_w[4];
double u[4];
double B_v[2][4] = {{1,1,1,1},{1,1,-L_S/R_W,L_S/R_W}};
double dJ[4] = {0,0,0,0}; // gradient of cost function
double A[4][4];  // quadatric term
double B[4]; // linear term
double c; // const term
double dBarrier[4]; // barrier gradient 
const double t = 8;
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

    for(int i = 0; i < 4; i++) {
        for(int j = i; j < 4; j++) {
            A[i][j] = W_vB_v[i]*W_vB_v[j]+B_w[i]*B_w[j];
            // adding identity matrix
            if(i == j) {
                A[i][j] += 1;
            }
            // symmetric matrix
            else{
                A[j][i] = A[i][j];
            }
        }
    }
    c = W_vc*W_vc;
    return;
}

// 4 torque input, 3 motor where RR and RL are added together
void dbarrierCalc(double * u , double * torqueLimit){
    dBarrier[0] = 1/(-u[0]+torqueLimit[0])-1/(u[0]+torqueLimit[0]);
    dBarrier[1] = 1/(-u[1]+torqueLimit[1])-1/(u[1]+torqueLimit[1]);
    dBarrier[2] = 1/(u[2]+u[3]-torqueLimit[2]);
    dBarrier[3] = 1/(u[2]+u[3]-torqueLimit[2]);
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

double JCalc(double * torque) {
    double temp[4];
    for(int i = 0; i < 4; i++) {
        temp[i] = dotProduct(A[i],torque,4); 
    }
    return dotProduct(torque,temp,4)-2*dotProduct(B,torque,4)+c;
}

// Gradient descent with backtracking line search
int gradientDescent(double *torque, double *torqueLimits, int maxIter, double tol) {
    double torque_new[4];
    double alpha = 0.000001;       // Initial step size
    double beta = 0.5;        // Step size reduction factor
    double c1 = 1e-4;         // Sufficient decrease parameter (Armijo condition)
    int iter = 0;
    double gradientNorm;
    double cost, new_cost;
    double tolerance = tol*tol;
    // Calculate initial gradient and cost
    dJCalc(torque, torqueLimits);
    cost = JCalc(torque);
    gradientNorm = (dotProduct(dJ, dJ, 4));
    
    // printf("Initial cost: %f, Initial gradient norm: %f\n", cost, gradientNorm);
    
    while (iter < maxIter && gradientNorm > tolerance) {
        // Store current gradient magnitude for comparison
        double gradientMagnitude = dotProduct(dJ, dJ, 4);
        
        // Backtracking line search
        double step = alpha;
        int backtrack_iter = 0;
        int feasible = 0;
        
        while (!feasible && backtrack_iter < 20) {
            // Try this step size
            for (int i = 0; i < 4; i++) {
                torque_new[i] = torque[i] - step * dJ[i];
            }
            
            // Check constraints (ensure feasibility)
            if (torque_new[0] > -torqueLimits[0] && torque_new[0] < torqueLimits[0] &&
                torque_new[1] > -torqueLimits[1] && torque_new[1] < torqueLimits[1] &&
                torque_new[2] + torque_new[3] < torqueLimits[2]) {
                
                // Point is feasible, now check sufficient decrease (Armijo condition)
                new_cost = JCalc(torque_new);
                double decrease = -step * c1 * gradientMagnitude;
                
                if (new_cost <= cost + decrease || backtrack_iter >= 15) {
                    feasible = 1; // Accept this step
                } else {
                    step *= beta; // Reduce step size
                }
            } else {
                step *= beta; // Reduce step size to maintain feasibility
            }
            
            backtrack_iter++;
        }
        
        // Check if we found a feasible step
        if (!feasible) {
            printf("Warning: Could not find feasible step, terminating.\n");
            break;
        }
        
        // Update solution
        for (int i = 0; i < 4; i++) {
            torque[i] = torque_new[i];
        }
        
        // Recalculate gradient and cost
        dJCalc(torque, torqueLimits);
        cost = new_cost;
        
        // Check for NaN
        int has_nan = 0;
        for (int i = 0; i < 4; i++) {
            if (isnan(dJ[i]) || isnan(torque[i])) {
                has_nan = 1;
                break;
            }
        }
        
        if (has_nan) {
            printf("Warning: NaN detected, terminating.\n");
            break;
        }
        
        gradientNorm = sqrt(dotProduct(dJ, dJ, 4));
        
        // Print progress occasionally
        // if (iter % 100 == 0) {
        //     printf("Iteration %d: cost = %f, gradient norm = %f\n", iter, cost, gradientNorm);
        // }
        
        iter++;
    }
    
    printf("Final cost: %f\n", cost);
    printf("Final gradient norm: %f\n", gradientNorm);
    printf("Iterations: %d\n", iter);
    
    return iter;
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
    double steering_angle = 3.14/12;
    initialize(steering_angle,normal);
    double total_torque = 647;
    constCalc(total_torque,3);
    double longVel = 30;
    double wheelBase = 3;
    double K_U = 0.025; 
	double yaw_ref = longVel/(wheelBase*(1+K_U*longVel*longVel))*steering_angle;
    double correctionTorqueF = 50*yaw_ref*R_W/2.0/1.6;

    double torqueLimits[] = {15*13.5,15*13.5,130*3.4};
    double torque[4] = {total_torque/4-correctionTorqueF,total_torque/4+correctionTorqueF,total_torque/4-correctionTorqueF,total_torque/4+correctionTorqueF};
    double torque_initial[4];
    for (int i = 0; i < 4; i++) {
        torque_initial[i] = torque[i];
    }
    double rate = 0.00000001;
    int r = 0;
    dJCalc(torque,torqueLimits);
    // printf("Magn: %f\n",dotProduct(dJ,dJ,4));
    printf("dJ: %f,%f,%f,%f\n",dJ[0],dJ[1],dJ[2],dJ[3]);
    double initalJ = JCalc(torque_initial);
    double maxIterations = 10000;
    double tolerance = 3;

    int iterations = gradientDescent(torque, torqueLimits, maxIterations, tolerance);
    QueryPerformanceCounter(&end);

    // Calculate the elapsed time in seconds
    elapsedTime = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    printf("Elapsed time: %f seconds\n\n", elapsedTime);
    printf("Intial Cost: %f\n",JCalc(torque_initial));
    printf("Magn: %f\n",dotProduct(dJ,dJ,4));
    // printf("Magnitude: %f\n",dotProduct(dJ,dJ,4));
    printf("Final Cost: %f\nIterations:%i\n",JCalc(torque),iterations);
    printf("Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);

    return 0;
}
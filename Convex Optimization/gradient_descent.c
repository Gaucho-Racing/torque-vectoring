#include<stdio.h>
#include<math.h>
#include<windows.h>
#define R_W 1
#define L_F 1
#define L_S 1

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
  
    double torqueLimits[] = {15*13.5,15*13.5,130*3.4};
    double torque[4] = {total_torque/4,total_torque/4,total_torque/4,total_torque/4};
    double torque_initial[4]= {total_torque/4,total_torque/4,total_torque/4,total_torque/4}; 
    double rate = 0.00000001;
    int r = 0;
    dJCalc(torque,torqueLimits);
    // printf("Magn: %f\n",dotProduct(dJ,dJ,4));
    // printf("dJ: %f,%f,%f,%f\n",dJ[0],dJ[1],dJ[2],dJ[3]);
    while(r < 20000 && dotProduct(dJ,dJ,4)>175000) {
        dJCalc(torque,torqueLimits);
        for(short j = 0; j < 4; j++) {
            torque[j] -= rate*dJ[j];
        }
        r += 1;
    }
    QueryPerformanceCounter(&end);

    // Calculate the elapsed time in seconds
    elapsedTime = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    printf("Elapsed time: %f seconds\n\n", elapsedTime);
    printf("Intial Cost: %f\n",JCalc(torque_initial));
    printf("Magn: %f\n",dotProduct(dJ,dJ,4));
    // printf("Magnitude: %f\n",dotProduct(dJ,dJ,4));
    printf("Final Cost: %f\nIterations:%i\n",JCalc(torque),r);
    printf("Torques: %f,%f,%f,%f\n",torque[0],torque[1],torque[2],torque[3]);

    return 0;
}
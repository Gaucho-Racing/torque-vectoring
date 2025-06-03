#ifndef _OPTIMI_H__
#define _OPTIMI_H__
#define slip_ratio_max (0.2)

struct tMyModel {
    double steerAngle;
	double longVel;
	double yawRate;
	double wheelBase;
	double halfWidth; // Vehicle Half Width
	double wheelRadius; // mean wheel radius
    double L_F; // distance from cog to front axle
	double K_U;	// under steer gradient
	// proportional and integral constant error for yaw rate ref
	double K_P;
	double K_I;
	double yawErrorI;
	double prevTime;
	double FL_Ratio;
	double FR_Ratio;
	double R_Ratio;
	double T_max;
	double K_PSR;
	double K_ISR;
	double slipRatioFLI;
	double slipRatioFRI;
};




double dotProduct(double * a, double * b, int size);

void initialize(double slip_ratio[4],double total_torque, double yaw_moment_request,struct tMyModel *mp);

void dJCalc(double u[4], double torqueLimit[6]); 

void dJd2JCalc(double u[4], double torqueLimit[6]);

double JCalc(double u[4],double torqueLimit[6]);

short solve_system(double A[4][4], double b[4], double x[4]);

void initialTorqueCorrection(double torque[4], double torqueLimits[6]);

short correctionTorques(double torque[4],double steering_angle,double total_torque, double yawMomentRequest, struct tMyModel *mp); 

void torqueLimitCorrection(double torqueLimit[6]);

short optimalTorque(double torque[4], struct tMyModel *mp, double yawMomentRequest, double total_torque, double slipRatio[4], double tractionLimits[6]);
#endif
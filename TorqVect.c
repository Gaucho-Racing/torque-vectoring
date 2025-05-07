/*
******************************************************************************
**  CarMaker - Version 13.0.1
**  Vehicle Dynamics Simulation Toolkit
**
**  Copyright (C)   IPG Automotive GmbH
**                  Bannwaldallee 60             Phone  +49.721.98520.0
**                  76185 Karlsruhe              Fax    +49.721.98520.99
**                  Germany                      WWW    www.ipg-automotive.com
******************************************************************************
**
** Simple powertrain control Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    PTControl_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.Control";
static const char ThisModelKind[]  = "TorqVect";
static const int  ThisVersionId    = 1;


struct tMyModel {
    double steerAngle;
	double longVel;
	double yawRate;
	double wheelBase;
	double halfWidth; // Vehicle Half Width
	double wheelRadius; // mean wheel radius
	double K_U;	// under steer gradient
	// proportional and integral constant error for yaw rate ref
	double K_P;
	double K_I;
	double yawErrorI;
	double prevTime;
	double FL_Ratio;
	double FR_Ratio;
	double R_Ratio;
	tDDictEntry * pK_U; 
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    /* static struct tMyModel MyModel_Dummy = { GBKind_NoGearBox };
    if (park)
	mp = &MyModel_Dummy; */

    /* Define here dict entries for dynamically allocated variables. */
}


static void
MyModel_DeclQuants (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
	// mp->first = 0;
    if (mp == NULL) {
	/* Define here dict entries for non-dynamically allocated (static) variables. */
		// double FrontLeftPos,RearLeftPos;
		// // e = DDefDouble(NULL,"Car.CRL.tx", "m", &FrontLeftPos, DVA_None);
		// // Log("Value:%f",FrontLeftPos);
		
		// f = DDefDouble(NULL, "Tr.CRL.C.t_1.x", "m", &RearLeftPos, DVA_None);
		// mp->wheelBase = FrontLeftPos-RearLeftPos;
		// if(e == NULL) {
		// 	LogErrF(EC_Init,"e, isn't a valid UAQ");
		// }
		// if(f == NULL) {
		// 	LogErrF(EC_Init,"f, isn't a valid UAQ");
		// }
		// Log("WheelBase: %f",mp->wheelBase);
    }
	else {
		// Calculating WheelBase using UAQ
		tDDictEntry* e= DDictGetEntry("Car.CFL.C.t_0.x");
		tDDictEntry* f= DDictGetEntry("Car.CRL.C.t_0.x");
		if(e!=NULL && f!=NULL){
			mp->wheelBase = e->GetFunc(e->Var) - f->GetFunc(f->Var);
			Log("Vehicle Wheel Base:%f\n",mp->wheelBase);
		}
		else {
			mp->wheelBase = 1;
		}
		// Calculating Vehicle Width using UAQ
		e= DDictGetEntry("Car.CFL.C.t_0.y");
		if(e!=NULL){
			mp->halfWidth = e->GetFunc(e->Var);
			Log("Vehicle Half Width:%f\n",mp->halfWidth);
		}
		else {
			mp->halfWidth = 0.5;
		}
		// tDDictEntry* pK_U = DDictGetEntry("Understeer Gradient");
		// if(pK_U == NULL) {
		// 	LogErrF(EC_General, "Understeer Gradient Access Point Failled");
		// }
		// mp->pK_U = pK_U;
		mp->K_U = 0.025;
		Log("Understeer Gradient:%f\n",mp->K_U);
		mp->K_P =95;
		mp->K_I = 1;
		mp->yawErrorI = 0;
		mp->prevTime = 0;
		MyModel_DeclQuants_dyn (mp, 0);
    }
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}

/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (struct tInfos  	    *Inf,
    	     struct tPTControlCfgIF *CfgIF,
             const char     	    *KindKey)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if (CfgIF->PTKind != PTKind_BEV) {
	LogErrF (EC_Init, "%s: supports only Electric powertrain", MsgPre);
	return NULL;
    }
	
	// LogErrF (EC_Init, "%s: supports only Electric powertrain", MsgPre);

	if (CfgIF->nMotor != 3) {
		LogErrF(EC_Init, "supports only 3 motors, %d motors",CfgIF->nMotor);
		return NULL;
	}

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTControl, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));
	mp->wheelRadius = CfgIF->WheelRadius;
	mp->FL_Ratio = CfgIF->Motor[0].Ratio;
	mp->FR_Ratio = CfgIF->Motor[1].Ratio;
	mp->R_Ratio = 3.4;

    /* get CfgIF parameters */
    if (PTControl_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	return NULL;

    /* CfgIF -> Model */
    // mp->GBKind = CfgIF->GearBox[0].GBKind;

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->StartEngineWithSST) {
	LogErrF (EC_Init, "%s: no support for using SST", MsgPre);
	return NULL;
    }

    return mp;
}


/* static int
PedalsReady2Start (struct tMyModel *mp, tPTControlIF *IF)
{
    int Ready = 1, TransmOK = 0;

    switch (mp->GBKind) {
	case (GBKind_NoGearBox):
	    TransmOK = 1;
	    break;
	case (GBKind_Manual):
	    if (IF->GearNoTrg==0 && IF->Clutch >= 0.9)
		TransmOK = 1;
	    break;
	case (GBKind_AutoWithManual):
	case (GBKind_AutoNoManual):
	    if (IF->SelectorCtrl == SelectorCtrl_N)
		TransmOK = 1;
	    break;
    }
    if (TransmOK && IF->Brake >= 0.5)
	Ready = 1;

    return Ready;
} */


static int
MyModel_Calc (void *MP, tPTControlIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
	mp->steerAngle = Steering.IF.Ang;
	mp->longVel = Car.ConBdy1.v_1[0];
	// mp->totalVel = Car.ConBdy1.vHori;
	mp->yawRate = Car.YawRate;
		
    /* Operation Error */
    IF->OperationError = No_WarnError;
	
    /* Strategy mode */
    IF->StrategyMode = Mode_ElecDrive;
	for(int i = 0; i <3; i++) {
		IF->MotorOut[i].Load = -99999;
		IF->MotorOut[i].rotv_trg = -99999;
	}
	/* IF->BattHVOut.Temp_trg = -99999;
	IF->BattLVOut.Temp_trg = -99999;
	IF->BattLVOut.MassFlowCool_trg = -99999;
	IF->BattHVOut.MassFlowCool_trg = -99999;
	IF->PwrSupplyOut.Pwr_HV1toLV_trg = -99999; */
	// IF->EngineOut.FuelCutOff = -99999; 
	// IF->EngineOut.Load = -99999;
	// IF->EngineOut.Trq_trg = -99999;
	// IF->EngineOut.rotv_trg = -99999;


    /* Actual & Target Operation State Handling */
    switch (IF->OperationState) {
	case (OperState_Absent):
	    /* Absent -> PowerOff */
	    if (IF->Key >= KeyPos_KeyIn_PowerOff)
		IF->OperationState = OperState_PowerOff;
	    break;

	case (OperState_PowerOff):
	    /* PowerOff -> Absent */
	    if (IF->Key == KeyPos_KeyOut) {
		IF->OperationState = OperState_Absent;
		goto OutOfOperState;
	    }

	    /* PowerOff -> PowerAcc */
	    if (IF->Key >= KeyPos_KeyIn_PowerAcc) {
		IF->OperationState = OperState_PowerAccessory;
		goto OutOfOperState;
	    }
	    break;

	case (OperState_PowerAccessory):
	    /* PowerAcc -> PowerOn */
	    if (IF->Key >= KeyPos_KeyIn_PowerOn) {
		IF->OperationState = OperState_PowerOn;
		goto OutOfOperState;
	    }
	    break;

	case (OperState_PowerOn):
	    /* PowerOn -> PowerOff */
	    if (IF->Key<=KeyPos_KeyIn_PowerOff) {
		IF->MotorOut[0].Trq_trg = 0.0;
		IF->MotorOut[1].Trq_trg = 0.0;
		IF->MotorOut[2].Trq_trg = 0.0;
		IF->OperationState = OperState_PowerOff;
		IF->Ignition = 0;
		goto OutOfOperState;
	    }
		
	    /* PowerOn -> Driving */
	    if (IF->Key == KeyPos_KeyIn_Starter) {
			for(int i = 0; i <3; i++) {
				IF->MotorOut[i].Load = -99999;
				IF->MotorOut[i].rotv_trg = -99999;
			}
			IF->Ignition = 1;
			IF->MotorOut[0].Trq_trg = IF->MotorIn[0].TrqMot_max;
			IF->MotorOut[1].Trq_trg = IF->MotorIn[1].TrqMot_max;
			IF->MotorOut[2].Trq_trg = IF->MotorIn[2].TrqMot_max;
			// mp->K_U = mp->pK_U->GetFunc(mp->pK_U->Var);
			Log("Understeer Gradient: %f\n",mp->K_U);
			IF->OperationState = OperState_Driving;
		goto OutOfOperState;
	    } 
		else {
			IF->Ignition = 0;
			IF->MotorOut[0].Trq_trg = 0.0;
			IF->MotorOut[1].Trq_trg = 0.0;
			IF->MotorOut[2].Trq_trg = 0.0;
	    }
	    break;

	case (OperState_Driving):
	    // /* Driving -> PowerOn */
	    // if (!IF->EngineIn.Engine_on) {
		// IF->EngineOut.Load    = 0.0;
		// IF->EngineOut.set_ISC = 0;
		// IF->OperationState = OperState_PowerOn;
		// goto OutOfOperState;
	    // }
		
	    /* Driving -> PowerOff */
	    if (IF->Key<=KeyPos_KeyIn_PowerOff) {
			for(int i =0;i < 3;i++ ) {
				IF->MotorOut[i].Trq_trg = 0;
			}
			IF->OperationState    = OperState_PowerOff;
			goto OutOfOperState;
	    }

		// Checking if Torque Vectoring should be turned on
		if(/*IF->UserSignal[4]*/1) {
			
			
			// mp->K_U = mp->pK_U->GetFunc(mp->pK_U->Var);
			
			double yaw_ref = mp->longVel/(mp->wheelBase*(1+mp->K_U*mp->longVel*mp->longVel))*mp->steerAngle;
			double yawError = yaw_ref- mp->yawRate;
			mp->yawErrorI += (yawError)*DeltaT;
			

			// double yaw_Moment = mp->K_P*(yawError)+mp->K_I*mp->yawErrorI;
			double yaw_Moment = mp->K_P*yawError;
			double correctionTorqueF = yaw_Moment*mp->wheelRadius/2.0/mp->halfWidth;
			/* Gas */
			double targetTorque = IF->Gas*IF->MotorIn[0].TrqMot_max;
			IF->MotorOut[0].Trq_trg = IF->Gas*IF->MotorIn[0].TrqMot_max-correctionTorqueF;
			IF->MotorOut[1].Trq_trg = IF->Gas*IF->MotorIn[1].TrqMot_max+correctionTorqueF;
			
			IF->MotorOut[2].Trq_trg = targetTorque*mp->FL_Ratio/mp->R_Ratio*2;
			if(SimCore.TimeWC-mp->prevTime > 1) {
				Log("YawError: %f\n",yawError);
				Log("YawErrorI: %f\n",mp->yawErrorI);
				// Log("Understeer: %f\n",mp->K_U);
				mp->prevTime = SimCore.TimeWC;
			}
		}
		else {
			mp->yawErrorI = 0;
			double targetTorque = IF->Gas*IF->MotorIn[0].TrqMot_max;
			IF->MotorOut[0].Trq_trg = targetTorque;
			IF->MotorOut[1].Trq_trg = targetTorque;
			IF->MotorOut[2].Trq_trg = targetTorque*mp->FL_Ratio/mp->R_Ratio*2;	
		}
	   
	    break;
    }
	
    OutOfOperState:
		// IF->MotorOut[0].Trq_trg = 0;
		// IF->MotorOut[1].Trq_trg = 0;
		// IF->MotorOut[2].Trq_trg = 0;

    return 0;
}


int
PTControl_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTControl.VersionId =		ThisVersionId;
    m.PTControl.New =			MyModel_New;
    m.PTControl.Calc =			MyModel_Calc;
    m.PTControl.Delete =		MyModel_Delete;
    m.PTControl.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTControl.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTControl, ThisModelKind, &m);
}

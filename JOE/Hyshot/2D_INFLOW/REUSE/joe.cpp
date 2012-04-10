#include "JoeWithModels.h"

#include "turbModels/TurbModel_KOM.h"
#include "turbModels/TurbModel_KEps.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "turbModels/TurbModel_SA.h"
#include "turbModels/TurbModel_V2F.h"

//#include "combModels/CombModel_Mixing.h"
//#include "combModels/CombModel_SteadyFlamelet.h"
#include "combModels/CombModel_VariableProperties.h"


// ###########################################################################################
// ------                                                                               ------
// ------                    Vanilla version of joe called MyJoe                        ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe : public JoeWithModels 
{
public:

  MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) 
  {
    if (mpi_rank == 0)      cout << "MyJoe() with Combustion" << endl;
    deta_dx = NULL;      registerScalar(deta_dx, "deta_dx",  CV_DATA);

    if (checkParam("MACH"))
    {
      Remove("inlet");

      double Ma  = getDoubleParam("MACH");
      double aoa = M_PI/180.0*getDoubleParam("AOA");

      double sos = sqrt(GAMMA*p_ref/rho_ref);
      double speed = Ma*sos;

      double velx = speed*cos(aoa);
      double vely = speed*sin(aoa);

      char name[200];
      sprintf(name, "inlet = CBC %.6lf %.6lf 0.0 %.6lf %.6lf", velx, vely, T_ref, p_ref);
      ParamMap::add(string(name));

      if (mpi_rank == 0) {
	cout << "Mach input detected, rewriting 'inlet' BC to:" << endl;
	cout << string(name) << endl;
      }
    }
  }

  virtual ~MyJoe()  {}

  double *deta_dx;
  
  double xKnick, comb_d, comb_complete, comb_k, comb_xc, comb_yc;
  double *src;
  double C1, C2, Lmix;
  double xCombIn, xCombOut, specificQ;
  double tStart, tEnd;

  void initialHook()
  {

    if (mpi_rank == 0)
        cout << "MyJoe InitialHook()" << endl;

    JoeWithModels::initialHook();

    src = new double[ncv];
    
    Param *p = getParam("EQUIVALENCE_RATIO");
    tStart = p->getDouble(2);
    tEnd = p->getDouble(3);

    double phi    = p->getDouble(1);
    double Hcomb  = getDoubleParam("H_COMB", "1.2e8");
    double f      = getDoubleParam("F_STOICHIOMETRIC", "0.029");
    comb_d        = getDoubleParam("COMB_D", "1.0");
    comb_xc        = getDoubleParam("COMB_XC", "0.0");
    comb_yc        = getDoubleParam("COMB_YC", "0.0");
    comb_complete = getDoubleParam("COMB_COMPLETION_RATE", "0.95");
    
    comb_k = pow(-log(1.0-comb_complete), comb_d);
    
    if (mpi_rank == 0)
      cout << "computed k = " << comb_k << endl;

    xCombIn  = 0.40615 + comb_xc;
    xKnick = 0.649873;
    xCombOut = 0.773120;
    Lmix  = xCombOut-xCombIn;

    double mDotOxidizer = computeMassFlowRate(0);

//  specificQ = phi*f*mDotOxidizer*Hcomb;
    specificQ = phi*f*Hcomb;
    
    if (mpi_rank == 0)
        cout << "specificQ k = " << specificQ << " mDotOxidizer = " << mDotOxidizer << endl;
//    cout <<specificQ<<" "<<phi<<" "<<f<<" "<<mDotOxidizer<<" "<<Hcomb<<" "<<Lmix<<" "<<Area<<endl;
    
    
    for (int icv=0; icv<ncv; icv++)
    {
      double xCoord = x_cv[icv][0];

      double g=-5.0+4.0*comb_yc+(x_cv[icv][1]-0.11952)/(0.12932-0.11952)*10.0;
      if (xCoord > xCombIn)
      {
        double xdl = (xCoord-xCombIn)/Lmix;
        deta_dx[icv] = comb_d*pow(comb_k*xdl, comb_d-1.0)*comb_k/Lmix*exp(-pow(comb_k*xdl, comb_d));
        deta_dx[icv] *= 1.0-tanh(g);
      }
      else                        
        deta_dx[icv] = 0.0;
    }
  }

  virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5])
  {
    if (getIntParam("COMBUSTION") == 1)
    {
      double mDotOxidizer = computeMassFlowRate(0);
      for (int icv=0; icv<ncv; icv++)
      {
        double xCoord = x_cv[icv][0];
        double area = getArea(xCoord);
        rhs_rhoE[icv] += mDotOxidizer*specificQ*deta_dx[icv]*cv_volume[icv]/area;
      }
    }
  }

  double computeMassFlowRate(int output)
  {
    double x0 = 0.4020;
    double x1 = 0.4022;
    double dx = x1 - x0;
    double mfCPU = 0.0;
    double area = 0.0;
      
    for (int icv = 0; icv < ncv; icv++)
    {
      if ((x_cv[icv][0] > x0) && (x_cv[icv][0] < x1))
      {
        int foc_f = faocv_i[icv];
        int foc_l = faocv_i[icv + 1] - 1;

        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          int ifa = faocv_v[foc];
          if (x_fa[ifa][0] > x_cv[icv][0] + 0.45 * dx)
            area = fabs(fa_normal[ifa][0]);
        }

        mfCPU += rhou[icv][0] * area;
      }
    }
    
    double massflow = 0.0;
    MPI_Allreduce(&mfCPU, &massflow, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    return massflow;
  }

  double getArea(double x)
  {
    double area = (0.129320007-0.119519997)*(0.009375-0.0);

    if (x>xKnick)
      area += (x-xKnick)/(xCombOut-xKnick)*(0.155597-0.129320007)*(0.009375-0.0);
    
    return area;
  }


  /*
   * write wall values to file in tecplot format
   */
  virtual void writeWallValues(string name)
  {
    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == name)
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            if (x_fa[ifa][1] > 0.1194)
              my_wall_faces++;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "%s.dat", name.c_str());
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"rho\" \"press\" \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"muLam\"\n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", tot_wall_faces);
    }
    else
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      if ( (fp=fopen(fname,"a"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
    }

    for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone != faZoneList.end(); faZone++)
      if (faZone->getKind() == FA_ZONE_BOUNDARY)
        if (faZone->getNameString() == name)
        {
          for (int ifa = faZone->ifa_f; ifa <= faZone->ifa_l; ifa++)
          if (x_fa[ifa][1] > 0.1194)
          {
            int icv = cvofa[ifa][0];

            double n[3], s_half[3], vel[3], velTang[3];
            
            normVec3d(n, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv]);

            vel[0] = rhou[icv][0]/rho[icv];
            vel[1] = rhou[icv][1]/rho[icv];
            vel[2] = rhou[icv][2]/rho[icv];
            double un = vecDotVec3d(n, vel);

            velTang[0] = vel[0] - un*n[0];
            velTang[1] = vel[1] - un*n[1];
            velTang[2] = vel[2] - un*n[2];

            double velMag = sqrt(vecDotVec3d(velTang, velTang));
            double walld = fabs(vecDotVec3d(s_half, n));

            double wallTemp = 300.0;
            double mulam = mul_fa[ifa];
            
            double tau = mulam*velMag/walld*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = walld*utau/(mulam/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mulam/Pr;
            double qDot = kTotal*(temp[icv]-wallTemp)/walld;

            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
                x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv], temp[icv], tau, qDot, yplus, mul_fa[ifa]);
          }
        }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }

  virtual void temporalHook()
  {
    if (step%10 == 0)
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
        if (zone->getKind() == FA_ZONE_BOUNDARY)
        {
          Param *param;
          if (getParam(param, zone->getName()))
            if (param->getString() == "WALL")
              writeWallValues(zone->getName());
        }
  }

  virtual void finalHook()
  {
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            writeWallValues(zone->getName());
      }
  }
};


/*
 * MyJoe with Spalart & Allmaras Model
 */
class MyJoeSA : public MyJoe, public RansTurbSA
{
public:
  MyJoeSA(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSA()" << endl;
  }

  virtual ~MyJoeSA() {}
};

/*
 * MyJoe with Menter SST Model
 */
class MyJoeSST : public MyJoe, public RansTurbKOmSST
{
public:
   double *SST_F1;            ///< blending function parameter
   double *trans_scale;       ///< transition scaling parameter
   double x_body, x_cowl, x_ramp;

public:

  MyJoeSST(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0) cout << "MyJoeSST()" << endl;

    if ( checkParam("TRANSITION_COWL") )
    {
      x_cowl = getParam("TRANSITION_COWL")->getDouble("X");
      if (mpi_rank == 0) cout << "COWL TRANSITION = " << x_cowl << endl;
    }
    if ( checkParam("TRANSITION_BODY") )
    {
      x_body = getParam("TRANSITION_BODY")->getDouble("X");
      if (mpi_rank == 0) cout << "BODY TRANSITION = " << x_body << endl;
    }
    if ( checkParam("TRANSITION_RAMP") )
    {
      x_ramp = getParam("TRANSITION_RAMP")->getDouble("X");
      if (mpi_rank == 0) cout << "RAMP TRANSITION = " << x_ramp << endl;
    }

    SST_F1 = NULL;         registerScalar(SST_F1,   "SST_F1",   CV_DATA);
    trans_scale = NULL;    registerScalar(trans_scale,   "trans_scale",   CV_DATA);
  }

  virtual ~MyJoeSST() {}

  void initialHook()
  {
/*    double dum1, dum2;

    for (int icv = 0; icv < ncv; icv++)
    {
      rho[icv] = press[icv] / (RoM[icv]*temp[icv]);
      rhou[icv][0] = rho[icv]*vel[icv][0];
      rhou[icv][1] = rho[icv]*vel[icv][1];
      rhou[icv][2] = rho[icv]*vel[icv][2];

      double H_inlet, dum1, dum2;
      calcThermoProp_T(press[icv], H_inlet, RoM[icv], gamma[icv], rho[icv], temp[icv], NULL, 0);

      double Ekin = 0.5 * (rhou[icv][0]*rhou[icv][0] + rhou[icv][1]*rhou[icv][1] + rhou[icv][2]*rhou[icv][2]) / rho[icv];
      rhoE[icv] = H_inlet * rho[icv] - press[icv] + Ekin;
    }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
*/

    // Set BL transition scaling flag only once - in initialHook
	for (int icv=0; icv<ncv; icv++)
	{
      trans_scale[icv] = 1.0;
      double wall_offset = 0.0;

      // Cowl chamber transition
      if ( checkParam("TRANSITION_COWL") )
      {
        // enforce transition on cowl (upper) combustor wall
        if ( (x_cv[icv][1] > 0.1173  && x_cv[icv][0] > 0.289 ) &&
             (x_cv[icv][1] > 0.12436 && x_cv[icv][0] < x_cowl) )
        {
          wall_offset = 0.00007;
          trans_scale[icv] = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
        }
      }

      if ( checkParam("TRANSITION_BODY") )
      {
        // enforce transition on body (lower) combustor wall
        if ( (x_cv[icv][1] > 0.1173  && x_cv[icv][0] > 0.289 ) &&
             (x_cv[icv][1] < 0.12436 && x_cv[icv][0] < x_body) )
        {
          wall_offset = 0.00008;
          trans_scale[icv] = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
        }
      }

      if ( checkParam("TRANSITION_RAMP") )
      {
        // enforce transition on intake ramp wall
        if ( (x_cv[icv][1] < 0.1173  && x_cv[icv][0] < x_ramp ) )
        {
          wall_offset = 0.0003;
          trans_scale[icv] = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
        }
      }
	}

	updateCvData(trans_scale, REPLACE_DATA);
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
    {
      double Pk = 0.0;
//      double scale = 1.0;       // scale factor for transition specification
//      double wall_offset = 0.0; //wall distance at which boundary layer is *defined*

      for (int icv=0; icv<ncv; icv++)
      {
        // compute turbulence production
	    double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
        double mut = min(max(rho[icv]*kine[icv]*zeta, 0.0), 1.0e5);

        Pk = mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];

        // production limiter of ...
        if ( SST_limitPK == 1 )  Pk = min(Pk, 20.0*betaStar*rho[icv]*kine[icv]*omega[icv]);

        // enforce production >= 0
        Pk = max(Pk, 0.0);

//        // Cowl chamber transition
//        if ( checkParam("TRANSITION_COWL") )
//        {
//	      // enforce transition on cowl (upper) combustor wall
//          if ( (x_cv[icv][1] > 0.1173  && x_cv[icv][0] > 0.289 ) &&
//               (x_cv[icv][1] > 0.12436 && x_cv[icv][0] < x_cowl) )
//          {
//            wall_offset = 0.00007;
//            scale = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
//          }
//        }
//
//        if ( checkParam("TRANSITION_BODY") )
//        {
//	      // enforce transition on body (lower) combustor wall
//          if ( (x_cv[icv][1] > 0.1173  && x_cv[icv][0] > 0.289 ) &&
//               (x_cv[icv][1] < 0.12436 && x_cv[icv][0] < x_body) )
//          {
//            wall_offset = 0.00008;
//            scale = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
//          }
//        }
//
//        if ( checkParam("TRANSITION_RAMP") )
//        {
//	      // enforce transition on intake ramp wall
//          if ( (x_cv[icv][1] < 0.1173  && x_cv[icv][0] < x_ramp ) )
//          {
//            wall_offset = 0.0003;
//            scale = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
//          }
//        }
//
//        trans_scale[icv] = scale;

        /* enforce transition according to Menter gamma-Re_theta criterion:
           - Pk in k-equation goes to zero in laminar BL
           - damping effect on dissipation term in k-equation [0.1-1]
           - F1 blending term (omega equation) goes to 1.0 in laminar BL
         */
        double src = trans_scale[icv]*Pk - (0.1+0.9*trans_scale[icv])*betaStar*rho[icv]*omega[icv]*kine[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - (0.1+0.9*trans_scale[icv])*betaStar*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }

    if (name == "omega")
    {
      double Pk = 0.0;
//      double scale = 1.0;       // scale factor for transition specification
//      double wall_offset = 0.0; //wall distance at which boundary layer is *defined*

      for (int icv=0; icv<ncv; icv++)
      {
        double F1 = blendFuncF1[icv];
        SST_F1[icv] = F1;

        double zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
        Pk = strMag[icv]*strMag[icv] - 2./3.*zeta*diverg[icv];
        Pk = max(Pk, 0.0);

//        // Cowl chamber transition
//        if ( checkParam("TRANSITION_COWL") )
//        {
//	      // enforce transition on cowl (upper) combustor wall
//          if ( (x_cv[icv][1] > 0.1173  && x_cv[icv][0] > 0.289 ) &&
//               (x_cv[icv][1] > 0.12436 && x_cv[icv][0] < x_cowl) )
//          {
//            wall_offset = 0.00007;
//            scale = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
//          }
//        }
//
//        if ( checkParam("TRANSITION_BODY") )
//        {
//	      // enforce transition on body (lower) combustor wall
//          if ( (x_cv[icv][1] > 0.1173  && x_cv[icv][0] > 0.289 ) &&
//               (x_cv[icv][1] < 0.12436 && x_cv[icv][0] < x_body) )
//          {
//            wall_offset = 0.00008;
//            scale = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
//          }
//        }
//
//        if ( checkParam("TRANSITION_RAMP") )
//        {
//	      // enforce transition on intake ramp wall
//          if ( (x_cv[icv][1] < 0.1173  && x_cv[icv][0] < x_ramp ) )
//          {
//            wall_offset = 0.0003;
//            scale = 0.5*(1.0 + tanh( (wallDist[icv] - wall_offset)/(0.1*wall_offset) ));
//          }
//        }

        double F1_diff = 1.0-F1;

        F1 = 1.0 - trans_scale[icv]*F1_diff;

	    double alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
	    double beta = F1*beta_1 + (1.0 - F1)*beta_2;

        double src = alfa*rho[icv]*Pk + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
  }

};

/*
 * MyJoe with WILCOX Model
 */
class MyJoeWCX : public MyJoe, public RansTurbKOm
{
public:
  MyJoeWCX(char *name) : MyJoe(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeWCX()" << endl;
  }

  virtual ~MyJoeWCX() {}

};


/*
 * MyJoe with SA Model and variable properties
 */
class MyJoeSA_VP : public MyJoeSA, public RansCombVarProperties
{
public:
  MyJoeSA_VP(char *name) : MyJoeSA(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSA_VP()" << endl;
  }

  virtual ~MyJoeSA_VP() {}

  void initialHook()
  {
    double dum1, dum2;

    for (int icv = 0; icv < ncv; icv++)
    {
      rho[icv] = press[icv] / (RoM[icv]*temp[icv]);
      rhou[icv][0] = rho[icv]*vel[icv][0];
      rhou[icv][1] = rho[icv]*vel[icv][1];
      rhou[icv][2] = rho[icv]*vel[icv][2];

      double H_inlet, dum1, dum2;
      calcThermoProp_T(press[icv], H_inlet, RoM[icv], gamma[icv], rho[icv], temp[icv], NULL, 0);

      double Ekin = 0.5 * (rhou[icv][0]*rhou[icv][0] + rhou[icv][1]*rhou[icv][1] + rhou[icv][2]*rhou[icv][2]) / rho[icv];
      rhoE[icv] = H_inlet * rho[icv] - press[icv] + Ekin;
    }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
  }
};


/*
 * MyJoe with Menter SST Model and variable properties
 */
class MyJoeSST_VP : public MyJoeSST, public RansCombVarProperties
{
public:
  MyJoeSST_VP(char *name) : MyJoeSST(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeSST_VP()" << endl;
  }
  
  virtual ~MyJoeSST_VP() {}
  
  void initialHook()
  {

    if (mpi_rank == 0)
      cout << "set JoeSST_VP InitialCondition()"<< endl;

    MyJoe::initialHook();
    MyJoeSST::initialHook();


    // only set rho, rhou[i],rhoE if starting from scratch, not if from restart
    if (!checkScalarFlag("RHO"))
      for (int icv = 0; icv < ncv; icv++)
      {
        rho[icv] = press[icv] / (RoM[icv]*temp[icv]);
        rhou[icv][0] = rho[icv]*vel[icv][0];
        rhou[icv][1] = rho[icv]*vel[icv][1];
        rhou[icv][2] = rho[icv]*vel[icv][2];

        double H_inlet;
        calcThermoProp_T(press[icv], H_inlet, RoM[icv], gamma[icv], rho[icv], temp[icv], NULL, 0);

        double Ekin = 0.5 * (rhou[icv][0]*rhou[icv][0] + rhou[icv][1]*rhou[icv][1] + rhou[icv][2]*rhou[icv][2]) / rho[icv];
        rhoE[icv] = H_inlet * rho[icv] - press[icv] + Ekin;
      }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
  }
};



/*
 * MyJoe with Wilcox k-omega Model and variable properties
 */
class MyJoeWCX_VP : public MyJoeWCX, public RansCombVarProperties
{
public:
  MyJoeWCX_VP(char *name) : MyJoeWCX(name), UgpWithCvCompFlow(name)
  {
    if (mpi_rank == 0)
      cout << "MyJoeWCX_VP()" << endl;
  }

  virtual ~MyJoeWCX_VP() {}

  void initialHook()
  {
    double dum1, dum2;

    for (int icv = 0; icv < ncv; icv++)
    {
      rho[icv] = press[icv] / (RoM[icv]*temp[icv]);
      rhou[icv][0] = rho[icv]*vel[icv][0];
      rhou[icv][1] = rho[icv]*vel[icv][1];
      rhou[icv][2] = rho[icv]*vel[icv][2];

      double H_inlet, dum1, dum2;
      calcThermoProp_T(press[icv], H_inlet, RoM[icv], gamma[icv], rho[icv], temp[icv], NULL, 0);

      double Ekin = 0.5 * (rhou[icv][0]*rhou[icv][0] + rhou[icv][1]*rhou[icv][1] + rhou[icv][2]*rhou[icv][2]) / rho[icv];
      rhoE[icv] = H_inlet * rho[icv] - press[icv] + Ekin;
    }

    // update all data
    updateCvData(rhou, REPLACE_ROTATE_DATA);
    updateCvData(rho, REPLACE_DATA);
    updateCvData(rhoE, REPLACE_DATA);
  }
};






int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  // the run number specifies the class which is going to be instantiated
  // set run to default value of 0 => instantiation of UgpWithCvCompFlow
  int run = 0;

  // set input name to default "Joe.in"
  char inputFileName[50];
  sprintf(inputFileName, "Joe.in");

  for (int i=1; i<argc; i++)
  {
    string str(argv[i]);
    if (from_string<int>(run, str, std::dec))
    {
      if (mpi_rank == 0)
        cout << "You have specified run number = " << run << endl;
    }
    else
      strcpy(inputFileName, argv[i]);
  }

  if (mpi_rank == 0)
  {
    cout << "SPECIFIED INPUT NAME = " << inputFileName << endl;
    cout << "SPECIFIED RUN = " << run << endl;
  }


  try {

    // declare pointer to JoeWithModels
    JoeWithModels *joe;
    
    switch (run)
    {
    case 0:   joe = new MyJoe(inputFileName);         break;
    case 1:   joe = new MyJoeSA(inputFileName);       break;
    case 2:   joe = new MyJoeSST(inputFileName);      break;
    case 3:   joe = new MyJoeWCX(inputFileName);      break;
    case 10:  joe = new MyJoeSA_VP(inputFileName);    break;
    case 11:  joe = new MyJoeSST_VP(inputFileName);   break;
    case 12:  joe = new MyJoeWCX_VP(inputFileName);   break;
    default: 
      if (mpi_rank == 0)
        cerr << "ERROR: run number not available!" << endl;
      throw(-1);
    }
    
    // provide total runtime 
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    // run joe
    joe->run();
    

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

    // delete joe (make sure memory is deallocated in destructors
    delete joe;
  }
  catch (int e) {
    cerr << "Exception: " << e << endl;
    MPI_Finalize();
    return(-1);
  }
  catch(...) {
    cerr << "unhandled Exception.\n" << endl;
    MPI_Finalize();
    return(-1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return (0);
}




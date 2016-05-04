
#include "ComputeIrrotational.hpp"
#include "Utils/DataMesh/Mesh.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "PointwiseFunctions/TensorFunctions/SpacePlusTimeDecomposition.hpp"
#include "Utils/StringParsing/StringUtils.hpp"
#include "Dust/Functionals/GlobalDifferentiator.hpp"
#include "Hydro/HydroStates/AttenuatedBinary.hpp"
#include "Dust/Domain/Subdomain/Subdomain.hpp"
#include "Utils/LowLevelUtils/ConvertNumberToString.hpp"

namespace ComputeItems {

  using std::string;


  
  ///////////////////////////////////////////////////////////////////////
  ////////////////////ROTATION TERM/////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ComputeRotationTerm::ComputeRotationTerm(const string& opts):
     mpResult(0) {
     OptionParser p(opts,Help());
     mOutput=p.Get<string>("Output");
     mPsi=p.Get<string>("ConformalFactor");
     mg=p.Get<string>("ConformalMetric");
     mCenterNS=p.Get<string>("CenterNS");
     mStarRotation=p.Get<string>("StarRotation");
     mSubd=p.Get<MyVector<std::string> >("Subdomains");
     mCFPower=p.Get<double>("CFPower");
   }

   void ComputeRotationTerm::RecomputeData(const DataBoxAccess& box)
   const {
     
     std::string mname=box.Get<Subdomain>("Subdomain").Name();
     bool compute=false;
     for (int i =0; i <mSubd.Size(); i++){
       if (mSubd[i]==mname){
 	compute=true;
       }
     }
    
     if (compute){
       const Dm& Psi = box.Get<TDm>(mPsi)();
       const TDm& g = box.Get<TDm>(mg);
       const MyVector<DataMesh>& coords=box.Get<MyVector<DataMesh> >("GlobalCoords");
       const MyVector<double>& CenterNS = box.Root().Get<MyVector<double> >(mCenterNS);
       const TDm& StarRotation = box.Get<TDm>(mStarRotation);
 
       const Mesh& sd=box.Get<Mesh>("Mesh");
       if(mpResult==0)
 	mpResult = new Tensor<DataMesh>(3, "1", sd);
       TDm temp(3,"1",sd,0.0);   
	  	  
       temp(0) = StarRotation(1)*(coords[2]-CenterNS[2]) - StarRotation(2)*(coords[1]-CenterNS[1]);
       temp(1) = StarRotation(2)*(coords[0]-CenterNS[0]) - StarRotation(0)*(coords[2]-CenterNS[2]);
       temp(2) = StarRotation(0)*(coords[1]-CenterNS[1]) - StarRotation(1)*(coords[0]-CenterNS[0]);
 
       for (int i=0; i<3; i++){
 	(*mpResult)(i)=0;
       }
     
       for (int i =0; i<3; i++){
 	for (int j=0; j<3; j++){
 	  (*mpResult)(i) += g(i,j)*temp(j)/sqr(Psi);
 	}
       }

       
       for (int i=0; i<3; i++){
	 (*mpResult)(i) *= pow(Psi, mCFPower);
       }
	 
     } else {
       const Mesh& sd=box.Get<Mesh>("Mesh");
       if(mpResult==0)
 	mpResult = new Tensor<DataMesh>(3, "1", sd);
       for (int i=0; i <3; i++){
 	(*mpResult)(i)=0.0;
       }
       
     }	    
   }




  ///////////////////////////////////////////////////////////////////////
  /////////////////////////INERTIAL VELOCITY/////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  ComputeInertialVelocity::ComputeInertialVelocity(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mShift=p.Get<string>("Shift");
    mdPot=p.Get<string>("dHydroPot","dHydroPot");
    mInvg=p.Get<string>("InvConformalMetric","InvConformalMetric");
    mOmega=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mC=p.Get<string>("EulerConstant","EulerConstant");
    mCenterW=p.Get<string>("CenterW","CenterW");
    mCenterNSName=p.Get<string>("CenterNS","CenterNS");
    mdtExpansionFactorName = p.Get<string>("dtExpansionFactor",
					   "dtExpansionFactor");
    mSolidBodyRotation=p.Get<MyVector< double> >("SolidBodyRotationFrequency",MyVector<double>(MV::Size(3),0));
    mRot = p.Get<string>("RotationTerm", "RotationTerm");
  }
  
  void ComputeInertialVelocity::RecomputeData(const DataBoxAccess& box) 
    const {
    Dm Psi =  box.Get<TDm>(mPsi)();
    Dm NPsi =  box.Get<TDm>(mNPsi)();
    TDm InShift =  box.Get<TDm>(mShift);
    TDm Invg =  box.Get<TDm>(mInvg);
    TDm dPot =  box.Get<TDm>(mdPot);
    const double CEuler=box.Root().Get<double>(mC);
    MyVector<double> W=box.Root().Get<MyVector<double> >(mCenterW);
    MyVector<double> mCenterNS=box.Root().Get<MyVector<double> >(mCenterNSName);
    TDm Rot = box.Get<TDm>(mRot);

    Dm Lapse=NPsi/Psi;

    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "1", sd);
    
    //get OmegaOrbit from the databox
    const MyVector<double> Omega = box.Root().Get<MyVector<double > >(mOmega);
    const MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");


    double mdtExpansionFactor = box.Root().Get<double >(mdtExpansionFactorName);
    TDm mdtExpansionFactorTimesrHi(3,"1",sd,0.);
    mdtExpansionFactorTimesrHi(0) = mdtExpansionFactor * coords[0];
    mdtExpansionFactorTimesrHi(1) = mdtExpansionFactor * coords[1];
    mdtExpansionFactorTimesrHi(2) = mdtExpansionFactor * coords[2];

    //Omegaxr
    TDm OmegaxrHi = InShift;
    OmegaxrHi(0) = Omega[1]*coords[2] - Omega[2]*coords[1];
    OmegaxrHi(1) = Omega[2]*coords[0] - Omega[0]*coords[2];
    OmegaxrHi(2) = Omega[0]*coords[1] - Omega[1]*coords[0];

    TDm SolidRotationV=InShift;
    SolidRotationV(0) = mSolidBodyRotation[1]*(coords[2]-mCenterNS[2])
      -mSolidBodyRotation[2]*(coords[1]-mCenterNS[1]);
    SolidRotationV(1) = mSolidBodyRotation[2]*(coords[0]-mCenterNS[0])
      -mSolidBodyRotation[0]*(coords[2]-mCenterNS[2]);
    SolidRotationV(2) = mSolidBodyRotation[0]*(coords[1]-mCenterNS[1])
      -mSolidBodyRotation[1]*(coords[0]-mCenterNS[0]);

    TDm Shift=InShift;
    for(int i=0;i<3;i++) Shift(i)+=OmegaxrHi(i)
			       +mdtExpansionFactorTimesrHi(i);

    Dm tempA=Lapse*0.0+CEuler;
    TDm tempB=Shift;
    for(int i=0;i<3;i++)
      {
	tempB(i)=0.0;
	tempA+=(dPot(i)+W[i]+Rot(i))*Shift(i);
	for(int j=0;j<3;j++){
	  tempB(i)+=Invg(i,j)*(dPot(j)+W[j]+Rot(j));
	}
	tempB(i)/=sqr(sqr(Psi));
      }
    for(int i=0;i<3;i++)
      {
	(*mpResult)(i)=tempB(i)/tempA*Lapse+SolidRotationV(i)/Lapse;
      }

  }

  
  ///////////////////////////////////////////////////////////////////////
  /////////////////////////Enthalpy//////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  
  ComputeIrrotationalEnthalpy::ComputeIrrotationalEnthalpy(const string& opts):
    mSubd(MV::Size(1),""),mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mShift=p.Get<string>("Shift");
    mU=p.Get<string>("InertialVelocity","InertialVelocity");
    mg=p.Get<string>("ConformalMetric","ConformalMetric");
    mOmega=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mC=p.Get<string>("EulerConstant","EulerConstant");
    mDrag=p.Get<string>("DragParam","DragParam");
    mSubd=p.Get<MyVector<std::string> >("Subdomains");
    mCenter=p.Get<std::string>("Center");
    mdtExpansionFactorName = p.Get<string>("dtExpansionFactor",
					   "dtExpansionFactor");
    mdPot = p.Get<string>("dHydroPot", "dHydroPot");
    mInvg = p.Get<string>("InvConformalMetric", "InvConformalMetric");
    mRot = p.Get<string>("RotationTerm", "RotationTerm");
    mCenterW = p.Get<string>("CenterW", "CenterW");
    }
  
  void ComputeIrrotationalEnthalpy::RecomputeData(const DataBoxAccess& box) 
    const {
       std::string mname=box.Get<Subdomain>("Subdomain").Name();
    bool compute=false;
    for (int i =0; i <mSubd.Size(); i++){
      if (mSubd[i]==mname){
	compute=true;
      }
    }
    const Mesh& sd=box.Get<Mesh>("Mesh");
    Dm rad(sd);
    if (compute)
      {

	const Dm& Psi = box.Get<TDm>(mPsi)();
	const Dm& NPsi = box.Get<TDm>(mNPsi)();
	const Dm Lapse = NPsi / Psi;
	const TDm& dPot = box.Get<TDm>(mdPot);
	const TDm& InShift = box.Get<TDm>(mShift);
	const TDm& Invg = box.Get<TDm>(mInvg);
	const MyVector<double> Omega=box.Root().Get<MyVector<double> >(mOmega);
	const MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");
	const double CEuler=box.Root().Get<double>(mC);
	const MyVector<double> W = box.Root().Get<MyVector<double> >(mCenterW);
	const MyVector<double> DragParam = box.Root().Get<MyVector<double> >(mDrag);
	const MyVector<double> Center = box.Root().Get<MyVector<double> >(mCenter);	
	const TDm& Rot = box.Get<TDm>(mRot);
	
	if (mpResult==0){
	  mpResult = new Tensor<DataMesh>(3, "", sd);
	}
	
	rad=0.0;
	for (int i =0; i<3; i++){
	  rad += sqr(coords[i]-Center[i]);
	}
	rad=sqrt(rad);

	double mdtExpansionFactor = box.Root().Get<double >(mdtExpansionFactorName);
	TDm mdtExpansionFactorTimesrHi(3, "1", sd, 0.);
	mdtExpansionFactorTimesrHi(0) = mdtExpansionFactor * coords[0];
	mdtExpansionFactorTimesrHi(1) = mdtExpansionFactor * coords[1];
	mdtExpansionFactorTimesrHi(2) = mdtExpansionFactor * coords[2];

	TDm OmegaxrHi = InShift;
	OmegaxrHi(0) = Omega[1]*coords[2] - Omega[2]*coords[1];
	OmegaxrHi(1) = Omega[2]*coords[0] - Omega[0]*coords[2];
	OmegaxrHi(2) = Omega[0]*coords[1] - Omega[1]*coords[0];

	TDm Shift = InShift;
	for (int i=0; i<3; i++){
	  Shift(i) += OmegaxrHi(i) +mdtExpansionFactorTimesrHi(i);
	}


	TDm moddPot = dPot;
	for (int i =0; i <3; i++){
	  moddPot(i) += W[i];
	}

	
	Dm b1(sd,0.0);
	for (int i=0; i<3; i++){
	  b1 += Shift(i)*moddPot(i);
	}
	b1 = sqr(b1 + CEuler);

	Dm b2(sd,0.0);
	for (int i=0; i<3; i++){
	  for (int j=0; j<3; j++){
	    b2 += Invg(i,j) / (sqr(sqr(Psi))) * (moddPot(i)+Rot(i)) * Rot(j);
	  }
	}
	b2 *= (2 * Lapse * Lapse);

	Dm L2(sd, 0.0);
	L2 = ((b1+b2) + sqrt(sqr(b1)+2*b1*b2))/(2*sqr(Lapse));
	
	Dm H = L2;
	for (int i=0; i <3; i++){
	  for (int j=0; j<3; j++){
	    H -= Invg(i,j)/(sqr(sqr(Psi))) * (moddPot(i)+Rot(i))*(moddPot(j)+Rot(j));
	  }
	}
	H = sqrt(H);
			    
	
			      

	/*

	 

	Dm b1(sd,0.0);
	for (int i=0;i<3;i++){
	  b1 += Shift(i)*moddPot(i);
	}
	b1 = sqr(b1 + CEuler);

	Dm b2(sd,0.0);
	for (int i=0;i<3;i++){
	  for (int j=0;j<3;j++){
	    b2 += Invg(i,j)/(sqr(sqr(Psi)))*(moddPot(i)+Rot(i))*Rot(j);
	  }
	}

	Dm H = sqr(b1) + 2*b1*b2;
	H = sqrt(H);
	H += (b1+b2);
	H = H / (2*sqr(Lapse));
	
	for (int i =0; i <3; i++){
	  for (int j=0; j<3; j++){
	    H -= Invg(i,j) / (sqr(sqr(Psi))) * (moddPot(i)) * (moddPot(j));
	  } 
	}
	H = sqrt(H);

	*/

	for (int dcur=0; dcur<3; dcur++){
	  H *= (1.0 - DragParam[dcur]*(coords[dcur]-Center[dcur]));
	}

	for (int i =0; i< H.Size(); i++){
	  if (H[i] >0.1){
	    (*mpResult)()[i] = H[i];
	  } else {
	    (*mpResult)()[i] = 0.1;
	  }
	  
	}
      } else {
      if (mpResult==0){
	mpResult = new Tensor<DataMesh>(3,"",sd);
      }
      Dm H(sd);
      H=0.1;
      (*mpResult)()=H;
    }
  }
  

  ///////////////////////////////////////////////////////////////////////
  /////////////////////////Conformal metric//////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  BinaryConformalMetric::BinaryConformalMetric(const string& opts):
    mSpin(MV::fill,0,0,0),mpResult(0) {
    OptionParser p(opts, Help());
    mMassBH=p.Get<double>("MassBH",1.0);
    mCenterBH=p.Get<string>("CenterBH","CenterBH");
    mCenterNS=p.Get<string>("CenterNS","CenterNS");
    mOmegaOrbit=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mSpin=p.Get<MyVector<double> >("Spin");
    mdta=p.Get<string>("dtExpansionFactor","dtExpansionFactor");
    rhoc=p.Get<double>("rhoc");
    mEOS=p.Get<string>("EOS");
    AttenuationRadiusBH=p.Get<double>("AttenuationRadiusBH");
    HorRadius=p.Get<double>("HorizonBH",0.0);
    Ratio=p.Get<double>("Ratio",1.0);
    doBoost=p.Get<bool>("Boost",true);
    AttenuationRadiusNS=p.Get<double>("AttenuationRadiusNS");
    mOutput=p.Get<string>("Output");
    }
  
  void BinaryConformalMetric::RecomputeData(const DataBoxAccess& box) 
    const {

    MyVector<double> CenterBH=box.Root().Get<MyVector<double> >(mCenterBH);
    MyVector<double> CenterNS=box.Root().Get<MyVector<double> >(mCenterNS);
    double OmegaOrbit=box.Root().Get<MyVector<double> >(mOmegaOrbit)[2];
    double dta=box.Root().Get<double>(mdta);

    string binopts="BH=KerrSchild(Mass=";
    binopts+=DoubleToStringDecimal(mMassBH,16)+string(";Center=");
    binopts+=DoubleToStringDecimal(CenterBH[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[2],16)+string("; Spin=");		  
    binopts+=DoubleToStringDecimal(mSpin[0],16)+string(",");		  
    binopts+=DoubleToStringDecimal(mSpin[1],16)+string(",");
    binopts+=DoubleToStringDecimal(mSpin[2],16)+string("; Velocity=");
    if(doBoost)
    {
        binopts+=DoubleToStringDecimal(dta*CenterBH[0]-OmegaOrbit*CenterBH[1],10)+string(",");
        binopts+=DoubleToStringDecimal(dta*CenterBH[1]+OmegaOrbit*CenterBH[0],10)+string(",0;); NS=Star(rho_ci=");
    }
    else
    {
        binopts+=string("0,0,0;); NS=Star(rho_ci=");
    }
    binopts+=DoubleToStringDecimal(rhoc,16)+string("; dimi=3; EOS =");
    binopts+=mEOS;
    binopts+=string("; Center=");
    binopts+=DoubleToStringDecimal(CenterNS[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[2],16)+string("; Velocity=");
    binopts+=DoubleToStringDecimal(dta*CenterNS[0]-OmegaOrbit*CenterNS[1],10)+string(",");
    binopts+=DoubleToStringDecimal(dta*CenterNS[1]+OmegaOrbit*CenterNS[0],10)+string(",0;); Attenuation=true; AttenuationCenterNS=");
    binopts+=DoubleToStringDecimal(CenterNS[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[2],16)+string("; AttenuationCenterBH=");
    binopts+=DoubleToStringDecimal(CenterBH[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[2],16)+string("; AttenuationRadiusNS=");
    binopts+=DoubleToStringDecimal(AttenuationRadiusNS,16)+string("; AttenuationRadiusBH=");
    binopts+=DoubleToStringDecimal(AttenuationRadiusBH,16)+string("; HorizonBH=");
    binopts+=DoubleToStringDecimal(HorRadius,16)+string("; Ratio=");
    binopts+=DoubleToStringDecimal(Ratio,16)+string("; Spin=");		  
    binopts+=DoubleToStringDecimal(mSpin[0],16)+string(",");		  
    binopts+=DoubleToStringDecimal(mSpin[1],16)+string(",");
    binopts+=DoubleToStringDecimal(mSpin[2],16)+string(";");

    AttenuatedBinary binary(binopts);
    MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");
    DataMesh Zero=coords[0]*0.0;
    MyVector<DataMesh> gcoords(MV::fill,Zero,coords[0],coords[1],coords[2]);
    binary.SetCoordinates(gcoords);

    const Mesh& sd=box.Get<Mesh>("Mesh");
    Tensor<DataMesh> g(3,"11",sd);
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "11", sd);


    binary.GetLowerMetric(g);
    for(int i=0;i<3;i++)
      {
	for(int j=0;j<3;j++)
	  (*mpResult)(i,j)=g(i,j);
      }
    }


  ///////////////////////////////////////////////////////////////////////
  /////////////////////////ExtrinsicCurvature////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  BinaryExtrinsicCurvature::BinaryExtrinsicCurvature(const string& opts):
    mSpin(MV::fill,0,0,0),mpResult(0) {
    OptionParser p(opts, Help());
    mMassBH=p.Get<double>("MassBH",1.0);
    mCenterBH=p.Get<string>("CenterBH","CenterBH");
    mCenterNS=p.Get<string>("CenterNS","CenterNS");
    mOmegaOrbit=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mSpin=p.Get<MyVector<double> >("Spin");
    mdta=p.Get<string>("dtExpansionFactor","dtExpansionFactor");
    rhoc=p.Get<double>("rhoc");
    mEOS=p.Get<string>("EOS");
    AttenuationRadiusBH=p.Get<double>("AttenuationRadiusBH");
    HorRadius=p.Get<double>("HorizonBH",0.0);
    Ratio=p.Get<double>("Ratio",0.0);
    doBoost=p.Get<bool>("Boost",true);
    AttenuationRadiusNS=p.Get<double>("AttenuationRadiusNS");
    mOutput=p.Get<string>("Output");
    }
  
  void BinaryExtrinsicCurvature::RecomputeData(const DataBoxAccess& box) 
    const {

    MyVector<double> CenterBH=box.Root().Get<MyVector<double> >(mCenterBH);
    MyVector<double> CenterNS=box.Root().Get<MyVector<double> >(mCenterNS);
    double OmegaOrbit=box.Root().Get<MyVector<double> >(mOmegaOrbit)[2];
    double dta=box.Root().Get<double>(mdta);

    string binopts="BH=KerrSchild(Mass=";
    binopts+=DoubleToStringDecimal(mMassBH,16)+string(";Center=");
    binopts+=DoubleToStringDecimal(CenterBH[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[2],16)+string("; Spin=");		  
    binopts+=DoubleToStringDecimal(mSpin[0],16)+string(",");		  
    binopts+=DoubleToStringDecimal(mSpin[1],16)+string(",");
    binopts+=DoubleToStringDecimal(mSpin[2],16)+string("; Velocity=");
    if(doBoost)
    {
        binopts+=DoubleToStringDecimal(dta*CenterBH[0]-OmegaOrbit*CenterBH[1],10)+string(",");
        binopts+=DoubleToStringDecimal(dta*CenterBH[1]+OmegaOrbit*CenterBH[0],10)+string(",0;); NS=Star(rho_ci=");
    }
    else
    {
        binopts+=string("0,0,0;); NS=Star(rho_ci=");
    }
    binopts+=DoubleToStringDecimal(rhoc,16)+string("; dimi=3; EOS =");
    binopts+=mEOS;
    binopts+=string("; Center=");
    binopts+=DoubleToStringDecimal(CenterNS[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[2],16)+string("; Velocity=");
    binopts+=DoubleToStringDecimal(dta*CenterNS[0]-OmegaOrbit*CenterNS[1],10)+string(",");
    binopts+=DoubleToStringDecimal(dta*CenterNS[1]+OmegaOrbit*CenterNS[0],10)+string(",0;); Attenuation=true; AttenuationCenterNS=");
    binopts+=DoubleToStringDecimal(CenterNS[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[2],16)+string("; AttenuationCenterBH=");
    binopts+=DoubleToStringDecimal(CenterBH[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[2],16)+string("; AttenuationRadiusNS=");
    binopts+=DoubleToStringDecimal(AttenuationRadiusNS,16)+string("; AttenuationRadiusBH=");
    binopts+=DoubleToStringDecimal(AttenuationRadiusBH,16)+string("; HorizonBH=");
    binopts+=DoubleToStringDecimal(HorRadius,16)+string("; Ratio=");
    binopts+=DoubleToStringDecimal(Ratio,16)+string("; Spin=");		  
    binopts+=DoubleToStringDecimal(mSpin[0],16)+string(",");		  
    binopts+=DoubleToStringDecimal(mSpin[1],16)+string(",");
    binopts+=DoubleToStringDecimal(mSpin[2],16)+string(";");

    AttenuatedBinary binary(binopts);
    MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");
    DataMesh Zero=coords[0]*0.0;
    MyVector<DataMesh> gcoords(MV::fill,Zero,coords[0],coords[1],coords[2]);
    binary.SetCoordinates(gcoords);

    const Mesh& sd=box.Get<Mesh>("Mesh");
    Tensor<DataMesh> K(3,"11",sd);
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "11", sd);


    binary.GetLowerExCurv(K);
    for(int i=0;i<3;i++)
      {
	for(int j=0;j<3;j++)
	  //(*mpResult)(i,j)=0.0;
	  (*mpResult)(i,j)=K(i,j);
      }

    }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////////Lapse////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  BinaryLapse::BinaryLapse(const string& opts):
    mSpin(MV::fill,0,0,0),mpResult(0) {
    OptionParser p(opts, Help());
    mMassBH=p.Get<double>("MassBH",1.0);
    mCenterBH=p.Get<string>("CenterBH","CenterBH");
    mCenterNS=p.Get<string>("CenterNS","CenterNS");
    mOmegaOrbit=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mSpin=p.Get<MyVector<double> >("Spin");
    mdta=p.Get<string>("dtExpansionFactor","dtExpansionFactor");
    rhoc=p.Get<double>("rhoc");
    mEOS=p.Get<string>("EOS");
    AttenuationRadiusBH=p.Get<double>("AttenuationRadiusBH");
    HorRadius=p.Get<double>("HorizonBH",0.0);
    Ratio=p.Get<double>("Ratio",0.0);
    doBoost=p.Get<bool>("Boost",true);
    AttenuationRadiusNS=p.Get<double>("AttenuationRadiusNS");
    mOutput=p.Get<string>("Output");
    }
  
  void BinaryLapse::RecomputeData(const DataBoxAccess& box) 
    const {

    MyVector<double> CenterBH=box.Root().Get<MyVector<double> >(mCenterBH);
    MyVector<double> CenterNS=box.Root().Get<MyVector<double> >(mCenterNS);
    double OmegaOrbit=box.Root().Get<MyVector<double> >(mOmegaOrbit)[2];
    double dta=box.Root().Get<double>(mdta);

    string binopts="BH=KerrSchild(Mass=";
    binopts+=DoubleToStringDecimal(mMassBH,16)+string(";Center=");
    binopts+=DoubleToStringDecimal(CenterBH[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[2],16)+string("; Spin=");		  
    binopts+=DoubleToStringDecimal(mSpin[0],16)+string(",");		  
    binopts+=DoubleToStringDecimal(mSpin[1],16)+string(",");
    binopts+=DoubleToStringDecimal(mSpin[2],16)+string("; Velocity=");
    if(doBoost)
    {
	binopts+=DoubleToStringDecimal(dta*CenterBH[0]-OmegaOrbit*CenterBH[1],10)+string(",");
	binopts+=DoubleToStringDecimal(dta*CenterBH[1]+OmegaOrbit*CenterBH[0],10)+string(",0;); NS=Star(rho_ci=");
    }    
    else
    {
	binopts+=string("0,0,0;); NS=Star(rho_ci=");
    }
    binopts+=DoubleToStringDecimal(rhoc,16)+string("; dimi=3; EOS =");
    binopts+=mEOS;
    binopts+=string("; Center=");
    binopts+=DoubleToStringDecimal(CenterNS[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[2],16)+string("; Velocity=");
    binopts+=DoubleToStringDecimal(dta*CenterNS[0]-OmegaOrbit*CenterNS[1],10)+string(",");
    binopts+=DoubleToStringDecimal(dta*CenterNS[1]+OmegaOrbit*CenterNS[0],10)+string(",0;); Attenuation=true; AttenuationCenterNS=");
    binopts+=DoubleToStringDecimal(CenterNS[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterNS[2],16)+string("; AttenuationCenterBH=");
    binopts+=DoubleToStringDecimal(CenterBH[0],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[1],16)+string(",");
    binopts+=DoubleToStringDecimal(CenterBH[2],16)+string("; AttenuationRadiusNS=");
    binopts+=DoubleToStringDecimal(AttenuationRadiusNS,16)+string("; AttenuationRadiusBH=");
    binopts+=DoubleToStringDecimal(AttenuationRadiusBH,16)+string("; HorizonBH=");
    binopts+=DoubleToStringDecimal(HorRadius,16)+string("; Ratio=");
    binopts+=DoubleToStringDecimal(Ratio,16)+string("; Spin=");		  
    binopts+=DoubleToStringDecimal(mSpin[0],16)+string(",");		  
    binopts+=DoubleToStringDecimal(mSpin[1],16)+string(",");
    binopts+=DoubleToStringDecimal(mSpin[2],16)+string(";");

    AttenuatedBinary binary(binopts);
    MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");
    DataMesh Zero=coords[0]*0.0;
    MyVector<DataMesh> gcoords(MV::fill,Zero,coords[0],coords[1],coords[2]);
    binary.SetCoordinates(gcoords);

    const Mesh& sd=box.Get<Mesh>("Mesh");
    Tensor<DataMesh> N(3,"",sd);
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "", sd);


    binary.GetPhysicalLapse(N);
    (*mpResult)()=N();

    }


  ///////////////////////////////////////////////////////////////////////
  ///  Tensors entering the equation to solve for the hydro potential ///
  ///////////////////////////////////////////////////////////////////////

  FIDPotential_BulkGlobalCoeff::FIDPotential_BulkGlobalCoeff(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mRho=p.Get<string>("BaryonDensity","BaryonDensity");
    rhoc=p.Get<double>("CentralDensity");
    mOutput=p.Get<string>("Output");
    mPowerLaw=p.Get<double>("PowerLaw",1.0);
    }
  
  void FIDPotential_BulkGlobalCoeff::RecomputeData(const DataBoxAccess& box) 
    const {
    Dm Rho =  box.Get<TDm>(mRho)();

    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "", sd);

    (*mpResult)()=Rho/rhoc;
    
    for(int i=0;i<(*mpResult)().Size();i++)
    {
	double sign=1.0;
	if((*mpResult)()[i]<0.0)
	{
	    sign=-1.0;
	    (*mpResult)()[i]*=-1.0;
	}
	(*mpResult)()[i]=sign*pow((*mpResult)()[i],mPowerLaw);
    }
    	
    }
  
  //####################################################################


  FIDPotential_SurfaceScalarCoeff::FIDPotential_SurfaceScalarCoeff(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mShift=p.Get<string>("Shift");
    mU=p.Get<string>("InertialVelocity","InertialVelocity");
    mOmega=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mH=p.Get<string>("Enthalpy","Enthalpy");
    rhoc=p.Get<double>("CentralDensity");
    mrho=p.Get<string>("BaryonDensity","BaryonDensity");
    mdrho=p.Get<string>("dBaryonDensity","dBaryonDensity");
    mg=p.Get<string>("ConformalMetric","ConformalMetric");
    mInvg=p.Get<string>("InvConformalMetric","InvConformalMetric");
    mCenterW=p.Get<string>("CenterW","CenterW");
    mdtExpansionFactorName = p.Get<string>("dtExpansionFactor",
					   "dtExpansionFactor");
    mPowerLaw=p.Get<double>("PowerLaw",1.0);
    mRot = p.Get<string>("RotationTerm", "RotationTerm");
    }
  
  void FIDPotential_SurfaceScalarCoeff::RecomputeData(const DataBoxAccess& box) 
    const {
    Dm Psi =  box.Get<TDm>(mPsi)();
    Dm NPsi =  box.Get<TDm>(mNPsi)();
    TDm InShift =  box.Get<TDm>(mShift);

    Dm Lapse=NPsi/Psi;
    TDm U =  box.Get<TDm>(mU);
    TDm H =  box.Get<TDm>(mH);
    TDm rho =  box.Get<TDm>(mrho);
    TDm drho =  box.Get<TTDm>(mdrho)();
    TDm Rot = box.Get<TDm>(mRot);

    const GlobalDifferentiator &Deriv = box.Get<GlobalDifferentiator>
	("GlobalDifferentiator");

    if(mPowerLaw==1.0)
	for(int i=0;i<drho.Dim();i++)
	    drho(i)/=rhoc;
    else
    {
	Dm rhoalpha=rho()/rhoc;
	for(int i=0;i<rhoalpha.Size();i++)
	{
	    double sign=1.0;
	    if(rhoalpha[i]<0.0)
	    {
		sign=-1.0;
		rhoalpha[i]*=-1.0;
	    }
	    rhoalpha[i]=sign*pow(rhoalpha[i],mPowerLaw)/mPowerLaw;
	}
	Deriv.Differentiate(rhoalpha,drho);
    }

    TDm g =  box.Get<TDm>(mg);
    TDm Invg =  box.Get<TDm>(mInvg);
    MyVector<double> W=box.Root().Get<MyVector<double> >(mCenterW);

    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "", sd);
    
    //get OmegaOrbit from the databox
    const MyVector<double> Omega = box.Root().Get<MyVector<double > >(mOmega);
    MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");

    double mdtExpansionFactor = box.Root().Get<double >(mdtExpansionFactorName);
    TDm mdtExpansionFactorTimesrHi(3,"1",sd);
    mdtExpansionFactorTimesrHi(0) = mdtExpansionFactor * coords[0];
    mdtExpansionFactorTimesrHi(1) = mdtExpansionFactor * coords[1];
    mdtExpansionFactorTimesrHi(2) = mdtExpansionFactor * coords[2];

    //Omegaxr
    TDm OmegaxrHi = InShift;
    OmegaxrHi(0) = Omega[1]*coords[2] - Omega[2]*coords[1];
    OmegaxrHi(1) = Omega[2]*coords[0] - Omega[0]*coords[2];
    OmegaxrHi(2) = Omega[0]*coords[1] - Omega[1]*coords[0];

    TDm Shift=InShift;
    for(int i=0;i<3;i++) Shift(i)+=OmegaxrHi(i)
			       +mdtExpansionFactorTimesrHi(i);

    Dm temp=H()*0.0+1.0;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	temp-=g(i,j)*sqr(sqr(Psi))*U(i)*U(j);
	}
    }
    for(int p=0;p<temp.Size();p++)
      if(temp[p]>0.)
	temp[p]=1.0/sqrt(temp[p]);
      else
	temp[p]=1.;

    TDm resint(3,"1",sd);

    for(int i=0;i<3;i++)
      {
	resint(i)=Shift(i)*H()*sqr(sqr(Psi))*temp/Lapse;
	for(int j=0;j<3;j++)
	  resint(i)-=Invg(i,j)*(W[j]+Rot(j));
	for(int j=0;j<H().Size();j++)
	  {
	    if(!(fabs(resint(i)[j])>0))
	      resint(i)[j]=0.0;
	  }
      }

    
    
    (*mpResult)()=0.0;
    for(int i=0;i<3;i++)
      (*mpResult)()-=resint(i)*drho(i);

    }

  
  //####################################################################

  FIDPotential_SurfaceFirstDerivCoeff::FIDPotential_SurfaceFirstDerivCoeff(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mH=p.Get<string>("Enthalpy","Enthalpy");
    rhoc=p.Get<double>("CentralDensity");
    mrho=p.Get<string>("BaryonDensity","BaryonDensity");
    mdrho=p.Get<string>("dBaryonDensity","dBaryonDensity");
    mPowerLaw=p.Get<double>("PowerLaw",1.0);
    mInvg=p.Get<string>("InvConformalMetric","InvConformalMetric");
    }
  
  void FIDPotential_SurfaceFirstDerivCoeff::RecomputeData(const DataBoxAccess& box) 
    const {
    TDm H =  box.Get<TDm>(mH);
    TDm Invg =  box.Get<TDm>(mInvg);

    TDm rho =  box.Get<TDm>(mrho);
    TDm drho =  box.Get<TTDm>(mdrho)();

    const GlobalDifferentiator &Deriv = box.Get<GlobalDifferentiator>
        ("GlobalDifferentiator");

    if(mPowerLaw==1.0)
        for(int i=0;i<drho.Dim();i++)
            drho(i)/=rhoc;
    else
    {
        Dm rhoalpha=rho()/rhoc;
        for(int i=0;i<rhoalpha.Size();i++)
	{
	    double sign=1.0;
            if(rhoalpha[i]<0.0)
            {
                sign=-1.0;
                rhoalpha[i]*=-1.0;
            }
            rhoalpha[i]=sign*pow(rhoalpha[i],mPowerLaw)/mPowerLaw;
	}
        Deriv.Differentiate(rhoalpha,drho);
    }


    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "1", sd);
    
    MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");
    
    for(int i=0;i<3;i++){
      (*mpResult)(i)=0.0;
      for(int j=0;j<3;j++)
	(*mpResult)(i)+=drho(j)*Invg(i,j);
    }
    }
  
  //####################################################################

  FIDPotential_BulkFirstDerivCoeff::FIDPotential_BulkFirstDerivCoeff(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mGamma=p.Get<string>("ConformalChristoffel2ndKind","ConformalChristoffel2ndKind");
    mH=p.Get<string>("Enthalpy","Enthalpy");
    mInvg=p.Get<string>("InvConformalMetric","InvConformalMetric");
    }
  
  void FIDPotential_BulkFirstDerivCoeff::RecomputeData(const DataBoxAccess& box) 
    const {
    Dm Psi =  box.Get<TDm>(mPsi)();
    Dm NPsi =  box.Get<TDm>(mNPsi)();

    Dm Lapse=NPsi/Psi;
    TDm H =  box.Get<TDm>(mH);
    TDm Invg =  box.Get<TDm>(mInvg);
    TDm Gamma =  box.Get<TDm>(mGamma);


    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "1", sd);
    
    Dm temp=Lapse*sqr(Psi)/H();

    const GlobalDifferentiator &Deriv = box.Get<GlobalDifferentiator>
                                               ("GlobalDifferentiator");
    
    TDm dtemp(3,"1",sd);
    Deriv.Differentiate(temp,dtemp);

    for(int i=0;i<3;i++)
      {
	dtemp(i)/=temp;
      }

    for(int i=0;i<3;i++)
      {
	(*mpResult)(i)=0.0;
	for(int j=0;j<3;j++)
	  {
	    (*mpResult)(i)-=Invg(i,j)*dtemp(j);
	    for(int k=0;k<3;k++)
	      (*mpResult)(i)+=Invg(j,k)*Gamma(i,j,k);
	  }
      }
    }
    

  //####################################################################

  FIDPotential_BulkScalarCoeff::FIDPotential_BulkScalarCoeff(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mShift=p.Get<string>("Shift");
    mH=p.Get<string>("Enthalpy","Enthalpy");
    mg=p.Get<string>("ConformalMetric","ConformalMetric");
    mU=p.Get<string>("InertialVelocity","InertialVelocity");
    mOmega=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mTrK=p.Get<string>("TrExtrinsicCurvature","TrExtrinsicCurvature");
    mCoeffC=p.Get<string>("BulkFirstDerivCoeff","BulkFirstDerivCoeff");
    mCenterW=p.Get<string>("CenterW","CenterW");
    mdtExpansionFactorName = p.Get<string>("dtExpansionFactor",
					   "dtExpansionFactor");
    mRot = p.Get<string>("RotationTerm", "RotationTerm");
    mInvg = p.Get<string>("InverseConformalMetric", "InverseConformalMetric");
   
    }
  
  void FIDPotential_BulkScalarCoeff::RecomputeData(const DataBoxAccess& box) 
    const {
    TDm Rot = box.Get<TDm>(mRot);
    TDm Invg = box.Get<TDm>(mInvg);
    Dm Psi =  box.Get<TDm>(mPsi)();
    Dm NPsi =  box.Get<TDm>(mNPsi)();
    TDm InShift =  box.Get<TDm>(mShift);

    Dm Lapse=NPsi/Psi;
    Dm H =  box.Get<TDm>(mH)();
    TDm U =  box.Get<TDm>(mU);
    TDm g =  box.Get<TDm>(mg);
    TDm A =  box.Get<TDm>(mCoeffC);
    MyVector<double> W=box.Root().Get<MyVector<double> >(mCenterW);
    Dm TrK=  box.Get<TDm>(mTrK)();
    


    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "", sd);
     
    //get OmegaOrbit from the databox
    const MyVector<double> Omega = box.Root().Get<MyVector<double > >(mOmega);
    const MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");

    double mdtExpansionFactor = box.Root().Get<double >(mdtExpansionFactorName);
    TDm mdtExpansionFactorTimesrHi(3,"1",sd,0.);
    mdtExpansionFactorTimesrHi(0) = mdtExpansionFactor * coords[0];
    mdtExpansionFactorTimesrHi(1) = mdtExpansionFactor * coords[1];
    mdtExpansionFactorTimesrHi(2) = mdtExpansionFactor * coords[2];

    //Omegaxr
    TDm OmegaxrHi = InShift;
    OmegaxrHi(0) = Omega[1]*coords[2] - Omega[2]*coords[1];
    OmegaxrHi(1) = Omega[2]*coords[0] - Omega[0]*coords[2];
    OmegaxrHi(2) = Omega[0]*coords[1] - Omega[1]*coords[0];

    TDm Shift=InShift;
    for(int i=0;i<3;i++) Shift(i)+=OmegaxrHi(i)
			       +mdtExpansionFactorTimesrHi(i);

    Dm gn(sd);
    gn=1.0;
    for(int i=0;i<3;i++)
      {
	for(int j=0;j<3;j++)
	  {
	    gn-=g(i,j)*U(i)*U(j)*sqr(sqr(Psi));
	  }
      }
    
    for(int p=0;p<gn.Size();p++)
      if(gn[p]>0.)
        gn[p]=1.0/sqrt(gn[p]);
      else
	gn[p]=1.;

    const GlobalDifferentiator &Deriv = box.Get<GlobalDifferentiator>
                                               ("GlobalDifferentiator");
    

    Tensor<TDm> dRot(3,"1",Rot);
    Deriv.Differentiate(Rot,dRot);
    
    TDm dgn(3,"1",sd);
    Deriv.Differentiate(gn,dgn);

    (*mpResult)()=H*TrK*gn*sqr(sqr(Psi));
	
    for(int i=0;i<3;i++)
      {
	(*mpResult)()+=H*sqr(sqr(Psi))/Lapse*Shift(i)*dgn(i);
      }
	
    for(int i=0;i<3;i++)
      {
	(*mpResult)()+=(W[i]+Rot(i))*A(i);
      }

    for (int i =0; i<3; i++){
      for (int j=0; j<3; j++){
	(*mpResult)() -= Invg(i,j)*dRot(i)(j);
      }
    }
    
    for(int j=0;j<H.Size();j++)
      {
        if(!(fabs((*mpResult)()[j])>0))
          (*mpResult)()[j]=0.0;
      }
    }
  

  //####################################################################

  TrKTimeInvariant::TrKTimeInvariant(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mShift=p.Get<string>("Shift");
    mOmega=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mdtExpansionFactorName = p.Get<string>("dtExpansionFactor",
					   "dtExpansionFactor");
    AttRad=p.Get<double>("AttenuationRadius");
    mGamma=p.Get<string>("ConformalChristoffel2ndKind","ConformalChristoffel2ndKind");
    }
  
  void TrKTimeInvariant::RecomputeData(const DataBoxAccess& box) 
    const {
    Dm Psi =  box.Get<TDm>(mPsi)();
    Dm NPsi =  box.Get<TDm>(mNPsi)();
    TDm InShift =  box.Get<TDm>(mShift);
    TDm Gamma=box.Get<TDm>(mGamma);

    Dm Lapse=NPsi/Psi;

    const Mesh& sd=box.Get<Mesh>("Mesh");
    if(mpResult==0) 
      mpResult = new Tensor<DataMesh>(3, "", sd);
    

    //get OmegaOrbit from the databox
    const MyVector<double> Omega = box.Root().Get<MyVector<double > >(mOmega);
    const MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");

    double mdtExpansionFactor = box.Root().Get<double >(mdtExpansionFactorName);
    TDm mdtExpansionFactorTimesrHi(3,"1",sd,0.);
    mdtExpansionFactorTimesrHi(0) = mdtExpansionFactor * coords[0];
    mdtExpansionFactorTimesrHi(1) = mdtExpansionFactor * coords[1];
    mdtExpansionFactorTimesrHi(2) = mdtExpansionFactor * coords[2];

    //Omegaxr
    TDm OmegaxrHi = InShift;
    OmegaxrHi(0) = Omega[1]*coords[2] - Omega[2]*coords[1];
    OmegaxrHi(1) = Omega[2]*coords[0] - Omega[0]*coords[2];
    OmegaxrHi(2) = Omega[0]*coords[1] - Omega[1]*coords[0];

    TDm Shift=InShift;
    for(int i=0;i<3;i++) Shift(i)+=OmegaxrHi(i)
			       +mdtExpansionFactorTimesrHi(i);

    const GlobalDifferentiator &Deriv = box.Get<GlobalDifferentiator>
                                               ("GlobalDifferentiator");
    
    TDm dtemp(3,"1",sd);

    Deriv.Differentiate(Psi,dtemp);

    (*mpResult)()=0.0;
    for(int i=0;i<3;i++)
      (*mpResult)()+=Shift(i)*dtemp(i);
    (*mpResult)()*=6.0/Psi;

    for(int i=0;i<3;i++)
      {
	Deriv.Differentiate(InShift(i),dtemp);
	(*mpResult)()+=dtemp(i);
	for(int k=0;k<3;k++)
	  (*mpResult)()+=Shift(i)*Gamma(k,i,k);
      }
    (*mpResult)()+=3.0*mdtExpansionFactor;
    (*mpResult)()/=Lapse;

    }
  //####################################################################

  dTrKTimeInvariant::dTrKTimeInvariant(const string& opts):
    mpResult(0) {
    OptionParser p(opts, Help());
    mOutput=p.Get<string>("Output");
    mPsi=p.Get<string>("ConformalFactor");
    mNPsi=p.Get<string>("LapseTimesConformalFactor");
    mShift=p.Get<string>("Shift");
    mK=p.Get<string>("TrExtrinsicCurvature","TrExtrinsicCurvature");
    mOmega=p.Get<string>("OmegaOrbit","OmegaOrbit");
    mdtExpansionFactorName = p.Get<string>("dtExpansionFactor",
					   "dtExpansionFactor");
    AttRad=p.Get<double>("AttenuationRadius");
    }
  
  void dTrKTimeInvariant::RecomputeData(const DataBoxAccess& box) 
    const {
    Dm Psi =  box.Get<TDm>(mPsi)();
    Dm NPsi =  box.Get<TDm>(mNPsi)();
    TDm InShift =  box.Get<TDm>(mShift);
    Dm K =  box.Get<TDm>(mK)();

    Dm Lapse=NPsi/Psi;

    const Mesh& sd=box.Get<Mesh>("Mesh");
    const TDm Ts(3, "1", sd);
    if(mpResult==0) 
      mpResult = new Tensor<Tensor<DataMesh> >(3, "", Ts);
     
    //get OmegaOrbit from the databox
    const MyVector<double> Omega = box.Root().Get<MyVector<double > >(mOmega);
    const MyVector<DataMesh> coords=box.Get<MyVector<DataMesh> >("GlobalCoords");

    double mdtExpansionFactor = box.Root().Get<double >(mdtExpansionFactorName);
    TDm mdtExpansionFactorTimesrHi(3,"1",sd,0.);
    mdtExpansionFactorTimesrHi(0) = mdtExpansionFactor * coords[0];
    mdtExpansionFactorTimesrHi(1) = mdtExpansionFactor * coords[1];
    mdtExpansionFactorTimesrHi(2) = mdtExpansionFactor * coords[2];

    //Omegaxr
    TDm OmegaxrHi = InShift;
    OmegaxrHi(0) = Omega[1]*coords[2] - Omega[2]*coords[1];
    OmegaxrHi(1) = Omega[2]*coords[0] - Omega[0]*coords[2];
    OmegaxrHi(2) = Omega[0]*coords[1] - Omega[1]*coords[0];

    TDm Shift=InShift;
    for(int i=0;i<3;i++) Shift(i)+=OmegaxrHi(i)
			       +mdtExpansionFactorTimesrHi(i);

    const GlobalDifferentiator &Deriv = box.Get<GlobalDifferentiator>
                                               ("GlobalDifferentiator");
    
    Dm temp(sd);
    TDm dtemp(3,"1",sd);

    Deriv.Differentiate(Lapse,dtemp);
    for(int i=0;i<3;i++)
      {
	(*mpResult)()(i)=-dtemp(i)*K;
      }
    
    temp=0.0;

    for(int i=0;i<3;i++)
      {
	Deriv.Differentiate(InShift(i),dtemp);
	temp+=dtemp(i);
      }
    Deriv.Differentiate(temp,dtemp);
    
    for(int i=0;i<3;i++)
      {
	(*mpResult)()(i)+=dtemp(i);
      }

    temp=0.0;

    for(int i=0;i<3;i++)
      {
	Deriv.Differentiate(Psi,dtemp);
	temp+=6.0*InShift(i)*dtemp(i)/Psi;
      }
    Deriv.Differentiate(temp,dtemp);
    
    for(int i=0;i<3;i++)
      {
	(*mpResult)()(i)+=dtemp(i);
      }
    
    Deriv.Differentiate(Psi,dtemp);
    for(int i=0;i<3;i++)
      {
	(*mpResult)()(i)+=6.0*mdtExpansionFactor*dtemp(i)/Psi;
      }

    (*mpResult)()(0)+=Omega[2]*dtemp(1)*6.0/Psi;
    (*mpResult)()(1)-=Omega[2]*dtemp(0)*6.0/Psi;

    for(int i=0;i<3;i++)
      {
	(*mpResult)()(i)/=Lapse;
      }

    Dm rad(sd);
    rad=0.0;
    for(int i=0;i<3;i++)
      rad+=sqr(coords[i]);
    rad=sqrt(rad);
    Dm mult=exp(-sqr(sqr(rad/AttRad)));
    for(int i=0;i<3;i++)
      {
	(*mpResult)()(i)*=mult;
	(*mpResult)()(i)-=4.0*sqr(rad)/sqr(sqr(AttRad))*coords[i]*K*mult;
      }
    }
  

}//namespace


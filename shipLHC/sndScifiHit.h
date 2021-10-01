#ifndef SNDSCIFIHIT_H
#define SNDSCIFIHIT_H 1

#include "SndlhcHit.h"
#include "ScifiPoint.h"
#include "TObject.h"
#include "TVector3.h"

class sndScifiHit : public SndlhcHit
{
  public:

    /** Default constructor **/
    sndScifiHit();
    sndScifiHit(Int_t detID);
    // make hit
    void makeHit(int detID,std::vector<ScifiPoint*>,std::vector<Float_t>);

 /** Destructor **/
    virtual ~sndScifiHit();

    /** Output to screen **/
	void Print() ;
	Float_t GetEnergy();

	void setInvalid() {flag = false;}
	bool isValid() const {return flag;}
	Int_t GetStation(){return int(fDetectorID/10000000);}
	bool isVertical(){
		if ( int(fDetectorID/100000)%10 == 1){return true;}
		else{return (false);}
	}
	Int_t  GetMat(){return ( int(fDetectorID/100000)%10);}
	Int_t  GetSiPM(){return ( int(fDetectorID/10000)%10);}
	Int_t GetSiPMChan(){return ( fDetectorID%1000);}
	Int_t GetChannelID(){return fDetectorID;}
	void SetThreshold(Float_t x){nphe_min = x;}
	Float_t GetThreshold(){return nphe_min;}
/*  
	from Guido (22.9.2021): A threshold of 3.5pe should be used, which corresponds to 0.031MeV.
	1 SiPM channel has 104 pixels, pixel can only see 0 or >0 photons.
*/
	Float_t nphe_min;
	Float_t nphe_max = 104;
  private:
    /** Copy constructor **/
    sndScifiHit(const sndScifiHit& hit);
    sndScifiHit operator=(const sndScifiHit& hit);
    Float_t ly_loss_mean(Float_t distance, Float_t* params);
    Float_t MeanAndRMS(Double_t ly);
    Float_t flag;   ///< flag

// get parameters from the Scifi detector for simulating the digitized information
	Float_t ly_loss_params[4] = {20.78, -0.26, 7.89, -3.06};

    ClassDef(sndScifiHit,1);

};

#endif
#include <ADC.h>
#include <ADC_util.h>
#include <SPI.h> 

#include "filterFunctions.h"
#include "BiquadFilter.h"

/*
 * Feedback v3
 * Do not consider resonant speaker. Compensation filter quality
 * factor is adjusted to stay under 1/sqrt(2).
 * v3.2 uses filter function
 */


//===================================================================
// Speaker TS parameters - Peerless

const float Mms = 10.0e-3; //#[kg]
const float Rms = 1.462; //#[kg/s]
const float Cms = 590e-6; //#[mm/N]
const float Bl  = 5.594; //#[T.m]
const float Re  = 7.0; //#[ohm]
const float Qs  = 1/(Rms+pow(Bl,2)/Re)*sqrt(Mms/Cms); //speaker quality factor

//Global variables
const float Fs = 48000.0;
float Xmax = 0.125e-3; // Maximum displacement definition (Xmax) [m]
float Gain = 0.5f;  // Amplifier gain + Teensy gain
volatile float x;    // Estimated displacement [m]
const float attack_smooth  = 0.02; // Attack..
const float release_smooth = 0.5;  // ... and release time for gain smoothing [s]
volatile float out = 0;    // Output voltage

// X/U filter
float b_xu[3], a_xu[3];   //analog  coeffs of X/U TF
float bd_xu[3], ad_xu[3]; //digital coeffs of X/U TF

BiquadFilterDF1 xuFilter; //tension to displacement (X/U) estimator


// Compensation filter & indicators
float b_comp[3], a_comp[3];   //analog  coeffs of Comp filter denumerator
float bd_comp[3], ad_comp[3]; //digital coeffs of Comp filter denumerator

BiquadFilterTDF2 CompFilter; //Compensation filter

const float C_min = min(1, Xmax*Re/(Gain*Bl*Cms)); // optimal minimal compliance to respect Xmax
float C_instant = 1.0f; // instantaneous ratio on compliance (Cms_comp/Cms)

const float Cms_min = 0.9*C_min*Cms;   // minimum compensation compliance with a margin factor
float Cms_comp = Cms; // variable Cms_compensation
float Cms_target = Cms; // Cms target to vary

const float Q0 = 0.707; // neutral quality factor
volatile float Qc = Qs;

volatile float Rms_comp = Rms; // Rms compensation to be adujsted real-time
float R_instant = 1.0f;        // instantaneous ratio of Rms (Rms/Rms_comp)


// Gain smoothing
const float attack_coeff  = 1 - exp(-2.2 / (attack_smooth  * Fs));
const float release_coeff = 1 - exp(-2.2 / (release_smooth * Fs));
float k;

//===================================================================
/*
 * IN1: A14
 * IN2: A15
 * IN3: A12
 * IN4: A13
 */
const int readPin_In1 = A14;
const int readPin_In3 = A12;
const int readPin_In4 = A13;

const int outPin_Aux1 = A21;
const int outPin_Aux2 = A22;

const int resolutionADC = 16;
const int resolutionDAC = 16; //MAX 5717 DAC
const int resolutionDACT = 12; //Teensy DAC

const float VrefADC = 3.3; // Teensy ADC VRef
const float VrefDACT = 3.3; // Teensy DAC VRef
const float VrefDAC = 4.096;  //Modified from 10.0 19/04/22

const int slaveSelectPin = 9;
const int LDAC = 10;

const int DIO = 2;

const float conversionConstADC = VrefADC/((1<<resolutionADC)-1);
const float conversionConstDAC = ((1<<resolutionDAC)-1)/(2.0*VrefDAC);
const float conversionConstDACT = ((1<<resolutionDACT)-1)/VrefDACT;
const int ADCAverages = 0;


//Timing
const float dt = 1.0/Fs; 
const float SPIMHz = 22.5;
const long SPIHz = SPIMHz*1000*1000;

//Input (ADC) initialisation
volatile uint16_t valADC0 = 0, valADC1 = 0;
volatile uint8_t VALADC0Ready = false, VALADC1Ready = false;
volatile float input1 = 0.0, input3 = 0.0, input4 = 0.0;


//Output (DAC) signal initialisation
volatile uint16_t valDAC = 0, valDACT = 0, valDACT2 = 0;
volatile float val4DACOut = 0.0;

// Using the Pedvide ADC Library on Github
ADC *adc = new ADC(); // adc object



void setup(void)
{
  // Declare needed pins as inputs or outputs
  pinMode(readPin_In1, INPUT);
  pinMode(readPin_In4, INPUT);
  pinMode (slaveSelectPin, OUTPUT);
  pinMode(DIO, OUTPUT);
  pinMode (LDAC, OUTPUT);

  /**
   * Here we initialise the ADCs
   */

  // When using the Pedvide ADC library we can set more options
  adc->adc0->setAveraging(ADCAverages); // set number of averages
  adc->adc0->setResolution(resolutionADC); // set bits of resolution
  adc->adc0->setReference(ADC_REFERENCE::REF_3V3); // Set voltage reference for ADC.
  adc->adc0->setSamplingSpeed(ADC_SAMPLING_SPEED::MED_SPEED); // change the sampling speed
  adc->adc0->setConversionSpeed(ADC_CONVERSION_SPEED::MED_SPEED); // change the conversion speed

  adc->adc0->stopPDB();
  adc->adc0->startSingleRead(readPin_In1); // call this to setup everything before the pdb starts, differential is also possible
  adc->adc0->enableInterrupts(adc0_isr);
  adc->adc0->startPDB(Fs); //frequency in Hz
  
  adc->adc1->setAveraging(ADCAverages); // set number of averages
  adc->adc1->setResolution(resolutionADC); // set bits of resolution
  adc->adc1->setReference(ADC_REFERENCE::REF_3V3); // Set voltage reference for ADC.
  adc->adc1->setSamplingSpeed(ADC_SAMPLING_SPEED::MED_SPEED); // change the sampling speed
  adc->adc1->setConversionSpeed(ADC_CONVERSION_SPEED::MED_SPEED); // change the conversion speed

  adc->adc1->stopPDB();
  adc->adc1->startSingleRead(readPin_In4); // call this to setup everything before the pdb starts, differential is also possible
  adc->adc1->enableInterrupts(adc1_isr);
  adc->adc1->startPDB(Fs); //frequency in Hz

  
  NVIC_SET_PRIORITY(IRQ_USBOTG, 200);
  
  // Here we initialise the SPI channel
  SPI.begin();
  SPI.beginTransaction(SPISettings(SPIHz, MSBFIRST, SPI_MODE0));

  // Serial monitor/plotter
  Serial.begin(115200);
  Serial.println("Compliance max ratio:");
  Serial.println(C_min, 3);

  // Here we set up the onboard DACsand resolution definition Teensy DACs
  analogWriteResolution(12);

  //===================================================================
  //===================================================================
  setXUSpeakerCoefficients(Re, Bl, Rms, Mms, Cms, Fs, bd_xu, ad_xu);
  setCompFilterCoefficients(Re, Bl, Rms, Rms_comp, Mms, Cms, Cms_comp, Fs, bd_comp, ad_comp);

  xuFilter.setCoefficients(bd_xu, ad_xu);
  CompFilter.setCoefficients(bd_comp, ad_comp);

  // Serial.println("coeff xu - b then a");
  // Serial.println(bd_xu, 10); Serial.println(ad_xu, 10);
  // Serial.println("Comp coeff - b then a");
  // Serial.println(bd_comp, 10); Serial.println(ad_comp, 10);
}



void Operations(void)
{
  // SPI output is done first.
  // This ensure the DAC output timing is synchronised with the ADC via the PDB (Programmable Delay Block).
  uint8_t SPIBuff[2] = {0};   //Initialise the SPI buffer.  Initialiser la buffer SPI.
  SPIBuff[0] = valDAC >> 8;   // Bit Splitting. 16 bit -> |8 bit||8 bit|
  SPIBuff[1] = valDAC & 0xFF; // Separation des bits fort et faible. 16 bit -> |8 bit||8 bit|

  // Sending data to MAX5717 Via SPI. See Data sheet for details
  digitalWrite(LDAC,HIGH);    
  // take the SS pin low to select the chip:
  digitalWrite(slaveSelectPin,LOW);
  //  send in the address and value via SPI:
  SPI.transfer(SPIBuff, 2);
  // take the SS pin high to de-select the chip:
  digitalWrite(slaveSelectPin,HIGH);
  digitalWrite(LDAC,LOW);

  //=======================Core algorithm here=========================
  //===================================================================
   
  //Read ADC Value and convert to voltage while removing any offsets
  input1 = (valADC0 * conversionConstADC - 1.65);
  // input4 = (valADC1 * conversionConstADC - 1.65);  

  // Estimate the displacement considering amp. gain and monitor it
  x = xuFilter.process(out*Gain);
  
  // valDACT2 = ((x*1000 + VrefDACT * 0.5) * conversionConstDACT);   //DAC Output - Teensy 12 bit
  // analogWrite(outPin_Aux2, valDACT2);

  // Compute Cms compensation value
  if (abs(x) > Xmax){
    Cms_target = Cms_min;
  } else {
    Cms_target = Cms;
  }

  if (Cms_target < Cms_comp) {
    k = attack_coeff;
  } else {
    k = release_coeff;
  }

  // Cms_comp = (1-k)*Cms_comp + k*Cms_target;
  Cms_comp = Cms;

  // Compute new equivalent quality factor of Hxu*Hcomp
  Qc = 1/(Rms+pow(Bl,2)/Re)*sqrt(Mms/Cms_comp);

  // Adjust Rms comp
  // if (Qc > Q0){
  //   Rms_comp = 1/Q0*sqrt(Mms/Cms_comp) - (Bl*Bl/Re);
  // } else {
  //   Rms_comp = Rms;
  // }

  // For test
  Rms_comp = Rms;

  // Monitor compliance ratio (C_intant is between 0 and 1)
  C_instant = Cms_comp/Cms;
  valDACT  = ((C_instant + VrefDACT * 0.5) * conversionConstDACT);    //DAC Output - Teensy 12 bit
  analogWrite(outPin_Aux1, valDACT);
  // Serial.println("Compliance instantaneous ratio:");
  // Serial.println(C_instant, 3);
  
  // and monitor ratio of Rms
  R_instant = Rms/Rms_comp;
  valDACT2  = ((R_instant + VrefDACT * 0.5) * conversionConstDACT);    //DAC Output - Teensy 12 bit
  analogWrite(outPin_Aux2, valDACT2); 

  // Update filter coefficient (A_2_comp only for now)
  setCompFilterCoefficients(Re, Bl, Rms, Rms_comp, Mms, Cms, Cms_comp, Fs, bd_comp, ad_comp);

  CompFilter.setCoefficients(bd_comp, ad_comp);

  // volatile float A_2_comp = 1/Cms_comp;

  // volatile float a_0_comp = A_0*4*pow(Fs,2) + A_1*2*Fs + A_2_comp;
  // volatile float a_1_comp = -2*A_0*4*pow(Fs,2) + 2*A_2_comp;
  // volatile float a_2_comp = A_0*4*pow(Fs,2) - A_1*2*Fs + A_2_comp;


  //Hcomp coefficients normalization by a_0_comp 
  // b_0_comp = b_0_comp/a_0_comp;
  // b_1_comp = b_1_comp/a_0_comp;
  // b_2_comp = b_2_comp/a_0_comp;
  // a_1_comp = a_1_comp/a_0_comp;
  // a_2_comp = a_2_comp/a_0_comp;
  // a_0_comp = 1.0f;

  // Apply filter
  out = CompFilter.process(input1);

  // val4DACOut = b_0_comp*input1 + d0_comp;

  // Update filter states
  // d0_comp = b_1_comp*input1 - a_1_comp*val4DACOut + d1_comp;
  // d1_comp = b_2_comp*input1 - a_2_comp*val4DACOut;

  val4DACOut = x*1e3;

  valDAC = (uint16_t)((val4DACOut + VrefDAC)*conversionConstDAC);

}



void loop(void)
{
  if (VALADC0Ready && VALADC1Ready){
    Operations();
    VALADC0Ready = false;
    VALADC1Ready = false;
  } 
}



void adc0_isr() {
    valADC0 = (uint16_t)adc->adc0->readSingle();
    VALADC0Ready = true;
}

void adc1_isr() {
    valADC1 = (uint16_t)adc->adc1->readSingle();
    VALADC1Ready = true;
}
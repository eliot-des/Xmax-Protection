#include <ADC.h>
#include <ADC_util.h>
#include <SPI.h> 

#include "filterFunctions.h"
#include "BiquadFilter.h"
#include "DelayLine.h"

// V3 with only Cms_comp + Rms_comp adjustment

//===================================================================
//Peerless SDS-P830656 driver params
const float fr  = 65.5;     //[Hz]
const float Re  = 7.00;      //[Ohm]
const float Bl  = 5.595;     //[N/A]
const float Rms = 1.4497;    //[N.s/m]
const float Mms = 10.0e-3;  //[Kg]
const float Cms = 595e-6;   //[N/m]
const float Qs  = 1/(Rms+pow(Bl,2)/Re)*sqrt(Mms/Cms); //speaker quality factor

//Global variables

BiquadFilterDF1 xuFilter;     //tension to displacement (X/U) estimator
BiquadFilterTDF2 CompFilter;  //Compensation filter
BiquadFilterTDF2 CompFilterDelayed;  //Compensation filter delayed
DelayLine delayLine;



const float Fs = 48000.0;
const float Xmax = 10.0e-3; // Maximum displacement definition (Xmax) [m]
const float Gain = 16.33f+0.1f;  // Amplifier gain
volatile float x = 0.0f;   // Estimated displacement [m]
const float attack_smooth  = 0.003; // 0.02; // Attack time for gain smoothing
const float release_smooth = 0.5;  // Release time for gain smoothing

const int nAttack        = 2.0f * Fs * attack_smooth;

volatile float outDelayed = 0; 
volatile float out = 0;    // Output voltage

// X/U filter
float b_xu[3], a_xu[3];   //analog  coeffs of X/U TF

// Compensation filter & indicators
float b_comp[3], a_comp[3];   //analog  coeffs of Comp filter denumerator


const float C_min = min(1, Xmax*Re*2.0/(Gain*Bl*Cms)); // optimal minimal compliance to respect Xmax
float C_instant = C_min;                           // instantaneous ratio on compliance (Cms_comp/Cms)

const float Cms_min = 0.9*C_min*Cms;    // minimum compensation compliance with a margin factor
float Cms_comp      = Cms;              // variable Cms_compensation
float Cms_target    = Cms;              // Cms target to vary

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
// volatile float val4DACOut = 0.0;

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
  Serial.println("Speaker Qs:");
  Serial.println(Qs, 3);

  // Here we set up the onboard DACsand resolution definition Teensy DACs
  analogWriteResolution(12);

  //===================================================================
  //===================================================================
  setXUSpeakerCoefficients(Re, Bl, Rms, Mms, Cms, Fs, b_xu, a_xu);
  updateCompensationCoefficients(Re, Bl, Rms, Mms, Cms, Rms_comp, Cms_comp, Fs, b_comp, a_comp);

  xuFilter.setCoefficients(b_xu, a_xu);
  CompFilter.setCoefficients(b_comp, a_comp);
  CompFilterDelayed.setCoefficients(b_comp, a_comp);


  delayLine.setMaximumDelayInSamples(nAttack);
  delayLine.reset();
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
  delayLine.write(input1);

  // Estimate the displacement considering amp. gain and monitor it
  x = xuFilter.process(out*Gain);

  // Compute Cms compensation value
  if (fabs(x) > Xmax){
    Cms_target = Cms_min;
  } else {
    Cms_target = Cms;
  }

  //smoothing process
  if (Cms_target < Cms_comp) { 
    k = attack_coeff;
  } else { 
    k = release_coeff;
  }
  Cms_comp = (1-k)*Cms_comp + k*Cms_target;

  // Compute new equivalent quality factor of Hxu*Hcomp
  Qc = 1/(Rms+pow(Bl,2)/Re)*sqrt(Mms/Cms_comp);

  // Adjust Rms comp
  if (Qc > Q0){
    Rms_comp = 1/Q0*sqrt(Mms/Cms_comp) - (Bl*Bl/Re);
  } else {
    Rms_comp = Rms;
  }

  // Update compensation filter coefficient (bd_comp, ad_comp)
  updateCompensationCoefficients(Re, Bl, Rms, Mms, Cms, Rms_comp, Cms_comp, Fs, b_comp, a_comp);
  CompFilter.setCoefficients(b_comp, a_comp);
  CompFilterDelayed.setCoefficients(b_comp, a_comp);

  // Apply filter
  out = CompFilter.process(input1);
  outDelayed = CompFilterDelayed.process(delayLine.read(nAttack));
  //val4DACOut = input4;

  valDAC = (uint16_t)((outDelayed + VrefDAC)*conversionConstDAC);

  valDACT  = ((Cms_comp/Cms + VrefDACT * 0.5) * conversionConstDACT);    //DAC Output - Teensy 12 bit
  valDACT2  = ((x*500 + VrefDACT * 0.5) * conversionConstDACT);    //DAC Output - Teensy 12 bit # x*500 pour loud music tests Serie 1
  
  analogWrite(outPin_Aux1, valDACT);
  analogWrite(outPin_Aux2, valDACT2); 
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
#include <ADC.h>
#include <ADC_util.h>
#include <SPI.h>

#include "LS_Model.h"

LS_Model lsModel_p;
LS_Model lsModel_np;


/*
 * Speaker TS parameters, Xmax definition & chain gain
 */

// full-range1 parameters
const float Mms = 2.72e-3; //#[kg]
const float Rms = 0.368907; //#[kg/s]
const float Cms = 0.61e-3; //#[mm/N]
const float Bl = 5.17; //#[T.m]
const float Re = 7.3; //#[ohm]
const float Qs = 1/(Rms+pow(Bl,2)/Re)*sqrt(Mms/Cms); //speaker quality factor

// Maximum displacement definition (Xmax)
const float Xmax = 1.5e-3; //#[m]
// Amplifier + Teensy board Gain
const float Gain = 8;


/*
 * Timing
 */
const float sampleRateHz = 48000.0;
const float dt = 1.0/sampleRateHz; 
const float SPIMHz = 22.5;
const long SPIHz = SPIMHz*1000*1000;


/*
 * X/U IIR filter
 */

// analog coefficients
const float B_0 = 0;
const float B_1 = 0;
const float B_2 = Bl/Re;
const float A_0 = Mms;
const float A_1 = Rms+pow(Bl,2)/Re;
const float A_2 = 1/Cms;

// digital coefficients from analytical bilinear transform
volatile float b_0 = B_2;
volatile float b_1 = 2*B_2;
volatile float b_2 = B_2;
volatile float a_0 = A_0*4*pow(sampleRateHz,2) + A_1*2*sampleRateHz + A_2;
volatile float a_1 = -2*A_0*4*pow(sampleRateHz,2) + 2*A_2;
volatile float a_2 = A_0*4*pow(sampleRateHz,2) - A_1*2*sampleRateHz + A_2;

// copying for compensation filter - step right after
volatile float b_0_comp = a_0;
volatile float b_1_comp = a_1;
volatile float b_2_comp = a_2;

// //Coefficients are normalized with a_0 -> not working like that, need to be in the 
// b_0 = b_0/a_0;
// b_1 = b_0/a_0;
// b_2 = b_0/a_0;
// a_1 = b_0/a_0;
// a_2 = b_0/a_0;
// a_0 = 1.0f;

// DFI buffer (useless?)
// volatile float d0_hp = 0, d1_hp = 0, d2_hp = 0, d3_hp = 0;


/*
 * Compensation filter
 */

// Optimal minimal compliance to respect Xmax
const float R = Xmax*Re/(Gain*Bl*Cms);
const float Cms_min = 0.9*R*Cms; //minimum compensation compliance with a margin factor

// Initialization
float Cms_comp = Cms;
float Cms_target;
float Rms_comp = Rms;
const float Q0 = 1/sqrt(2) ; // neutral quality factor
float Qc;

// Coefficients of the compensation filter
// B_i_comp coefficients are the A_i ones...
// float A_0_comp = A_0;
float A_1_comp = Rms_comp+pow(Bl,2)/Re;
float A_2_comp = 1/Cms_comp;

// digital coefficients of the initial compensation filter
volatile float a_0_comp = A_0*4*pow(sampleRateHz,2) + A_1_comp*2*sampleRateHz + A_2_comp;
volatile float a_1_comp = -2*A_0*4*pow(sampleRateHz,2) + 2*A_2_comp;
volatile float a_2_comp = A_0*4*pow(sampleRateHz,2) - A_1_comp*2*sampleRateHz + A_2_comp;

// //Coefficients are normalized with a_0_comp
// b_0_comp /= a_0_comp;
// b_1_comp /= a_0_comp;
// b_2_comp /= a_0_comp;
// a_1_comp /= a_0_comp;
// a_2_comp /= a_0_comp;
// a_0_comp = 1.0f;

// TDFII buffer (compensation filter)
volatile float d0_comp = 0, d1_comp = 0;


/*
 * Gain smoothing setting
 */

// Attack and release time in seconds for the gain smoothing
const float attack_smooth  = 0.02;
const float release_smooth = 0.5;
const float attack_coeff  = 1 - exp(-2.2 / (attack_smooth  * sampleRateHz));
const float release_coeff = 1 - exp(-2.2 / (release_smooth * sampleRateHz));
float k;

/*
 * Outputs
 */

// Estimated displacement with protection
volatile float x_est;
// Estimated displacemet without protection
volatile float x_est_nc;
// Protected signal - limited voltage signal sent to the speaker
volatile float out;

/*
 * First sample deal
 */

bool itsFirst = true;
// bool itsOn = true;
// int incomignByte = 0;
// bool modRelease = false;
// bool modAttack = false;

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


/**
 * Initialisation I/O
 */
//Input (ADC) initialisation
volatile uint16_t valADC0 = 0, valADC1 = 0;
volatile uint8_t VALADC0Ready = false, VALADC1Ready = false;
volatile float input1 = 0.0, input3 = 0.0, input4 = 0.0;


//Output (DAC) signal initialisation
volatile uint16_t valDAC = 0, valDACT = 0;
volatile float val4DACOut = 0.0;


// Using the Pedvide ADC Library on Github
ADC *adc = new ADC(); // adc object

void setup(void)
{
  // Normalize coefficients
  b_0 = b_0/a_0;
  b_1 = b_1/a_0;
  b_2 = b_2/a_0;
  a_1 = a_1/a_0;
  a_2 = a_2/a_0;
  a_0 = 1.0f;

  //Coefficients are normalized with a_0_comp
  b_0_comp /= a_0_comp;
  b_1_comp /= a_0_comp;
  b_2_comp /= a_0_comp;
  a_1_comp /= a_0_comp;
  a_2_comp /= a_0_comp;
  a_0_comp = 1.0f;

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
  adc->adc0->startPDB(sampleRateHz); //frequency in Hz
  
  adc->adc1->setAveraging(ADCAverages); // set number of averages
  adc->adc1->setResolution(resolutionADC); // set bits of resolution
  adc->adc1->setReference(ADC_REFERENCE::REF_3V3); // Set voltage reference for ADC.
  adc->adc1->setSamplingSpeed(ADC_SAMPLING_SPEED::MED_SPEED); // change the sampling speed
  adc->adc1->setConversionSpeed(ADC_CONVERSION_SPEED::MED_SPEED); // change the conversion speed

  adc->adc1->stopPDB();
  adc->adc1->startSingleRead(readPin_In4); // call this to setup everything before the pdb starts, differential is also possible
  adc->adc1->enableInterrupts(adc1_isr);
  adc->adc1->startPDB(sampleRateHz); //frequency in Hz

  
  NVIC_SET_PRIORITY(IRQ_USBOTG, 200);
  
  /**
   * Here we initialise the SPI channel
   */

  SPI.begin();
  SPI.beginTransaction(SPISettings(SPIHz, MSBFIRST, SPI_MODE0));

  /**
   * Incase you need to use the serial monitor/plotter
   */
  //Serial.begin(9600);

  /**
   * Here we set up the onboard DACs
   * 
   * Definition de la resolution des DAC Teensy (si besoin)
   */
  //analogWriteResolution(12);

}


void Operations(void)
{
  /**
   * SPI output is done first.
   * This ensure the DAC output timing is synchronised with the ADC via the PDB.
   * 
   * La sortie via SPI est fait en premier.
   * C'est pour synchroniser la sortie DAC SPI avec l'ADC. Le frequence d'échantillonage du DAC et ADC est donc gerer par le PDB (Programmable Delay Block)
   */
  uint8_t SPIBuff[2] = {0};   //Initialise the SPI buffer.  Initialiser la buffer SPI.
  SPIBuff[0] = valDAC >> 8;   // Bit Splitting. 16 bit -> |8 bit||8 bit|
  SPIBuff[1] = valDAC & 0xFF; // Separation des bits fort et faible. 16 bit -> |8 bit||8 bit|

  /**
   * Sending data to MAX5717 Via SPI. See Data sheet for details
   */
  digitalWrite(LDAC,HIGH);    
  // take the SS pin low to select the chip:
  digitalWrite(slaveSelectPin,LOW);
  //  send in the address and value via SPI:
  SPI.transfer(SPIBuff, 2);
  // take the SS pin high to de-select the chip:
  digitalWrite(slaveSelectPin,HIGH);
  digitalWrite(LDAC,LOW);


  /**
   * Here we read the ADC values and convert into a voltage, while removing any offsets
   */
  
  //Read ADC Value and convert to voltage
  input1 = (valADC0 * conversionConstADC - 1.65);
  // input4 = (valADC1 * conversionConstADC - 1.65);
  //Serial.println(input2);
  /**
   * Preparing value for DAC Output - SPI
   */

  if (itsFirst) {
    // First displacement sample is calculated with the input signal (input1)
    x_est = lsModel_p.x_predictor(b_0, b_1, b_2, a_1, a_2, input1, Gain);
    itsFirst = false; }
  else {
    // Other displacement samples are calculated with the filtered signal (out)
    x_est = lsModel_p.x_predictor(b_0, b_1, b_2, a_1, a_2, out, Gain); }
  

  // Estimation of the not corrected displacement
  // Not necessary but useful for comparing displacement estimation with the original signal.
  // x_est_nc = lsModel_np.x_predictor(b_0, b_1, b_2, a_1, a_2, input1, Gain);

  /*
  Feedback algorithm
  */

  // compare the peak of the displacement signal with Xmax
  if (abs(x_est) > Xmax) {
    Cms_target = Cms_min; }
  else {
    Cms_target = Cms; }    

  // we then apply the gain smoothing function to Cms_comp
  if (Cms_target < Cms_comp) {
    k = attack_coeff; }
  else {
    k = release_coeff; }
    
  Cms_comp = (1 - k) * Cms_comp + k * Cms_target;
  Cms_comp = 
  // Quality factor for the compensation filter - c stands for compensation
  Qc = 1/(Rms+pow(Bl,2)/Re)*sqrt(Mms/Cms_comp);

  // Modifiy Rms if nessecary to get Qc lower than Q0 (1/sqrt(2))
  if (Qc > Q0) {
    Rms_comp = 1/Q0*sqrt(Mms/Cms_comp) - (pow(Bl,2)/Re); }
  else {
    Rms_comp = Rms; }

  // Update compensation filter analog coefficient,
  A_1_comp = Rms_comp+pow(Bl,2)/Re;
  A_2_comp = 1/Cms_comp;

  // discretize them,
  b_0_comp = a_0;
  b_1_comp = a_1;
  b_2_comp = a_2;
  a_0_comp = A_0*4*pow(sampleRateHz,2) + A_1_comp*2*sampleRateHz + A_2_comp;
  a_1_comp = -2*A_0*4*pow(sampleRateHz,2) + 2*A_2_comp;
  a_2_comp = A_0*4*pow(sampleRateHz,2) - A_1_comp*2*sampleRateHz + A_2_comp;

  // and normalize by a_0_comp.
  b_0_comp /= a_0_comp;
  b_1_comp /= a_0_comp;
  b_2_comp /= a_0_comp;
  a_1_comp /= a_0_comp;
  a_2_comp /= a_0_comp;
  a_0_comp /= a_0_comp;

  // Transposed Direct Form II (TDFII) is used to apply the compensation filter
  out = b_0_comp*input1 + d0_comp;

  // Update buffer
  d0_comp = b_1_comp*input1 - a_1_comp*out + d1_comp;
  d1_comp = b_2_comp*input1 - a_2_comp*out;



  // Limiter voltage
  val4DACOut = input1;
  //val4DACOut = input4;
  valDAC = (uint16_t)((val4DACOut + VrefDAC)*conversionConstDAC);

  //Serial.println(valDAC);

  /*
   * DAC Output - Teensy 12 bit
   */
  //valDACT = ((val4DACOut + VrefDACT*0.5)*conversionConstDACT);
  
  //analogWrite(outPin_Aux1, valDACT);
  //analogWrite(outPin_Aux2, valDACT);
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

// pdb interrupt is enabled in case you need it.
//void pdb_isr(void) {
//        PDB0_SC &=~PDB_SC_PDBIF; // clear interrupt
//        //digitalWriteFast(LED_BUILTIN, !digitalReadFast(LED_BUILTIN) );
//}
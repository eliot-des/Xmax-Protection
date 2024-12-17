#include <ADC.h>
#include <ADC_util.h>
#include <Arduino.h>
#include <SPI.h>

#include "filterFunctions.h"
#include "DelayLine.h"
#include "MovingMax.h"
#include "BiquadFilter.h"

//===================================================================
//global variables

DelayLine delayLine;
DelayLine gainLine;

BiquadFilterDF1 xuFilter;           //tension to displacement (X/U) estimator
BiquadFilterTDF2 lowShelfFilter;    //adaptative low-shelf filter

const float Fs = 48000.0;
float b_xu[3], a_xu[3];       //digital coeffs of X/U TF
float b_shelf[3], a_shelf[3]; //digital coeffs of adaptative low shelf filter

//linear gain amplification ratio between output of the DAC and loudspeaker
const float ampFactor = 10.0f; 

//Peerless SDS-P830656 driver params
const float fr  = 72.0;     //[Hz]
const float Rec = 6.4;      //[Ohm]
const float Bl  = 5.74;     //[N/A]
const float Rms = 1.897;    //[N.s/m]
const float Mms = 8.88e-3;  //[Kg]
const float Cms = 560e-6;   //[N/m]


const float threshold   = 0.5e-3;  // Xmax threshold [mm]
const float attackTime  = 0.003;  // [s]
const float holdTime    = 0.002;  // [s]
const float releaseTime = 0.1;    // [s]


const int nAttack        = Fs * attackTime;
const int nAttackHold    = Fs * (attackTime+holdTime);
const float releaseCoeff = 1 - exp(-2.2 / (Fs * releaseTime));


MovingMax<float, nAttackHold> movingMax; // Create an instance of MovingMax


volatile float x = 0.0;  //predicted displacement
volatile float out = 0.0;

//side-chain variables
volatile float gainComp = 0.0; //gain computer output
volatile float gainMin  = 0.0; //gainComp filtered by the minimum filter
volatile float gainRel  = 0.0; //gainMin passing through release block
volatile float gainRelPrev = 0.0;
volatile float g = 0.0;        //gainRel passing through averaging filter -> gain applied to delayed input

//low-shelf coefficients
volatile float shelfGain = 0.0;//gain in dB applied to the low-shelf filter
const float shelfQ  = 0.707;     //Q factor of the shelf filter
const float Fc = fr;
const float wc = 2 * PI * Fc / Fs;
const float alpha = sin(wc) / (2 * shelfQ);
const float beta = cos(wc);

float lastShelfGain = -60.0; // Initialize to a default value

//===================================================================
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
const float VrefDAC = 4.096; //Modified from 10.0 19/04/22

const int slaveSelectPin = 9;
const int LDAC = 10;

const int DIO = 2;

const float conversionConstADC = VrefADC / ((1 << resolutionADC) - 1);
const float conversionConstDAC = ((1 << resolutionDAC) - 1) / (2.0 * VrefDAC);
const float conversionConstDACT = ((1 << resolutionDACT) - 1) / VrefDACT;
const int ADCAverages = 0;


//Timing
const float dt = 1.0 / Fs;
const float SPIMHz = 22.5;
const long SPIHz = SPIMHz * 1000 * 1000;

//Input (ADC) initialisation
volatile uint16_t valADC0 = 0, valADC1 = 0;
volatile uint8_t VALADC0Ready = false, VALADC1Ready = false;
volatile float input1 = 0.0, input3 = 0.0, input4 = 0.0;

//Output (DAC) signal initialisation
volatile uint16_t valDAC = 0, valDACT = 0, valDACT2 = 0;

// Using the Pedvide ADC Library on Github
ADC *adc = new ADC(); // adc object



void setup(void) {
    // Declare needed pins as inputs or outputs
    pinMode(readPin_In1, INPUT);
    pinMode(readPin_In4, INPUT);
    pinMode(slaveSelectPin, OUTPUT);
    pinMode(DIO, OUTPUT);
    pinMode(LDAC, OUTPUT);

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
    adc->adc0->startPDB(Fs); //frequency in Hz 165

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

    //Initialise the SPI channel
    SPI.begin();
    SPI.beginTransaction(SPISettings(SPIHz, MSBFIRST, SPI_MODE0));

    //Use serial plotter
    Serial.begin(115200);

    //Set up the onboard DACs -> define resolution
    analogWriteResolution(12);

    Serial.println("Begin of processing");


    //===================================================================
    //===================================================================
    setXUSpeakerCoefficients(Rec, Bl, Rms, Mms, Cms, Fs, b_xu, a_xu);
    setLowShelfCoefficients(shelfGain);

    xuFilter.setCoefficients(b_xu, a_xu);
    lowShelfFilter.setCoefficients(b_shelf, a_shelf);

    delayLine.setMaximumDelayInSamples(nAttack);
    delayLine.reset();

    gainLine.setMaximumDelayInSamples(nAttack);
    gainLine.reset();

    Serial.println("Setup done.");
    Serial.println(nAttack);
    Serial.println(nAttackHold);
}



void Operations(void) {
    //===================================================================
    // SPI output is done first.
    // This ensure the DAC output timing is synchronised with the ADC via the PDB (Programmable Delay Block)
    uint8_t SPIBuff[2] = { 0 }; //Initialise the SPI buffer. Initialiser la buffer SPI.
    SPIBuff[0] = valDAC >> 8; // Bit Splitting. 16 bit -> |8 bit||8 bit|
    SPIBuff[1] = valDAC & 0xFF; // Separation des bits fort et faible. 16 bit -> |8 bit||8 bit|

    //Sending data to MAX5717 Via SPI. See Data sheet for details
    digitalWrite(LDAC, HIGH);
    // take the SS pin low to select the chip:
    digitalWrite(slaveSelectPin, LOW);
    // send in the address and value via SPI:
    SPI.transfer(SPIBuff, 2);
    // take the SS pin high to de-select the chip:
    digitalWrite(slaveSelectPin, HIGH);
    digitalWrite(LDAC, LOW);
    //===================================================================
    //===================================================================
    // CORE OF THE ALGORITHM

    //Read ADC Value and convert to voltage while removing any offsets
    input1 = (valADC0 * conversionConstADC - 1.65);
    
    //Put the input signal in the delay line
    delayLine.write(input1);
  
    //Estimate the displacement from the given voltage
    x = xuFilter.process(input1);
    
    //Estimate the gain that should be apply on the low-shelf to not exceed the threshold displacement
    gainComp = min(1.0, threshold/abs(x * ampFactor));
  
    //Apply minimum moving filter
    gainMin = -movingMax.add(-gainComp);

    //Apply the release
    gainRel = min(gainMin, (1.0 - releaseCoeff) * gainRel + releaseCoeff * gainMin);
    if (gainRel > 0.99) gainRel = 1.0;
    
    //Average filtering with recursive implementation
    gainRelPrev = gainLine.read(nAttack);
    gainLine.write(gainRel);
    g = g + (gainRel - gainRelPrev)/nAttack;

    //Set adaptative low-shelf coefficients
    shelfGain = 20.0 * log10(g);
    if (shelfGain != lastShelfGain) {
      setLowShelfCoefficients(shelfGain);
      lowShelfFilter.setCoefficients(b_shelf, a_shelf);
      lastShelfGain = shelfGain;
    }


    //apply the low-shelf filter on the delayed input
    out = lowShelfFilter.process(delayLine.read(nAttack));
    //out = input1;
    //value for the "true" DAC output
    valDAC = (uint16_t)((out + VrefDAC) * conversionConstDAC);

    //values to tracks can be monitored using valDACT, valDACT2 (use auxilary DAC)
    valDACT  = ((gainMin + VrefDACT * 0.5) * conversionConstDACT);    //DAC Output - Teensy 12 bit
    valDACT2 = ((g + VrefDACT * 0.5) * conversionConstDACT);   //DAC Output - Teensy 12 bit
    analogWrite(outPin_Aux1, valDACT);
    analogWrite(outPin_Aux2, valDACT2);
}

void loop(void) {
    if (VALADC0Ready && VALADC1Ready) {
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



// Function to compute low-shelf filter coefficients
void setLowShelfCoefficients(float dBgain) {
    float A = pow(10, dBgain / 40);

    float sqrtA = sqrt(A);
    float Aplus1 = A + 1;
    float Aminus1 = A - 1;


    b_shelf[0] = A * (Aplus1 - Aminus1 * beta + 2 * sqrtA * alpha);
    b_shelf[1] = 2 * A * (Aminus1 - Aplus1 * beta);
    b_shelf[2] = A * (Aplus1 - Aminus1 * beta - 2 * sqrtA * alpha);

    a_shelf[0] = Aplus1 + Aminus1 * beta + 2 * sqrtA * alpha;
    a_shelf[1] = -2 * (Aminus1 + Aplus1 * beta);
    a_shelf[2] = Aplus1 + Aminus1 * beta - 2 * sqrtA * alpha;

    normalize(b_shelf, a_shelf);
}



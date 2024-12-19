#include <ADC.h>
#include <ADC_util.h>
#include <SPI.h>

#include "LS_Model.h"

LS_Model lsModel_p;
LS_Model lsModel_np;



// ---------------------------------------------------
// Loudspeaker
// ---------------------------------------------------
// #Parameters - Basic Linear Model
// At low frequencies Le is approximately 0

// Previous parameters
// const float Mms = 8.897e-3; //#[Kg]
// const float Rms = 1.75; //#[N.s/m]
// const float Kms = 1.73 * 1000; //#[N/m]
// const float Bl = 5.97; //#[N/A]
// const float Re = 5.72; //#[ohm]
// const float Le = 0.29e-3; // [H]

// Peerless Loudspeaker
const float Re  = 6.4;      //[Ohm]
const float Bl  = 5.74;     //[N/A]
const float Rms = 1.897;    //[N.s/m]
const float Mms = 8.88e-3;  //[Kg]
const float Cms = 560e-6;   //[N/m]
const float Kms = 1/Cms;

// Audax Loudspeaker
// const float Mms = 20.4e-3; //#[Kg]
// const float Rms = 1.379; //#[N.s/m]
// const float Kms = 520; //#[N/m]
// const float Bl = 9.2; //#[N/A]
// const float Re = 6.06; //#[ohm]

// Amplifier + Teensy board Gain
float Gain = 8.0f;
// Maximum displacement definition (Xmax)
float Xmax = 0.5e-3;
// Estimated displacement with protection
volatile float x_est;
// Estimated displacemet without protection
volatile float x_est_nc;
// Protected signal, output of the adaptive filter
volatile float out;

///////////////////////////
// Filter initialization //
///////////////////////////

// Linear coefficients of the discretized transfer
// function of the loudspeaker for implementing the
// IIR filter
float b_0, b_1, b_2, a_0, a_1, a_2;

// Coefficients for the adaptive high-pass filter
volatile float b0_hp, b1_hp, b2_hp, a0_hp, a1_hp, a2_hp;

// Adaptive filter cut-off frequency initialization
volatile float fc = 5;

// Attack and release time in seconds for the adaptive filter.
// tao_a: Attack time, tao_r: Release time
float tao_a = 0.01, tao_r = 0.1;
// float tao_a = 0.001, tao_r = 0.01;

// Maximum and minimum cut-off frequency for the adaptive filter
float fcmax = 200;
float fcmin = 5;

// High-pass filter buffer
volatile float d0_hp = 0, d1_hp = 0, d2_hp = 0, d3_hp = 0;

bool itsOn = true, itsFirst = true;
int incomignByte = 0;
bool modRelease = false;
bool modAttack = false;

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
const float VrefDAC = 4.096; //Modified from 10.0 19/04/22

const int slaveSelectPin = 9;
const int LDAC = 10;

const int DIO = 2;

const float conversionConstADC = VrefADC / ((1 << resolutionADC) - 1);
const float conversionConstDAC = ((1 << resolutionDAC) - 1) / (2.0 * VrefDAC);
const float conversionConstDACT = ((1 << resolutionDACT) - 1) / VrefDACT;
const int ADCAverages = 0;


/*
* Timing
*/
const float sampleRateHz = 48000.0;
const float dt = 1.0 / sampleRateHz;
const float SPIMHz = 22.5;
const long SPIHz = SPIMHz * 1000 * 1000;


/**
* Initialisation I/O
*/
//Input (ADC) initialisation
volatile uint16_t valADC0 = 0, valADC1 = 0;
volatile uint8_t VALADC0Ready = false, VALADC1Ready = false;
volatile float input1 = 0.0, input3 = 0.0, input4 = 0.0;


//Output (DAC) signal initialisation
volatile uint16_t valDAC = 0, valDACT = 0, valDACT2 = 0;
volatile float val4DACOut = 0.0;

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
    adc->adc0->startPDB(sampleRateHz); //frequency in Hz 165

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
    Serial.begin(115200);

    // Here the loudspeaker filter model coefficients are computed
    Bilinear_coeff(Re, Bl, Mms, Kms, Rms, sampleRateHz);

    /**
    * Here we set up the onboard DACs
    *
    * Definition de la resolution des DAC Teensy (si besoin)
    */
    analogWriteResolution(12);
}


 void Operations(void) {
    /**
    * SPI output is done first.
    * This ensure the DAC output timing is synchronised with the ADC via the PDB.
    *
    * La sortie via SPI est fait en premier.
    * C’est pour synchroniser la sortie DAC SPI avec l’ADC. Le frequence d’echantillonage du DAC et ADC est donc gerer par le PDB (Programmable Delay Block)
    */
    uint8_t SPIBuff[2] = { 0 }; //Initialise the SPI buffer. Initialiser la buffer SPI.
    SPIBuff[0] = valDAC >> 8; // Bit Splitting. 16 bit -> |8 bit||8 bit|
    SPIBuff[1] = valDAC & 0xFF; // Separation des bits fort et faible. 16 bit -> |8 bit||8 bit|

    /**
    * Sending data to MAX5717 Via SPI. See Data sheet for details
    */
    digitalWrite(LDAC, HIGH);
    // take the SS pin low to select the chip:
    digitalWrite(slaveSelectPin, LOW);
    // send in the address and value via SPI:
    SPI.transfer(SPIBuff, 2);
    // take the SS pin high to de-select the chip:
    digitalWrite(slaveSelectPin, HIGH);
    digitalWrite(LDAC, LOW);


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
      itsFirst = false;
    } else {
      // Other displacement samples are calculated with the filtered signal (out)
      x_est = lsModel_p.x_predictor(b_0, b_1, b_2, a_1, a_2, out, Gain);
    }

    // Estimation of the not corrected displacement

    // Not necessary but useful for comparing displacement estimation with
    // the original signal.
    x_est_nc = lsModel_np.x_predictor(b_0, b_1, b_2, a_1, a_2, input1, Gain);
    
    // For changing the cut-off frequency of the adaptive filter
    DynaFreq(x_est, Xmax);
    // Here is where the signal is filtered
    out = High_Pass_BW(fc, sampleRateHz, input1);
    
    // This is for monitoring the cut-off frequency through a digital output
    valDACT2 = ((fc / 200 + VrefDACT * 0.5) * conversionConstDACT);

    if (itsOn) {
      // Activate protection

      // To output estimated displacement according to the output
      valDACT = ((x_est * 500 + VrefDACT * 0.5) * conversionConstDACT);
      val4DACOut = out;
    } else {
      // Bypass protection

      // To output the estimated displacement without protection
      valDACT = ((x_est_nc * 500 + VrefDACT * 0.5) * conversionConstDACT);
      val4DACOut = input1;
    }

    valDAC = (uint16_t)((val4DACOut + VrefDAC) * conversionConstDAC);

    /**
    * DAC Output - Teensy 12 bit
    */
    // valDACT = ((val4DACOut + VrefDACT*0.5)*conversionConstDACT);


    analogWrite(outPin_Aux1, valDACT);
    analogWrite(outPin_Aux2, valDACT2);
}

void loop(void) {
    if (VALADC0Ready && VALADC1Ready) {
        Operations();
        VALADC0Ready = false;
        VALADC1Ready = false;

        /////////////////////////////////////////////////////////////
        // For controlling parameters through serial communication //
        /////////////////////////////////////////////////////////////

        // Send serial monitor messages to modify some parameters:
        // p: To activate Xmax protection
        // n: To desactivate protection
        // a: To modify attack time
        // After sending a, send the attack value in seconds
        // r: To modigy release time
        // After sending r, send the release value in seconds

        if (Serial.available() > 0) {
            if (modRelease) {
                float val = Serial.readString(4).toFloat();
                if (val != 0) {
                    tao_r = val;
                    Serial.print("Release time: ");
                    Serial.print(val, 3);
                    Serial.println(" s");
                    modRelease = false;
                }
            } else if (modAttack) {
                float val = Serial.readString(5).toFloat();
                if (val != 0) {
                    tao_a = val;
                    Serial.print("Attack time: ");
                    Serial.print(val, 4);
                    Serial.println(" s");
                    modAttack = false;
                }
            }
        
            incomignByte = Serial.read();
            if ((char)incomignByte == 'p' && (!modAttack && !modRelease)) {
                itsOn = true;
                Serial.println("Protection on");
            } else if ((char)incomignByte == 'n' && (!modAttack && !modRelease)) {
                itsOn = false;
                Serial.println("Protection off");
            } else if ((char)incomignByte == 'r' && !modRelease) {
                modRelease = true;
                Serial.println("Type the new release time in seconds");
            } else if ((char)incomignByte == 'a' && !modAttack) {
                modAttack = true;
                Serial.println("Type the new attack time in seconds");
            }
        }
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
//      PDB0_SC &=~PDB_SC_PDBIF; // clear interrupt
//      //digitalWriteFast(LED_BUILTIN, !digitalReadFast(LED_BUILTIN) );
//}

void Bilinear_coeff(float Re, float Bl, float Mms, float Kms, float Rms, float fs) {
    // This function discretizes the loudspeaker transfer function using the bilinear transform

    // Linear model assuming Le = 0
    // A and B are the analog coefficients of the transfer function of the loudspeaker
    // a and b are the discrete coefficients

    float Ts = 1 / fs;
    float B_0 = 1;
    float A_0 = (Re * Mms) / Bl;
    float A_1 = (Re * Rms + pow(Bl, 2)) / Bl;
    float A_2 = Re * Kms / Bl;

    /////////////////////////////
    // Bilinear discretization //
    // //
    // s = 2*fs*(z-1)/(z+1) //
    /////////////////////////////

    b_0 = B_0;
    b_1 = 2 * B_0;
    b_2 = B_0;

    a_0 = A_2 + 2 * A_1 / Ts + 4 * A_0 / pow(Ts, 2);
    a_1 = 2 * A_2 - 8 * A_0 / pow(Ts, 2);
    a_2 = A_2 - 2 * A_1 / Ts + 4 * A_0 / pow(Ts, 2);

    // Here all the coefficients are normalized with a_0
    b_0 /= a_0;
    b_1 /= a_0;
    b_2 /= a_0;

    a_1 /= a_0;
    a_2 /= a_0;
    a_0 /= a_0;
};

// '''Compute cut-off frequency for the adaptive filter'''
void DynaFreq(float x_est, float Xmax) {
    float M = .7; // Xmax safety margin
    float Na = tao_a * sampleRateHz;
    float Nr = tao_r * sampleRateHz;
    if (abs(x_est) > M * Xmax) {
        fc = (((fc - fcmin) / (fcmax - fcmin) - 1) * exp(-1.0 / Na) + 1) * (fcmax - fcmin) + fcmin;
    // Attack
    } else {
        fc = ((fc - fcmin) / (fcmax - fcmin) * exp(-1.0 / Nr) * (fcmax - fcmin) + fcmin);
    // Release
    }
}

///////////////////////////////
// Adaptive High-pass Filter //
///////////////////////////////
float High_Pass_BW(float fc, float fs, float input) {
    // This a second-order Butterworth high-pass filter with frequency warping
    float y;
    float wc = 2 * PI * fc / fs;
    float Q = 1 / sqrt(2);
    float alpha = sin(wc) / (2 * Q);

    b0_hp = (1 + cos(wc)) / 2;
    b1_hp = -1 - cos(wc);
    b2_hp = (1 + cos(wc)) / 2;
    a0_hp = 1 + alpha;
    a1_hp = -2 * cos(wc);
    a2_hp = 1 - alpha;

    // Coefficients normalization
    b0_hp /= a0_hp;
    b1_hp /= a0_hp;
    b2_hp /= a0_hp;

    a1_hp /= a0_hp;
    a2_hp /= a0_hp;
    a0_hp /= a0_hp;

    y = b0_hp * input + b1_hp * d0_hp + b2_hp * d1_hp - a1_hp * d2_hp - a2_hp * d3_hp;

    // Update filter states
    d1_hp = d0_hp;
    d0_hp = input;
    d3_hp = d2_hp;
    d2_hp = y;

    return y;
};
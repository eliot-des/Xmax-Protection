#include <ADC.h>
#include <ADC_util.h>
#include <SPI.h>  


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


/*
 * Timing
 */
const float sampleRateHz = 48000.0;
const float dt = 1.0/sampleRateHz; 
const float SPIMHz = 22.5;
const long SPIHz = SPIMHz*1000*1000;


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
  input4 = (valADC1 * conversionConstADC - 1.65);
  //Serial.println(input2);
  
  /**
   * Preparing value for DAC Output - SPI
   */  
  val4DACOut = input1;
  valDAC = (uint16_t)((val4DACOut + VrefDAC)*conversionConstDAC);

  //Serial.println(valDAC);
  /**
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


#ifndef LS_Model_H
#define LS_Model_H

#include <Arduino.h>

class LS_Model {
private:
    volatile float d0 = 0, d1 = 0, d2 = 0, d3 = 0;

public:
    volatile float x_predictor(float b_0, float b_1, float b_2, float a_1, float a_2, float input1, float Gain);

};

#endif
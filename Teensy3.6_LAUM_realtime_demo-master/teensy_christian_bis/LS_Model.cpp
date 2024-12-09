
#include "LS_Model.h"

// x_predictor is a method of this class that can estimate, using an IIR    filter structure (direct form I),
// the loudspeaker displacement based on the a and b coefficientes of the discretized transfer function,
// the input signal, and the gain.

volatile float LS_Model::x_predictor(float b_0, float b_1, float b_2, float a_1, float a_2, float input1, float Gain) {

    // Serial.print(input1);
    float y;
    // Apply filter
    input1 *= Gain;
    y = b_0 * input1 + b_1 * d0 + b_2 * d1 - a_1 * d2 - a_2 * d3;

    // Update filter states
    d1 = d0;
    d0 = input1;
    d3 = d2;
    d2 = y;

    return y;
}
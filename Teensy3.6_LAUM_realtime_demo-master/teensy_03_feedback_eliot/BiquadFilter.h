#ifndef BIQUADFILTER_H
#define BIQUADFILTER_H

#include <Arduino.h> // For compatibility with Teensy and Arduino


class BiquadFilterDF1 {
public:
    // Sets the coefficients for the filter
    void setCoefficients(const float* b, const float* a) {
        b0 = b[0];
        b1 = b[1];
        b2 = b[2];
        a0 = a[0];
        a1 = a[1];
        a2 = a[2];
    }

    // Processes a single float through the filter
    float process(float x) {
        float y = b0 * x + b1 * d0 + b2 * d1 - a1 * d2 - a2 * d3;

        d1 = d0;
        d0 = x;
        d3 = d2;
        d2 = y;

        return y;
    }

    // Resets the internal state of the filter
    void reset() {
        d0 = d1 = d2 = d3 = 0;
    }

private:
    float a0 = 1, a1 = 0, a2 = 0;
    float b0 = 1, b1 = 0, b2 = 0;
    float d0 = 0, d1 = 0, d2 = 0, d3 = 0; // Delay elements
};


class BiquadFilterTDF2 {
public:
    // Sets the coefficients for the filter
    void setCoefficients(const float* b, const float* a) {
        b0 = b[0];
        b1 = b[1];
        b2 = b[2];
        a0 = a[0];
        a1 = a[1];
        a2 = a[2];
    }

    // Processes a single float through the filter
    float process(float x) {
        float y = b0 * x + d0;

        d0 = b1 * x - a1 * y + d1;
        d1 = b2 * x - a2 * y;

        return y;
    }

    // Resets the internal state of the filter
    void reset() {
        d0 = d1 = 0;
    }

private:
    float a0 = 1, a1 = 0, a2 = 0; // Default coefficients for an identity filter
    float b0 = 1, b1 = 0, b2 = 0;
    float d0 = 0, d1 = 0; // Delay elements
};

#endif // BIQUADFILTER_H

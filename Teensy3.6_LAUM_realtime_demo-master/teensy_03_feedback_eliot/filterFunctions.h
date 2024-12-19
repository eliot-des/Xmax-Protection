#ifndef FILTERFUNCTIONS_H
#define FILTERFUNCTIONS_H

#include <cmath> // Include math functions like sin, pow, and sqrt

// Function to normalize coefficients
inline void normalize(float* bd, float* ad) {
    float ad0 = ad[0];
    bd[0] /= ad0;
    bd[1] /= ad0;
    bd[2] /= ad0;
    ad[0] /= ad0;
    ad[1] /= ad0;
    ad[2] /= ad0;
}

// Function to compute bilinear transform for second-order coefficients
inline void bilinear2ndOrder(float* b, float* a, float Fs, float* bd, float* ad) {
    bd[0] = b[0] * 4 * Fs * Fs + b[1] * 2 * Fs + b[2];
    bd[1] = -2 * b[0] * 4 * Fs * Fs + 2 * b[2];
    bd[2] = b[0] * 4 * Fs * Fs - b[1] * 2 * Fs + b[2];

    ad[0] = a[0] * 4 * Fs * Fs + a[1] * 2 * Fs + a[2];
    ad[1] = -2 * a[0] * 4 * Fs * Fs + 2 * a[2];
    ad[2] = a[0] * 4 * Fs * Fs - a[1] * 2 * Fs + a[2];

    normalize(bd, ad);
}

// Function to set X/U speaker coefficients
inline void setXUSpeakerCoefficients(float Rec, float Bl, float Rms, float Mms, float Cms, float Fs, float* bd_xu, float* ad_xu) {
    float b_xu[3] = { 0.0, 0.0, Bl / Rec };
    float a_xu[3] = { Mms, Rms + (Bl * Bl) / Rec, 1.0 / Cms };

    bilinear2ndOrder(b_xu, a_xu, Fs, bd_xu, ad_xu);
}

// Function to update compensation filter coefficients
inline void updateCompensationCoefficients(
    float Re, float Bl, float Rms, float Mms, float Cms, float Cms_comp,
    float Fs, float* b_comp, float* a_comp) 
{
    float B_comp[3] = {Mms, Rms + (Bl * Bl) / Re, 1.0 / Cms};
    float A_comp[3] = {Mms, Rms + (Bl * Bl) / Re, 1.0 / Cms_comp};
  
    bilinear2ndOrder(B_comp, A_comp, Fs, b_comp, a_comp);
}

#endif // FILTERFUNCTIONS_H

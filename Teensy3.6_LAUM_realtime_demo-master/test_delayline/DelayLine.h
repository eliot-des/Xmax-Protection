#ifndef DELAYLINE_H
#define DELAYLINE_H

#include <Arduino.h> // Include Arduino core for Teensy compatibility

class DelayLine {
public:
    // Constructor
    DelayLine() {}

    // Sets the maximum delay length in samples
    void setMaximumDelayInSamples(int maxLengthInSamples) {
        if (maxLengthInSamples > 0) {
            // Add 2 samples for future possible interpolation (can be omitted if not needed)
            int paddedLength = maxLengthInSamples + 2;
            if (bufferLength < paddedLength) {
                bufferLength = paddedLength;
                delete[] buffer; // Free old memory if it exists
                buffer = new float[bufferLength]; // Allocate new memory
                reset();
            }
        }
    }

    // Resets the delay line
    void reset() {
        writeIndex = bufferLength - 1;
        if (buffer != nullptr) {
            for (int i = 0; i < bufferLength; ++i) {
                buffer[i] = 0.0f;
            }
        }
    }

    // Writes an input sample to the delay line
    void write(float input) {
        if (bufferLength > 0 && buffer != nullptr) {
            writeIndex += 1;
            if (writeIndex >= bufferLength) {
                writeIndex = 0;
            }
            buffer[writeIndex] = input;
        }
    }

    // Reads a sample from the delay line with a specified delay
    float read(int delayInSamples) const {
        if (bufferLength > 0 && buffer != nullptr) {
            if (delayInSamples >= 0 && delayInSamples <= bufferLength - 1) {
                int readIndex = writeIndex - delayInSamples;
                if (readIndex < 0) {
                    readIndex += bufferLength;
                }
                return buffer[readIndex];
            }
        }
        return 0.0f; // Return 0 if delay is out of bounds or buffer is invalid
    }

    // Returns the length of the buffer
    int getBufferLength() const {
        return bufferLength;
    }

    // Destructor
    ~DelayLine() {
        delete[] buffer; // Free allocated memory
    }

private:
    float* buffer = nullptr; // Pointer to the delay buffer
    int bufferLength = 0;    // Length of the buffer
    int writeIndex = 0;      // Index of the most recent value written
};

#endif // DELAYLINE_H

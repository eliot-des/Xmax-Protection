#ifndef MINFILTER_H
#define MINFILTER_H

#include <Arduino.h>

class MinFilter {
public:
    // Default constructor
    MinFilter() : buffer(nullptr), bufferLength(0), currentMin(1.0f) {}

    // Destructor to clean up dynamically allocated memory
    ~MinFilter() {
        delete[] buffer;
    }

    // Resizes the filter buffer and resets the filter
    void resize(int maxSize) {
        if (maxSize <= 0) return;

        if (buffer != nullptr) {
            delete[] buffer; // Free the old buffer
        }

        bufferLength = maxSize;
        buffer = new float[maxSize];
        reset();
    }

    // Sets a new size for the window
    void set(int newSize) {
        if (newSize <= 0 || newSize > bufferLength) return;

        if (newSize < windowSize) {
            // Move tail forward to fit the new window size
            int elementsToSkip = windowSize - newSize;
            tail = (tail + elementsToSkip) % bufferLength;
            currentSize = min(currentSize, newSize);
            recalculateMin();
        }
        windowSize = newSize;
    }

    // Adds a new value and updates the minimum
    void add(float value) {
        if (currentSize < windowSize) {
            currentSize++;
        } else {
            if (buffer[tail] == currentMin && currentMin < max) {
                needsRecalculation = true;
            }
            tail = (tail + 1) % bufferLength;
        }

        buffer[head] = value;
        head = (head + 1) % bufferLength;

        if (value < currentMin) {
            currentMin = value;
        } else if (needsRecalculation) {
            recalculateMin();
        }
    }

    // Returns the current minimum of the window
    float getMinimum() const {
        return currentMin;
    }

    // Resets the filter
    void reset() {
        if (buffer != nullptr) {
            for (int i = 0; i < bufferLength; ++i) {
                buffer[i] = max;
            }
        }
        head = tail = 0;
        currentMin = max;
        currentSize = 0;
        needsRecalculation = false;
    }

    // Sets the maximum possible value for the buffer elements
    void setMax(float newMax) {
        max = newMax;
        recalculateMin();
    }

private:
    // Recalculates the minimum in the window
    void recalculateMin() {
        currentMin = max;
        for (int i = 0; i < currentSize; ++i) {
            int index = (tail + i) % bufferLength;
            if (buffer[index] < currentMin) {
                currentMin = buffer[index];
            }
        }
        needsRecalculation = false;
    }

    float* buffer = nullptr;  // Circular buffer for storing values
    int bufferLength = 0;     // Total buffer length
    int windowSize = 0;       // Current window size
    int head = 0;             // Index of the newest element
    int tail = 0;             // Index of the oldest element
    float currentMin = 1.0f;  // Current minimum value in the window
    float max = 1.0f;         // Maximum value possible for input
    int currentSize = 0;      // Current number of elements in the window
    bool needsRecalculation = false; // Flag for recalculating the minimum
};

#endif // MINFILTER_H



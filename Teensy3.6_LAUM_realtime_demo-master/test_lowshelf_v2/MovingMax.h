#ifndef MOVINGMAX_H
#define MOVINGMAX_H

#include <Arduino.h>

template<class T, int N>
class MovingMax {
public:
    MovingMax() {
        m_scanmax = m_inmax = T(-999999); // Set to a very low value
        for (int i = 0; i < N; i++) {
            m_data[i] = m_inmax;
        }
        m_writepos = 0;
        m_scanend = 0;
        m_scanpos = 0;
        m_scanstart = (N - 1) / 2;
    }

    T add(T x) {
        // Step 1: Update scanning range
        if (--m_scanpos >= m_scanend) {
            m_inmax = max(m_inmax, x);
            m_data[m_scanpos] = max(m_data[m_scanpos], m_data[m_scanpos + 1]);
        } else {
            m_scanmax = m_inmax;
            m_inmax = x;
            m_scanend = m_scanend ^ ((N + 1) / 2);
            m_scanstart = m_scanstart ^ ((N - 1) / 2 ^ (N - 1));
            m_scanpos = m_scanstart;
        }

        // Step 2: Write to the circular buffer
        m_data[m_writepos] = x;
        if (++m_writepos >= N) {
            m_writepos = 0;
        }

        // Step 3: Calculate moving max
        T outmax = m_data[m_writepos];
        return max(m_inmax, max(m_scanmax, outmax));
    }

private:
    int m_writepos;      // Write position in the circular buffer
    int m_scanpos;       // Current scan position
    int m_scanend;       // End of the scan range
    int m_scanstart;     // Start of the scan range
    T m_inmax;           // Maximum value in the current input window
    T m_scanmax;         // Maximum value in the scanning range
    T m_data[N];         // Circular buffer to store input values
};

#endif // MOVINGMAX_H



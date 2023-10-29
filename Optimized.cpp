// SIMD

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <ctime>
#include <emmintrin.h>

const int N = 5000; // Size of the 2D square array

// Define a complex number
using Complex = std::complex<double>;

// Cooley-Tukey FFT algorithm
void fft(std::vector<Complex>& data) {
    const int n = data.size();
    if (n <= 1) return;

    std::vector<Complex> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < n / 2; k++) {
        Complex t = std::polar(1.0, -2.0 * M_PI * k / n) * odd[k];
        data[k] = even[k] + t;
        data[k + n / 2] = even[k] - t;
    }
}

// Perform FFT using SSE intrinsics
void fft_sse(std::vector<Complex>& data) {
    const int n = data.size();
    if (n <= 1) return;

    std::vector<Complex> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    fft_sse(even);
    fft_sse(odd);

    // Use SSE intrinsics for vectorized complex number operations
    __m128d w = _mm_setr_pd(1.0, -2.0 * M_PI / n);
    for (int k = 0; k < n / 2; k += 2) {
        __m128d odd_vec = _mm_loadu_pd(reinterpret_cast<double*>(&odd[k]));
        __m128d even_vec = _mm_loadu_pd(reinterpret_cast<double*>(&even[k]));

        __m128d t_real = _mm_mul_pd(_mm_shuffle_pd(w, w, 1), odd_vec);
        __m128d t_imag = _mm_mul_pd(w, odd_vec);

        __m128d new_data_real = _mm_add_pd(even_vec, t_real);
        __m128d new_data_imag = _mm_sub_pd(even_vec, t_real);

        _mm_storeu_pd(reinterpret_cast<double*>(&data[k]), new_data_real);
        _mm_storeu_pd(reinterpret_cast<double*>(&data[k + n / 2]), new_data_imag);

        w = _mm_add_pd(w, _mm_mul_pd(_mm_set1_pd(2.0), _mm_mul_pd(w, _mm_setr_pd(-0.5, 0.5))));
    }
}

int main() {
    // Create a large 2D array (N x N)
    std::vector<std::vector<Complex>> data(N, std::vector<Complex>(N));

    // Initialize data with some values (e.g., a pattern)

    // Start the timer
    std::time_t start_time = std::time(nullptr);

    // Perform 2D FFT on the data using SSE
    for (int i = 0; i < N; i++) {
        fft_sse(data[i]);
    }
    for (int j = 0; j < N; j++) {
        std::vector<Complex> column(N);
        for (int i = 0; i < N; i++) {
            column[i] = data[i][j];
        }
        fft_sse(column);
        for (int i = 0; i < N; i++) {
            data[i][j] = column[i];
        }
    }

    // End the timer
    std::time_t end_time = std::time(nullptr);

    std::cout << "2D FFT completed in " << difftime(end_time, start_time) << " seconds." << std::endl;

    return 0;
}
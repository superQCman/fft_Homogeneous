#include <iostream>
#include <complex>
#include <vector>
#include <thread>
#include <cmath>
#include "apis_c.h"

const double PI = 3.141592653589793238460;

using Complex = std::complex<double>;
using CArray = std::vector<Complex>;

// FFT 算法
void fft(CArray& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    // 分解为偶数和奇数部分
    CArray even(N / 2);
    CArray odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // 递归调用 FFT
    fft(even);
    fft(odd);

    // 组合结果
    for (size_t k = 0; k < N / 2; ++k) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

// 多线程版本的 FFT 调用
void parallel_fft(CArray& x, size_t num_threads, int idX, int idY) {
    const size_t N = x.size();
    if (N <= 1) return;

    if (num_threads <= 1 || N <= 16) {
        // 如果线程数少于等于 1，或者数据量很小，则使用单线程 FFT
        fft(x);
        return;
    }

    // 分解为偶数和奇数部分
    CArray even(N / 2);
    CArray odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    double *even_first = new double[even.size()];
    double *even_second = new double[even.size()];
    double *odd_first = new double[odd.size()];
    double *odd_second = new double[odd.size()];
    for (size_t i = 0; i < N / 2; ++i) {
        even_first[i] = even[i].real();
        even_second[i] = even[i].imag();
        odd_first[i] = odd[i].real();
        odd_second[i] = odd[i].imag();
    }
    int *num_threads_even = new int;
    int *num_threads_odd = new int;
    *num_threads_even = num_threads / 2;
    *num_threads_odd = num_threads / 2;
    double* even_receive_first = new double[even.size()];
    double* even_receive_second = new double[even.size()];
    double* odd_receive_first = new double[odd.size()];
    double* odd_receive_second = new double[odd.size()];

    // 递归调用 FFT，使用多线程
    InterChiplet::sendMessage(0, 1, idX, idY, even_first, even.size() * sizeof(double));
    InterChiplet::sendMessage(0, 1, idX, idY, even_second, even.size() * sizeof(double));
    InterChiplet::sendMessage(0, 2, idX, idY, odd_first, odd.size() * sizeof(double));
    InterChiplet::sendMessage(0, 2, idX, idY, odd_second, odd.size() * sizeof(double));

    InterChiplet::sendMessage(0, 1, idX, idY, num_threads_even, sizeof(int));
    InterChiplet::sendMessage(0, 2, idX, idY, num_threads_odd, sizeof(int));

    InterChiplet::receiveMessage(idX, idY, 0, 1, even_receive_first, even.size() * sizeof(double));
    InterChiplet::receiveMessage(idX, idY, 0, 1, even_receive_second, even.size() * sizeof(double));
    InterChiplet::receiveMessage(idX, idY, 0, 2, odd_receive_first, odd.size() * sizeof(double));
    InterChiplet::receiveMessage(idX, idY, 0, 2, odd_receive_second, odd.size() * sizeof(double));

    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = Complex(even_receive_first[i], even_receive_second[i]);
        odd[i] = Complex(odd_receive_first[i], odd_receive_second[i]);
    }

    // std::thread even_thread(parallel_fft, std::ref(even), num_threads / 2);
    // std::thread odd_thread(parallel_fft, std::ref(odd), num_threads / 2);

    // // 等待子线程完成
    // even_thread.join();
    // odd_thread.join();

    // 组合结果
    for (size_t k = 0; k < N / 2; ++k) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }

    delete[] even_first;
    delete[] even_second;
    delete[] odd_first;
    delete[] odd_second;
    delete num_threads_even;
}

int main(int argc, char **argv) {
    int idX = atoi(argv[1]);
    int idY = atoi(argv[2]);
    // 输入数据
    const size_t N = std::pow(2,30);  // 输入长度必须是 2 的幂
    CArray data(N);

    // 初始化数据
    for (size_t i = 0; i < N; ++i) {
        data[i] = Complex(i, 0);  // 实部为 i，虚部为 0
    }

    // 打印输入
    std::cout << "Input data:" << std::endl;
    for (size_t i = 0; i < 8; ++i) {  // 仅打印前8个元素
        std::cout << data[i] << std::endl;
    }

    // 使用 4 个线程执行 FFT
    size_t num_threads = 4;
    parallel_fft(data, num_threads, idX, idY);

    // 打印输出
    std::cout << "\nFFT result:" << std::endl;
    for (size_t i = 0; i < 8; ++i) {  // 仅打印前8个元素
        std::cout << data[i] << std::endl;
    }

    return 0;
}

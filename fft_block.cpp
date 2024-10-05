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
void parallel_fft(CArray& x, size_t num_threads,int idX, int idY) {
    const size_t N = x.size();
    
    if (N <= 1) return;

    if (num_threads <= 1 || N <= 16) {
        // 如果线程数少于等于 1，或者数据量很小，则使用单线程 FFT
        fft(x);
        std::cout << "##################### N_1: " << N << std::endl;
        // 发送结果回 fft_main.cpp
        double* result_real = new double[N];
        double* result_imag = new double[N];
        for (size_t i = 0; i < N; ++i) {
            result_real[i] = x[i].real();
            result_imag[i] = x[i].imag();
        }
        if (idY == 1 || idY == 2){
            InterChiplet::sendMessage(0, 1, idX, idY, result_real, N * sizeof(double));
            InterChiplet::sendMessage(0, 1, idX, idY, result_imag, N * sizeof(double));
        }
        else if(idY == 4 || idY == 3){
            InterChiplet::sendMessage(0, 2, idX, idY, result_real, N * sizeof(double));
            InterChiplet::sendMessage(0, 2, idX, idY, result_imag, N * sizeof(double));
        }
        
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
    if (idX == 0 && idY == 1){
        InterChiplet::sendMessage(1, 1, idX, idY, even_first, even.size() * sizeof(double));
        InterChiplet::sendMessage(1, 1, idX, idY, even_second, even.size() * sizeof(double));
        InterChiplet::sendMessage(1, 2, idX, idY, odd_first, odd.size() * sizeof(double));
        InterChiplet::sendMessage(1, 2, idX, idY, odd_second, odd.size() * sizeof(double));

        InterChiplet::sendMessage(1, 1, idX, idY, num_threads_even, sizeof(int));
        InterChiplet::sendMessage(1, 2, idX, idY, num_threads_odd, sizeof(int));

        InterChiplet::receiveMessage(idX, idY, 1, 1, even_receive_first, even.size() * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 1, 1, even_receive_second, even.size() * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 1, 2, odd_receive_first, odd.size() * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 1, 2, odd_receive_second, odd.size() * sizeof(double));
    }
    else if(idX == 0 && idY == 2){
        InterChiplet::sendMessage(1, 3, idX, idY, even_first, even.size() * sizeof(double));
        InterChiplet::sendMessage(1, 3, idX, idY, even_second, even.size() * sizeof(double));
        InterChiplet::sendMessage(1, 4, idX, idY, odd_first, odd.size() * sizeof(double));
        InterChiplet::sendMessage(1, 4, idX, idY, odd_second, odd.size() * sizeof(double));

        InterChiplet::sendMessage(1, 3, idX, idY, num_threads_even, sizeof(int));
        InterChiplet::sendMessage(1, 4, idX, idY, num_threads_odd, sizeof(int));

        InterChiplet::receiveMessage(idX, idY, 1, 3, even_receive_first, even.size() * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 1, 3, even_receive_second, even.size() * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 1, 4, odd_receive_first, odd.size() * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 1, 4, odd_receive_second, odd.size() * sizeof(double));
    }
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = Complex(even_receive_first[i], even_receive_second[i]);
        odd[i] = Complex(odd_receive_first[i], odd_receive_second[i]);
    }

    // 组合结果
    for (size_t k = 0; k < N / 2; ++k) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
    // 发送结果回 fft_main.cpp
    double* result_real = new double[N];
    double* result_imag = new double[N];
    for (size_t i = 0; i < N; ++i) {
        result_real[i] = x[i].real();
        result_imag[i] = x[i].imag();
    }

    InterChiplet::sendMessage(0, 0, idX, idY, result_real, N * sizeof(double));
    InterChiplet::sendMessage(0, 0, idX, idY, result_imag, N * sizeof(double));

    delete[] result_real;
    delete[] result_imag;

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
    const size_t N = std::pow(2,30)/std::pow(2,idX+1);  // 数据量
    CArray data(N);

    // 接收从 fft_main.cpp 发送过来的数据
    double* real_part = new double[N];
    double* imag_part = new double[N];
    if (idX == 0){
        InterChiplet::receiveMessage(idX, idY, 0, 0, real_part, N * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 0, 0, imag_part, N * sizeof(double));
    }else if(idY == 1 || idY == 2){
        InterChiplet::receiveMessage(idX, idY, 0, 1, real_part, N * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 0, 1, imag_part, N * sizeof(double));
    }else if(idY == 3 || idY == 4){
        InterChiplet::receiveMessage(idX, idY, 0, 2, real_part, N * sizeof(double));
        InterChiplet::receiveMessage(idX, idY, 0, 2, imag_part, N * sizeof(double));
    }

    for (size_t i = 0; i < N; ++i) {
        data[i] = Complex(real_part[i], imag_part[i]);
    }

    delete[] real_part;
    delete[] imag_part;

    // 打印输入
    std::cout << "------------------------------------------Input data:------------------------------------------" << std::endl;
    for (size_t i = 0; i < 8; ++i) {  // 仅打印前8个元素
        std::cout << data[i] << std::endl;
    }

    // 接受从 fft_main.cpp 发送过来的线程数
    int *num_threads = new int;
    if (idX == 0){
        InterChiplet::receiveMessage(idX, idY, 0, 0, num_threads, sizeof(int));
    }else if (idY == 1 || idY == 2){
        InterChiplet::receiveMessage(idX, idY, 0, 1, num_threads, sizeof(int));
    }else if (idY == 3 || idY == 4){
        InterChiplet::receiveMessage(idX, idY, 0, 2, num_threads, sizeof(int));
    }
    parallel_fft(data, *num_threads,idX,idY);

    // 打印输出
    std::cout << "\n------------------------------------------FFT result:------------------------------------------" << std::endl;
    for (size_t i = 0; i < 8; ++i) {  // 仅打印前8个元素
        std::cout << data[i] << std::endl;
    }

    return 0;
}

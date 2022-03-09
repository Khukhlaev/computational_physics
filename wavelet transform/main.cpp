#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>


std::vector<double> waveletHaarTransform(std::vector<double> source) {

    if (source.size() == 1) return source;

    std::vector<double> result, tmp;
    for (int i = 0; i < source.size() - 1; i += 2) {
        result.push_back((source[i] - source[i + 1]) / 2.0);
        tmp.push_back((source[i] - source[i + 1]) / 2.0);
    }

    tmp = waveletHaarTransform(tmp);
    result.insert(result.end(), tmp.begin(), tmp.end());

    return result;
}


void DaubechiesTransformIteration(std::vector<double>& signal, int upper_bound,
                                  const std::vector<double>& h, const std::vector<double>& g) {

    int i, j;
    int n_transform = h.size();

    if (upper_bound < n_transform) return;

    int half = int(upper_bound / 2);

    std::vector<double> tmp(upper_bound);

    i = 0;
    for (j = 0; j < upper_bound; j += 2) {
        for (int k = 0; k < n_transform; k++) {
            tmp[i] += signal[(j + k) % upper_bound] * h[k];
            tmp[i + half] += signal[(j + k) % upper_bound] * g[k];
        }
        i++;
    }

    for (int k = 0; k < upper_bound; k++) {
        signal[k] = tmp[k];
    }

}

std::vector<double> waveletDaubechiesTransform(const std::vector<double>& signal, const std::vector<double>& h,
                                               const std::vector<double>& g, int level = -1) {

    std::vector<double> result = signal;

    if (signal.size() < h.size()) {
        std::cout << "ERROR, signal size must be bigger than coefficient size" << std::endl;
        return result;
    }

    for (int i = signal.size(); i >= h.size(); i /= 2) {
        level--;
        DaubechiesTransformIteration(result, i, h, g);

        if (level == 0) break;
    }

    return result;
}


void invDaubechiesTransformIteration(std::vector<double>& signal, int upper_bound,
                                  const std::vector<double>& h, const std::vector<double>& g) {

    int i, j;
    int n_transform = h.size();

    if (upper_bound < n_transform) return;

    int half = int(upper_bound / 2);

    std::vector<double> tmp(upper_bound);

    for (int k = 0; k < n_transform / 2; k++) {
        tmp[0] += signal[(half - n_transform / 4 + k) % half] * h[n_transform - 2 - 2 * k];
        tmp[0] += signal[(half - n_transform / 4 + k) % half + half] * g[n_transform - 2 - 2 * k];
        tmp[1] += signal[(half - n_transform / 4 + k) % half] * h[n_transform - 1 - 2 * k];
        tmp[1] += signal[(half - n_transform / 4 + k) % half + half] * g[n_transform - 1 - 2 * k];

    }

    j = 2;
    for (i = 0; i < half - 1; i++) {
        for (int k = 0; k < n_transform / 2; k++) {
            tmp[j] += signal[(i + k) % upper_bound] * h[n_transform - 2 - 2 * k];
            tmp[j] += signal[(i + k + half) % upper_bound] * g[n_transform - 2 - 2 * k];
            tmp[j + 1] += signal[(i + k) % upper_bound] * h[n_transform - 1 - 2 * k];
            tmp[j + 1] += signal[(i + k + half) % upper_bound] * g[n_transform - 1 - 2 * k];
        }
        j += 2;
    }

    for (int k = 0; k < upper_bound; k++) {
        signal[k] = tmp[k];
    }

}

std::vector<double> invWaveletDaubechiesTransform(const std::vector<double>& signal, const std::vector<double>& h,
                                               const std::vector<double>& g, int level = -1) {

    std::vector<double> result = signal;

    if (signal.size() < h.size()) {
        std::cout << "ERROR, signal size must be bigger than coefficient size" << std::endl;
        return result;
    }
    int start_i = h.size();

    if (level == 1) start_i = signal.size();

    for (int i = start_i; i <= signal.size(); i *= 2) {
        invDaubechiesTransformIteration(result, i, h, g);
    }

    return result;
}




void fillCoeffs (std::vector<std::pair<std::vector<double>, std::vector<double>>>& coeffs) {
    std::ifstream in("dw_coeff.txt");

    std::string buffer;
    int current_n;
    double coef;

    while (std::getline(in, buffer)) {
        if (buffer.length() == 0) continue;

        std::stringstream stream (buffer);
        buffer.clear();
        stream >> buffer;
        if (buffer[0] == 'N')  {
            stream >> current_n;
            continue;
        }

        stream >> coef;

        coeffs[current_n - 2].first.push_back(coef);

    }

    double sign = 1;
    for (int i = 1; i < coeffs.size(); i++) {
        coeffs[i].second.resize(coeffs[i].first.size());

        for (int j = 0; j < coeffs[i].first.size(); j++) {
            coeffs[i].second[j] = sign * coeffs[i].first[coeffs[i].first.size() - j - 1];

            sign *= -1;
        }
    }

}




int main() {

    // Daubechies Transform

    std::vector<std::pair<std::vector<double>, std::vector<double>>> coeffs(9);

    coeffs[0].first.push_back(0.4829629131445341433748715998644486838169524195042022752011715382);
    coeffs[0].first.push_back(0.836516303737807905575293780916873203459370388348439293495341473);
    coeffs[0].first.push_back(0.224143868042013381025972762240400355467883518184271761387168331);
    coeffs[0].first.push_back(-0.1294095225512603811744494188120241641745344506599652569070016037);

    coeffs[0].second.resize(4);
    coeffs[0].second[0] = coeffs[0].first[3];
    coeffs[0].second[1] = -coeffs[0].first[2];
    coeffs[0].second[2] = coeffs[0].first[1];
    coeffs[0].second[3] = -coeffs[0].first[0];

    fillCoeffs(coeffs);

    std::vector<double> signal = {0,0,0,0,1,0,0,0};
    std::vector<double> res = waveletDaubechiesTransform(signal, coeffs[0].first, coeffs[0].second);

    std::vector<double> res1 = invWaveletDaubechiesTransform(res, coeffs[0].first, coeffs[0].second);

    for (auto x : res1) {
        std::cout << x << " ";
    }

    std::cout << std::endl;

//    std::ofstream out("out.txt");
//
//    std::vector<double> signal1(4096, 0);
//    for (int i = 0; i < signal1.size(); i++) {
//        signal1[i] = sin(250 * 3.14 * double(i) / 4096 * double (i) / 4096) ;
//    }
//
//
//    std::vector<double> res3 = waveletDaubechiesTransform(signal1, coeffs[0].first, coeffs[0].second, 1);
//
//    for (auto x : res3) {
//        out << x << " ";
//    }

    std::ofstream out1("out1.txt");

    std::vector<double> signal2(4096, 0);
    signal2[2048] = 1;

    std::vector<double> res4 = signal2;
    for (int i = 0; i < 5; i++) {
        res4 = invWaveletDaubechiesTransform(res4, coeffs[2].first, coeffs[2].second, 1);
    }

    for (auto x : res4) {
        out1 << x << " ";
    }

    return 0;
}

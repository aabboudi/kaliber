#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "kaliber.h"

namespace Kaliber {

/***********************************************/
/***********************************************/
/****************** UTILITIES ******************/
/***********************************************/
/***********************************************/

const double PI = 3.14159265358979323846;

double erf_inv(double x, double epsilon = 1e-8, int max_it = 100) {
    double guess = 0.0;

    for (int i = 0; i < max_it; ++i) {
        double error = std::erf(guess) - x;
        if(std::fabs(error) < epsilon) { return guess; }

        // Update guess using Newton-Raphson method
        guess -= error / (std::sqrt(M_PI) * std::exp(-guess * guess));
    }
    return guess;
}

unsigned long long factorial(int k) {
    return ( k==0 || k==1 ) ? 1 : k*factorial(k-1);
}

/***********************************************/
/***********************************************/
/*********** DISTRIBUTION PRINCIPLES ***********/
/***********************************************/
/***********************************************/

double Gaussian::probability(double x) const { return 1.0 / std::sqrt(2.0 * PI) * std::exp(-0.5 * x * x); }

double Gaussian::cumulativeProbability(double x) const { return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0))); }

Binomial::Binomial(int n, double p) : n(n), p(p) {}

double Binomial::probability(double k) const {
    return std::exp(std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1)) * std::pow(p, k) * std::pow(1.0 - p, n - k);
}

double Binomial::cumulativeProbability(double k) const {
    double result = 0.0;
    for (int i = 0; i <= k; ++i) { result += probability(i); }
    return result;
}

Bernoulli::Bernoulli(double p) : Binomial(1,p) {}

double Bernoulli::probability(double x) const {
    if(x == 0 || x == 1) { return std::pow(1 - p, 1 - x) * std::pow(p, x); }
    else { return 0.0; }
}

double Bernoulli::cumulativeProbability(double x) const {
    if(x < 0)          { return 0.0; }
    else if(x < 1)     { return 1 - p; }
    else if(x >= 1)    { return 1.0; }
    else               { return 0.0; }
}

Poisson::Poisson(double lambda) : lambda(lambda) {}

double Poisson::probability(double k) const {
    return std::exp(-lambda) * std::pow(lambda, k) / factorial(k);
}

double Poisson::cumulativeProbability(double k) const {
    double result = 0.0;
    for (int i = 0; i <= k; ++i) { result += probability(i); }
    return result;
}

Exponential::Exponential(double lambda) : lambda(lambda) {}

double Exponential::probability(double x) const {
    return lambda * std::exp(-lambda * x);
}

double Exponential::cumulativeProbability(double x) const {
    return 1.0 - std::exp(-lambda * x);
}

Geometric::Geometric(double p) : p(p) {}

double Geometric::probability(double k) const {
    if(k < 0 || p <= 0.0 || p >= 1.0) {
        return 0.0;
    }

    return std::pow(1.0 - p, k) * p;
}

double Geometric::cumulativeProbability(double k) const {
    if(k < 0 || p <= 0.0 || p >= 1.0) { return 0.0; }
    return 1.0 - std::pow(1.0 - p, k + 1);
}

Uniform::Uniform(double a, double b) : a(a), b(b) {}

double Uniform::probability(double x) const {
    if(x < a || x > b || a >= b) { return 0.0; }
    return 1.0 / (b - a);
}

double Uniform::cumulativeProbability(double x) const {
    if(x < a) { return 0.0; }
    else if(x >= b) { return 1.0; }
    return (x - a) / (b - a);
}

LogNormal::LogNormal(double mean, double stddev) : mean(mean), stddev(stddev) {}

double LogNormal::probability(double x) const {
  
    if(x <= 0.0 || stddev <= 0.0) { return 0.0; }

    double exponent = -0.5 * std::pow((std::log(x) - mean) / stddev, 2);
    return (1.0 / (x * stddev * std::sqrt(2.0 * PI))) * std::exp(exponent);
}

double LogNormal::cumulativeProbability(double x) const {

    if(x <= 0.0 || stddev <= 0.0) { return 0.0; }
    double z = (std::log(x) - mean) / stddev;
    return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
}

/***********************************************/
/***********************************************/
/*********** REGRESSIONAL COMPUTATION **********/
/***********************************************/
/***********************************************/

Regression::Regression(const std::vector<double>& x, const std::vector<double>& y)
    : x_data(x), y_data(y) {}

// Implementation of LinearRegression class
LinearRegression::LinearRegression(const std::vector<double>& x, const std::vector<double>& y)
    : Regression(x, y), slope(0), intercept(0) {}

void LinearRegression::fit() {
    if(x_data.size() != y_data.size() || x_data.empty()) {
        std::cerr << "Error: Invalid input data for linear regression fitting." << std::endl;
        return;
    }

    double mean_x = 0, mean_y = 0;
    for (size_t i = 0; i < x_data.size(); ++i) {
        mean_x += x_data[i] / x_data.size();
        mean_y += y_data[i] / y_data.size();
    }

    double numerator = 0, denominator = 0;
    for (size_t i = 0; i < x_data.size(); ++i) {
        numerator += (x_data[i] - mean_x) * (y_data[i] - mean_y);
        denominator += (x_data[i] - mean_x) * (x_data[i] - mean_x);
    }

    if(denominator != 0) {
        slope = numerator / denominator;
        intercept = mean_y - slope * mean_x;
    } else {
        std::cerr << "Error: Denominator is zero. Unable to fit the linear regression model." << std::endl;
    }
}

double LinearRegression::predict(double input) const { return slope * input + intercept; }
double LinearRegression::getSlope() const { return slope; }
double LinearRegression::getIntercept() const { return intercept; }

// Cubic
CubicRegression::CubicRegression(const std::vector<double>& x, const std::vector<double>& y)
    : Regression(x, y), coefficients(std::vector<double>(4, 0.0)) {}

void CubicRegression::fit() {
    if(x_data.size() != y_data.size() || x_data.empty()) {
        std::cerr << "Error: Invalid input data for cubic regression fitting." << std::endl;
        return;
    }

    int degree = coefficients.size() - 1;

    // Create augmented matrix [X | Y]
    std::vector<std::vector<double>> augmentedMatrix(x_data.size(), std::vector<double>(degree + 2, 0.0));

    for (size_t i = 0; i < x_data.size(); ++i) {
        for (int j = 0; j <= degree; ++j) {
            augmentedMatrix[i][j] = pow(x_data[i], degree - j);
        }
        augmentedMatrix[i][degree + 1] = y_data[i];
    }

    // Gaussian elimination
    for (int i = 0; i < degree; ++i) {
        for (int j = i + 1; j < degree + 1; ++j) {
            double ratio = augmentedMatrix[j][i] / augmentedMatrix[i][i];
            for (int k = 0; k < degree + 2; ++k) {
                augmentedMatrix[j][k] -= ratio * augmentedMatrix[i][k];
            }
        }
    }

    // Back substitution
    for (int i = degree; i >= 0; --i) {
        coefficients[i] = augmentedMatrix[i][degree + 1];
        for (int j = i + 1; j < degree + 1; ++j) {
            coefficients[i] -= augmentedMatrix[i][j] * coefficients[j];
        }
        coefficients[i] /= augmentedMatrix[i][i];
    }
}

double CubicRegression::predict(double input) const {
    double result = 0.0;
    int degree = coefficients.size() - 1;
    for (int i = 0; i <= degree; ++i) { result += coefficients[i] * pow(input, degree - i); }
    return result;
}

const std::vector<double>& CubicRegression::getCoefficients() const { return coefficients; }

// Polynomial
PolynomialRegression::PolynomialRegression(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& coeffs)
    : Regression(x, y), coefficients(coeffs) {}

void PolynomialRegression::fit() {
    if(x_data.size() != y_data.size() || x_data.empty() || coefficients.size() == 0) {
        std::cerr << "Error: Invalid input data or coefficients for polynomial regression fitting." << std::endl;
        return;
    }

    int degree = coefficients.size() - 1;

    // Create augmented matrix [X | Y]
    std::vector<std::vector<double>> augmentedMatrix(x_data.size(), std::vector<double>(degree + 2, 0.0));

    for (size_t i = 0; i < x_data.size(); ++i) {
        for (int j = 0; j <= degree; ++j) {
            augmentedMatrix[i][j] = pow(x_data[i], degree - j);
        }
        augmentedMatrix[i][degree + 1] = y_data[i];
    }

    // Gaussian elimination
    for (int i = 0; i < degree; ++i) {
        for (int j = i + 1; j < degree + 1; ++j) {
            double ratio = augmentedMatrix[j][i] / augmentedMatrix[i][i];
            for (int k = 0; k < degree + 2; ++k) {
                augmentedMatrix[j][k] -= ratio * augmentedMatrix[i][k];
            }
        }
    }

    // Back substitution
    for (int i = degree; i >= 0; --i) {
        coefficients[i] = augmentedMatrix[i][degree + 1];
        for (int j = i + 1; j < degree + 1; ++j) {
            coefficients[i] -= augmentedMatrix[i][j] * coefficients[j];
        }
        coefficients[i] /= augmentedMatrix[i][i];
    }
}

double PolynomialRegression::predict(double input) const {
    double result = 0.0;
    int degree = coefficients.size() - 1;
    for (int i = 0; i <= degree; ++i) { result += coefficients[i] * pow(input, degree - i); }
    return result;
}

const std::vector<double>& PolynomialRegression::getCoefficients() const { return coefficients; }
}

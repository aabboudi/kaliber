#include "kaliber.h"
#include <iostream>
#include <vector>

using namespace Kaliber;

int main() {

    // DISTRIBUTION
    // Gaussian
    Gaussian gaussianDist;
    double x1 = 1.5;
    std::cout << "Gaussian Probability at x=" << x1 << " is: " << gaussianDist.probability(x1) << std::endl;
    std::cout << "Gaussian Cumulative Probability at x=" << x1 << " is: " << gaussianDist.cumulativeProbability(x1) << std::endl;

    // Binomial Distribution
    Binomial binomialDist(5, 0.5);
    double k2 = 2;
    std::cout << "Binomial Probability at k=" << k2 << ": " << binomialDist.probability(k2) << std::endl;
    std::cout << "Binomial Cumulative Probability at k=" << k2 << " is: " << binomialDist.cumulativeProbability(k2) << std::endl;

    // Poisson Distribution
    Poisson poissonDist(3.0);
    double k3 = 2;
    std::cout << "Poisson Probability at k=" << k3 << ": " << poissonDist.probability(k3) << std::endl;
    std::cout << "Poisson Cumulative Probability at k=" << k3 << " is: " << poissonDist.cumulativeProbability(k3) << std::endl;



    // REGRESSION
    std::vector<double> x_data = {1, 2, 3, 4, 5, 6};
    std::vector<double> y_data = {2, 4, 6, 8, 10};
    double input = 6;


    // Linear Regression
    LinearRegression linearRegression(x_data, y_data);
    linearRegression.fit();
    std::cout << "Linear Regression Prediction at input=" << input << " is: " << linearRegression.predict(input) << std::endl;

    // Cubic Regression
    CubicRegression cubicRegression(x_data, y_data);
    cubicRegression.fit();
    std::cout << "Cubic Regression Prediction at input=" << input << " is: " << cubicRegression.predict(input) << std::endl;

    // Polynomial Regression
    std::vector<double> polynomialCoeffs = {1, -1, 2};
    PolynomialRegression polyRegression(x_data, y_data, polynomialCoeffs);
    polyRegression.fit();
    std::cout << "Polynomial Regression Prediction at input=" << input << " is: " << polyRegression.predict(input) << std::endl;

    return 0;
}

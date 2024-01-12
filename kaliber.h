#ifndef KALIBER_H
#define KALIBER_H

#include <cmath>
#include <vector>
#include <iostream>

namespace Kaliber {

double erf_inv(double, double, int);
unsigned long long factorial(int k);

class Dist {
public:
    virtual double probability(double x) const = 0;
    virtual double cumulativeProbability(double x) const = 0;
    virtual ~Dist() = default;
};

// Gaussian Dist
class Gaussian : public Dist {
public:
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;
};

// Binomial Dist
class Binomial : virtual public Dist {
protected:
    int n;
    double p;

public:
    Binomial(int n, double p);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;
};

class Bernoulli : virtual public Binomial {

public:
    Bernoulli(double p);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;
};

// Poisson Dist
class Poisson : public Dist {
public:
    Poisson(double lambda);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    double lambda;
};

// Exponential Dist
class Exponential : public Dist {
public:
    Exponential(double lambda);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    double lambda;
};

// Chi-squared Dist
class ChiSquared : public Dist {
public:
    ChiSquared(int degreesOfFreedom);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    int degreesOfFreedom;
};

// Normal Dist Quantiles
class NormalQuantile : public Dist {
public:
    NormalQuantile(double mean, double stddev);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    double mean;
    double stddev;
};

// Geometric Dist
class Geometric : public Dist {
public:
    Geometric(double p);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    double p;
};

// Uniform Dist
class Uniform : public Dist {
public:
    Uniform(double a, double b);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    double a;
    double b;
};

// Log-Normal Dist
class LogNormal : public Dist {
public:
    LogNormal(double mean, double stddev);
    double probability(double x) const override;
    double cumulativeProbability(double x) const override;

private:
    double mean;
    double stddev;
};


// REGRESSION
class Regression {
protected:
    std::vector<double> x_data;
    std::vector<double> y_data;

public:
    Regression(const std::vector<double>& x, const std::vector<double>& y);
    virtual void fit() = 0;
    virtual double predict(double input) const = 0;
    virtual ~Regression() = default;
};

// Linear Regression
class LinearRegression : public Regression {
private:
    double slope;
    double intercept;

public:
    LinearRegression(const std::vector<double>& x, const std::vector<double>& y);
    void fit() override;
    double predict(double input) const override;
    double getSlope() const;
    double getIntercept() const;
};

// Cubic Regression
class CubicRegression : public Regression {
private:
    std::vector<double> coefficients;

public:
    CubicRegression(const std::vector<double>& x, const std::vector<double>& y);
    void fit() override;
    double predict(double input) const override;
    const std::vector<double>& getCoefficients() const;
};

// Polynomial Regression
class PolynomialRegression : public Regression {
private:
    std::vector<double> coefficients;

public:
    PolynomialRegression(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& coeffs);
    void fit() override;
    double predict(double input) const override;
    const std::vector<double>& getCoefficients() const;
};

class HypothesisTest {
private:

};
}
#endif // KALIBER_H

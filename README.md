# Kaliber C++ Library

Kaliber is a C++ library that provides tools for probability computation and regression analysis. It includes various probability distribution functions and regression models, allowing users to perform statistical calculations and fit data with regression models.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
  - [Manual Installation](#manual-installation)
  - [Vcpkg Installation](#vcpkg-installation)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Features

### Probability Distributions
- Gaussian Distribution
- Binomial Distribution
- Bernoulli Distribution
- Poisson Distribution
- Exponential Distribution
- Chi-squared Distribution
- Normal Distribution Quantiles
- Geometric Distribution
- Uniform Distribution
- Log-Normal Distribution

### Regression Models
- Linear Regression
- Cubic Regression
- Polynomial Regression

### Utilities
- Error function inverse (`erf_inv`)
- Factorial computation (`factorial`)

## Installation

### Manual Installation

To use Kaliber in your C++ project manually, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/aabboudi/kaliber.git
   ```
2. Include the kaliber.h header file in your project and compile your project with kaliber.cpp:
    ```bash
    g++ -o your_program your_program.cpp kaliber.cpp
    ```

### Vcpkg Installation

Kaliber supports Vcpkg for easy integration into your project. Follow these steps:

1. Install Vcpkg by following the instructions [here](https://example.com).

2. Open a terminal and navigate to your project directory.

3. Install Kaliber using Vcpkg:

    ```bash
    vcpkg install kaliber
    ```

4. Integrate Kaliber into your CMake project:

    ```cmake
    find_package(kaliber CONFIG REQUIRED)
    target_link_libraries(your_target PRIVATE kaliber)
    ```

## Usage
To use Kaliber in your C++ code, include the `kaliber.h` header file in your project and make use of the provided classes and functions.

    ```cpp
    #include "kaliber.h"

    int main() {
        // Example usage of Gaussian Distribution
        Kaliber::Gaussian gaussian;
        double x = 1.5;
        double probability = gaussian.probability(x);
        double cumulativeProbability = gaussian.cumulativeProbability(x);

        // Example usage of Linear Regression
        std::vector<double> x_data = {1.0, 2.0, 3.0};
        std::vector<double> y_data = {2.0, 4.0, 5.0};
        Kaliber::LinearRegression linearRegression(x_data, y_data);
        linearRegression.fit();
        double prediction = linearRegression.predict(4.0);

        return 0;
    }

    ```

## Contributing

If you would like to contribute to Kaliber, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Make your changes and commit them with clear messages.
4. Create a pull request.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for further details.
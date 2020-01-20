#include <iostream>
#include <random>
#include <cmath>
#include <chrono>

using namespace std;

template<class T>
struct distribution {
    mt19937_64 engine;
    T distr;

    distribution(double a, double b, int seed) : engine(seed), distr(a, b) {}

    double operator()() {
        return distr(engine);
    }
};

struct old_implementation {
    distribution<uniform_real_distribution<>> uniform;
    distribution<normal_distribution<>> normal;
    double shape;

    old_implementation(double _shape, int uniform_seed, int normal_seed) 
        : uniform(0.0, 1.0, uniform_seed), normal(0.0, 1.0, normal_seed), shape(_shape) {}

    double operator()() {
        double b, c;
        double U, V, X, Y;

        b = shape - 1. / 3.;
        c = 1. / sqrt(9 * b);
        for (;;) {
            do {
                X = normal();
                V = 1.0 + c * X;
            } while (V <= 0.0);

            V = V * V * V;
            U = uniform();
            if (U < 1.0 - 0.0331 * (X * X) * (X * X)) {
                return b * V;
            }
            if (log(U) < 0.5 * X * X + b * (1. - V + log(V))) {
                return b * V;
            }
        }
    }
};

struct new_implementation {
    distribution<uniform_real_distribution<>> uniform;
    double alpha;

    new_implementation(double _shape, int uniform_seed) 
        : uniform(0.0, 1.0, uniform_seed), alpha(_shape) {}
    
    double operator()() {
        double alpha_p, beta_p, K_p, x_p;
        double prod, h;
        
        alpha_p = floor(alpha);
        if (alpha < 2.0) {
            beta_p = 1.0 / alpha;
            K_p = exp(1.0 - alpha) * pow(alpha, alpha - 1.0);
        }
        else {
            beta_p = (alpha_p - 1.0) / (alpha - 1.0);
            K_p = exp(alpha_p - alpha) * pow(alpha - 1.0, alpha - alpha_p);
        }
        for (;;) {
            prod = 1.0;
            for (int i = 0; i < alpha_p; i++) {
                prod *= uniform();
            }
            x_p = -1.0 / beta_p * log(prod);
            h = (pow(x_p, alpha - 1.0) * exp(-x_p)) / (K_p * pow(x_p, alpha_p - 1.0) * exp(-beta_p * x_p));
            if (uniform() <= h) {
                return x_p;
            }
        }
    }
};

struct stdlib_implementation {
    distribution<gamma_distribution<>> gamma;
    double shape;

    stdlib_implementation(double _shape, int seed) 
        : shape(_shape), gamma(shape, 1.0, seed) {}

    double operator()() {
        return gamma();
    }
};

template<class T>
int measure(T& implementation, int iterations) {
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        implementation();
    }
    auto t2 = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
}

void test_old(int iterations, double shape, int uniform_seed, int normal_seed) {
    old_implementation implementation(shape, uniform_seed, normal_seed);
    cout << "Old implementation: " << endl;
    int time = measure(implementation, iterations);
    cout << "time(mcrs): " << time << endl;
    cout << endl;
}

void test_new(int iterations, double shape, int uniform_seed) {
    new_implementation implementation(shape, uniform_seed);
    cout << "New implementation: " << endl;
    int time = measure(implementation, iterations);
    cout << "time(mcrs): " << time << endl;
    cout << endl;
}

void test_stdlib(int iterations, double shape, int gamma_seed) {
    stdlib_implementation implementation(shape, gamma_seed);
    cout << "Stdlib implementation: " << endl;
    int time = measure(implementation, iterations);
    cout << "time(mcrs): " << time << endl;
    cout << endl;
}

int main() {
    int iterations;
    double shape;
    cin >> shape >> iterations;
    int uniform_seed = 23531262;
    int normal_seed = 698536203;
    int gamma_seed = 205827155;
    test_old(iterations, shape, uniform_seed, normal_seed);
    test_new(iterations, shape, uniform_seed);
    test_stdlib(iterations, shape, gamma_seed);
    return 0;
}

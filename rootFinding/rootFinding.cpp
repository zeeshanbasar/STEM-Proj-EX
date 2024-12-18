#include <iostream>
#include <string>
#include <cmath>
#include "exprtk.hpp" // Make sure to download this header library

#include <math.h>
#include <algorithm>

using namespace std;

float const abs_tol = 1.0E-6;
float const rel_tol = 1.0E-4;
int const NMAX = 1000;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


// template <typename Derived>
template <typename T>
void newtonRhapson(string& fx, string& grad_fx, T& x, float lb, float ub){

    float f_a, f_b, f, grad_fa, grad_fb, grad_f;
    int N = 1;

    // Type definitions for the f(x)
    typedef exprtk::symbol_table<T> symbol_table_t;
    typedef exprtk::expression<T> expression_t;
    typedef exprtk::parser<T> parser_t;

    // Create a symbol table and add the variable for f(x)
    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x);
    symbol_table.add_constants();

    // Create the expression
    expression_t expr_f;
    expr_f.register_symbol_table(symbol_table);

    // Create the expression
    expression_t expr_gradf;
    expr_gradf.register_symbol_table(symbol_table);

    // Compile the expression
    parser_t parser;
    if (!parser.compile(fx, expr_f)) {
        std::cerr << "Error: Invalid f(x)." << std::endl;
        std::cerr << "Parser error: " << parser.error() << std::endl;
        return;
    }

    if (!parser.compile(grad_fx, expr_gradf)) {
        std::cerr << "Error: Invalid grad(f(x))." << std::endl;
        std::cerr << "Parser error: " << parser.error() << std::endl;
        return;
    }

    // Evaluate the expression and return the result
    float temp_x = x;

    x = lb;
    f_a = expr_f.value();
    grad_fa = expr_gradf.value();

    x = ub;
    f_b = expr_f.value();
    grad_fb = expr_gradf.value();

    x = temp_x;


    if ((f_a < 0 && f_b < 0) || (f_a > 0 && f_b > 0)) {
        if (sgn(grad_fa) == sgn(grad_fb)){
            cout << "Root does not exist in the domain. Quitting...";
            return;
        }
    }

    while (N < NMAX) {

        f = expr_f.value();        
        grad_f = expr_gradf.value();


        /* DEBUG STUFF */
        // cout << "x = " << x << "\n";
        // cout << "grad(f(x)) = " << grad_f << "\n";
        // cout << "f(x) = " << f << "\n";

        x -= f/grad_f;

        if (abs(f) < abs_tol){

            cout << "Solution found at x = " << x << "\n";
            cout << "Residual value f(x) = " << f << "\n";

            return;
        }

        if (abs(grad_f) < abs_tol){

            cout << "Gradient is zero at x = " << x << "\n";
            cout << "Method failed. Last solution at x = " << x << ", evaluated to f(x) = " << f << "\n";

            return;
        }

        N += 1;

    }

    cout << "Method failed. Last solution at x = " << x << ", evaluated to f(x) = " << f;

}


int main() {

    float lb, ub, x;    

    // Prompt for f(x)
    string fx;
    cout << "f(x) = ";
    getline(cin, fx);

    // Prompt for grad(f(x))
    string grad_fx;
    cout << "grad(f(x)) = ";
    getline(cin, grad_fx);

    while (true) {

        bool caseFlag = false;

        cout << "Domain [lb, ub] = \nlb = ";
        cin >> lb;
        cout << "ub = ";
        cin >> ub;

        while (!caseFlag){
            try {
                string input;
                cout << "\nEnter a value for x0 (or 'q' to quit): ";
                getline(cin, input);

                // Check if user wants to quit
                if (input == "q" || input == "quit") {
                    break;
                }

                // Try to convert input to a number
                try {
                    x = stod(input);
                } catch (const invalid_argument&) {
                    cerr << "Error: Invalid number. Please enter a numeric value." << endl;
                    continue;
                } catch (const out_of_range&) {
                    cerr << "Error: Number out of range." << endl;
                    continue;
                }
                if (x < lb || x > ub){
                    cout << "x0 out of bounds, try again!" << "\n";
                    cout << "x0 = ";
                    caseFlag = false;

                } else {
                    caseFlag = true;
                }
            }
            catch (bool caseFlag) {
                caseFlag = false;
            }
        }

        // Evaluate the root
        newtonRhapson(fx, grad_fx, x, lb, ub);

        break;
    }

    return 0;
}
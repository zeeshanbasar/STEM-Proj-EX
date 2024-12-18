#include <iostream>
#include <math.h>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

using Eigen::MatrixXf;
using Eigen::VectorXf;

float const abs_tol = 1.0E-6;
float const rel_tol = 1.0E-4;
int const NMAX = 1000;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template <typename Derived>
void bisection(const EigenBase<Derived>& coeffs, float x, float lb, float ub){

    VectorXf const poly_coeffs = coeffs.derived();

    float f, c, f_a, f_b, f_c;
    int N = 1;

    VectorXf const powers = VectorXf::LinSpaced(poly_coeffs.size(), 0, poly_coeffs.size()-1).reverse();

    VectorXf a = pow(lb*VectorXf::Ones(powers.size()).array(), powers.array());
    VectorXf b = pow(ub*VectorXf::Ones(powers.size()).array(), powers.array());

    f_a = poly_coeffs.dot(a);
    f_b = poly_coeffs.dot(b);

    if ((f_a < 0 && f_b < 0) || (f_a > 0 && f_b > 0)) {
        cout << "Root does not exist in the domain. Quitting...";
        return;
    }


    while (N < NMAX) {
        
        a = pow(lb*VectorXf::Ones(powers.size()).array(), powers.array());
        b = pow(ub*VectorXf::Ones(powers.size()).array(), powers.array());
        VectorXf c = pow(0.5*(lb+ub)*VectorXf::Ones(powers.size()).array(), powers.array());

        f_a = poly_coeffs.dot(a);
        f_b = poly_coeffs.dot(b);
        f_c = poly_coeffs.dot(c);

        if ((f_a < 0 && f_b < 0) || (f_a > 0 && f_b > 0)) {
            cout << "Root does not exist in the domain. Quitting...";
            return;
        }

        if (abs(f_c) < abs_tol || 0.5*(ub - lb) < rel_tol) {

            cout << "Solution found at x = " << 0.5*(lb + ub) << "\n";
            cout << "Residual f(x) = " << f_c << "\n";
            return;

        } else {
            
            N += 1;

            if (sgn(f_c) == sgn(f_a)) {
                lb = 0.5*(lb + ub);
            } else {
                ub = 0.5*(lb + ub);
            }
        }
    }

    cout << "Method failed. Last solution at x = " << c << ", evaluated to f(x) = " << f_c;

}

template <typename Derived>
void newtonRhapson(const EigenBase<Derived>& coeffs, float x, float lb, float ub){

    float f_a, f_b, f, f_prime;
    int N = 1;

    VectorXf const poly_coeffs = coeffs.derived();
    VectorXf const grad_poly_coeffs = poly_coeffs(seq(0,poly_coeffs.size()-2));

    VectorXf const powers = VectorXf::LinSpaced(poly_coeffs.size(), 0, poly_coeffs.size()-1).reverse();
    VectorXf const grad_powers = VectorXf::LinSpaced(poly_coeffs.size()-1, 0, poly_coeffs.size()-2).reverse();

    VectorXf a = pow(lb*VectorXf::Ones(powers.size()).array(), powers.array());
    VectorXf b = pow(ub*VectorXf::Ones(powers.size()).array(), powers.array());

    f_a = poly_coeffs.dot(a);
    f_b = poly_coeffs.dot(b);

    if ((f_a < 0 && f_b < 0) || (f_a > 0 && f_b > 0)) {
        cout << "Root does not exist in the domain. Quitting...";
        return;
    }

    VectorXf x_list, grad_x_list, temp;

    while (N < NMAX) {

        x -= f/f_prime;

        x_list = pow(x*VectorXf::Ones(powers.size()).array(), powers.array());

        grad_x_list = pow(x*VectorXf::Ones(grad_powers.size()).array(), grad_powers.array());

        temp = grad_x_list.array()*powers(seq(0, powers.size()-2)).array();

        f = poly_coeffs.dot(x_list);

        f_prime = grad_poly_coeffs.dot(temp);

        if (abs(f) < abs_tol || x < rel_tol){

            cout << "Solution found at x = " << x << "\n";
            cout << "Residual value f(x) = " << f;

            return;
        }
        
    }

    cout << "Method failed. Last solution at x = " << x << ", evaluated to f(x) = " << f;

}


int main(){

    int order;

    cout << "Order of polynomial, N = ";
    cin >> order;
    float in[order+1];

    for (int i = 0; i < order+1; i++) {
        cin >> in[i];
    }

    Map<VectorXf,0,InnerStride<1>> coeffs(in,order+1);

    cout << "Domain [lb, ub] = ";
    float lb, ub;
    cin >> lb >> ub;

    cout << "Choose method: 1. Bisection, 2. Newton-Rhapson: \n";
    bool caseFlag = false;
    int caseNum;
    // Get input case
    while (!caseFlag) {
        try {
            cin >> caseNum;
            if (caseNum < 1 || caseNum > 2){
                cout << "Wrong choice." << "\n";
                cout << "Enter choice: "<< "\n";
                caseFlag = false;

            } else {
                cout << "Selection: " << caseNum << "\n";
                cout << "\n";
                caseFlag = true;
            }
        }
        catch (bool caseFlag) {
            caseFlag = false;
        }
    }

    
    cout << "x0 = ";
    float x;
    caseFlag = false;

    while (!caseFlag){
        try {
            cin >> x;
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
    
    switch (caseNum)
    {
    case 1:

        bisection(coeffs, x, lb, ub);

        break;

    case 2:

        newtonRhapson(coeffs, x, lb, ub);
    
    default:
        break;
    }


    return 0;
}
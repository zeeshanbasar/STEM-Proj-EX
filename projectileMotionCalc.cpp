#include <iostream>
#include <math.h>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const float g = -9.80665;
const int N = 1000;
const double PI = 3.141592653589793;


using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(){

    cout << "Projectile Motion Calculator" << "\n";
    cout << "Assuming initial height (y0) and final height (yf) known, specify the other known variables:" << "\n";
    cout << "\t1. Initial Velocity (V0), Launch Angle (α) " << "\n";
    cout << "\t2. Horizontal Distance (l), Max. Height (h) " << "\n";
    cout << "\t3. Initial Velocity (V0), Flight Duration (t) " << "\n";
    cout << "\t4. Initial Velocity (V0), Maximum Height (h) " << "\n";
    cout << "\t5. Flight Duration (t), Horizontal Distance (l) " << "\n";
    cout << "\t6. Initial Velocity (V0), Horizontal Distance (l) " << "\n";
    cout << "\t7. Launch Angle (α), Horizontal Distance (l) " << "\n";
    cout << "\t8. Launch Angle (α), Flight Duration (t) " << "\n";
    cout << "\t9. Launch Angle (α), Maximum Height (h) " << "\n";
    cout << "Enter choice (1-9): ";
    
    int caseNum;
    bool caseFlag = false;

    // Get input case
    while (!caseFlag) {
        try {
            cin >> caseNum;
            if (caseNum < 1 || caseNum > 9){
                cout << "Wrong choice. Enter a number between 1 to 9 to specify the known parameters." << "\n";
                cout << "Enter choice (1-9): "<< "\n";
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

    float tf1, tf2, tf;

    float x0;
    cout << "Enter x0 (m):";
    cin >> x0;

    float y0;
    cout << "Enter y0 (m):";
    cin >> y0;

    float yf;
    cout << "Enter yf (m):";
    cin >> yf; 

    // Solve for each case
    switch (caseNum)
    {
    case 1: // 1. Initial Velocity (V0), Launch Angle (α)
    {

        float V0;
        cout << "Enter V0 (m/s):";
        cin >> V0;
        
        float alpha;
        cout << "Enter α (deg):";
        cin >> alpha;
        alpha *= PI/180.0;

        const float vx0 = cos(alpha)*V0;
        const float vy0 = sin(alpha)*V0;

        // solve for tf
        tf1 = (-vy0 + sqrt(vy0*vy0 - 2*g*(y0 - yf)))/g;
        tf2 = (-vy0 - sqrt(vy0*vy0 - 2*g*(y0 - yf)))/g;


        if ((tf1 <= 0) && (tf2 > 0)){
            tf = tf2;
            cout << "tf = " << tf << "\n";
        } else if ((tf1 > 0) && (tf2 <= 0)) {
            tf = tf1;
            cout << "tf = " << tf << "\n";
        } else {
            cout << "No valid solutions found for final time (tf). Quitting...";
            tf = 0;
            // cout << "tf1 = " << tf1;
            // cout << "tf2 = " << tf2;
        }

        VectorXf t = VectorXf::LinSpaced(N, 0, tf);

        VectorXf vx = VectorXf::Ones(N)*vx0;

        VectorXf vy = VectorXf::Ones(N)*vy0 + g*t;

        VectorXf x = vx0*t;

        VectorXf squared_time = t.array().square();

        VectorXf y = VectorXf::Ones(t.size())*y0 + vy0*t + 0.5*g*squared_time;

        VectorXf temp = vy.array()*vx.array().inverse();
        VectorXf alpha_t = (temp.array().atan());

        cout << "xf = " << x(lastN(1)) << "\n";
        cout << "yf = " << y(lastN(1)) << "\n";
        cout << "vxf = " << vx(lastN(1)) << "\n";
        cout << "vyf = " << vy(lastN(1)) << "\n";
        cout << "alpha_f = " << alpha_t((lastN(1)))*180.0/PI << "\n";
        Index maxRow_h, maxCol_h;
        float h_max = y.maxCoeff(&maxRow_h, &maxCol_h);   
        cout << "h_max = " << h_max << "\n";
        float t_hMax = t[maxRow_h];
        cout << "time to max h = " << t_hMax << "\n";


}
        break;

    case 2: // 2. Horizontal Distance (l), Max. Height (h)
{
        float l;
        cout << "Enter l (m):";
        cin >> l;

        float hMax;
        cout << "Enter h_max (m):";
        cin >> hMax;

        // solve for vy0
        float vy0 = sqrt(abs(2*g*(y0 - hMax)));
        

        // solve for tf
        tf1 = 2*(-vy0 + sqrt(vy0*vy0 - 2*g*(y0 - hMax)))/g;
        tf2 = 2*(-vy0 - sqrt(vy0*vy0 - 2*g*(y0 - hMax)))/g;

        if ((tf1 <= 0) && (tf2 > 0)){
            tf = tf2;
            cout << "tf = " << tf << "\n";
        } else if ((tf1 > 0) && (tf2 <= 0)) {
            tf = tf1;
            cout << "tf = " << tf << "\n";
        } else if ((tf1 > 0) && (tf2 > 0) && (tf1 == tf2)){
            tf = tf1;
            cout << "tf = " << tf << "\n";
        } else {
            cout << "No valid solutions found for final time (tf). Quitting...";
            tf = 0;
            cout << "tf1 = " << tf1;
            cout << "tf2 = " << tf2;
        }

        // solve for vx0
        float vx0 = (l - x0)/tf;

        float alpha = atan2(vy0, vx0);

        VectorXf t = VectorXf::LinSpaced(N, 0, tf);

        VectorXf vx = VectorXf::Ones(N)*vx0;

        VectorXf vy = VectorXf::Ones(N)*vy0 + g*t;

        VectorXf x = vx0*t;

        VectorXf squared_time = t.array().square();

        VectorXf y = VectorXf::Ones(t.size())*y0 + vy0*t + 0.5*g*squared_time;

        VectorXf temp = vy.array()*vx.array().inverse();
        VectorXf alpha_t = (temp.array().atan());

        cout << "xf = " << x(lastN(1)) << "\n";
        cout << "yf = " << y(lastN(1)) << "\n";
        cout << "vx0 = " << vx0 << "\n";
        cout << "vy0 = " << vy0 << "\n";
        cout << "vxf = " << vx(lastN(1)) << "\n";
        cout << "vyf = " << vy(lastN(1)) << "\n";
        cout << "alpha_f = " << alpha_t((lastN(1)))*180.0/PI << "\n";
        Index maxRow_h, maxCol_h;
        float h_max = y.maxCoeff(&maxRow_h, &maxCol_h);   
        cout << "h_max = " << h_max << "\n";
        float t_hMax = t[maxRow_h];
        cout << "time to max h = " << t_hMax << "\n";

}

        break;

    // case 3:
    //     /* code */
    //     break;

    // case 4:
    //     /* code */
    //     break;
    
    // case 5:
    //     /* code */
    //     break;

    // case 6:
    //     /* code */
    //     break;

    // case 7:
    //     /* code */
    //     break;

    // case 8:
    //     /* code */
    //     break;
    
    // case 9:
    //     /* code */
    //     break;

    // default:
    //     break;
    }

    return 0;
}



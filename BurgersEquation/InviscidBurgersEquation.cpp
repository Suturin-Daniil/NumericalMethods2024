#include <iostream>
#include <functional>
#include <type_traits>
#include <eigen-3.4.0/Eigen/Dense> 
#include <fstream>

float u0(float x){ // initial condition
    return cos(x*M_PI);
}
float f(float u){
    return pow(u,2)/2;
}
// Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> timeDerivative(float t, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> u){
//     unsigned int N = u.size();
// }
int main(){
    const unsigned int M = 100; // number of nodes
    float h = 1./100.; // spatial step
    float t = 0.1; // time step
    float T = 50.;
    float LeftBC = 0;
    float RightBC = 0;

    Eigen::Matrix<float, M+1, 1> x; // Grid points
    x.setZero();
    for (unsigned int i = 0; i < M+1; i++){
        x(i) = h*i;
    }

    Eigen::Matrix<float, M, 1> xc; // Central points
    xc.setZero();
    for (unsigned int i = 0; i < M; i++){
        xc(i) = (x(i+1) + x(i))/2;
    }

    float time = 0;
    Eigen::Matrix<float, M+2, 1> result;
    result.setZero();
    Eigen::Matrix<float, M+2, 1> previous;
    previous.setZero();
    Eigen::Matrix<float, M+2, 1> U0;
    for (unsigned int i = 1; i < M+1; i++){
        U0(i) = u0(x(i));
    }
    U0(0) = 0;
    U0(M+1) = 0;

    Eigen::Matrix<float, M + 1, 1> f_interior;
    f_interior.setZero();

    std::ofstream Solution("BurgersSolution.txt");
    Solution << "StepX " << h << " StepT " << t << std::endl; 

    while (time <= T){ // time loop
        if (time == 0){
            for (unsigned int i = 0; i < M+1; i++){
                if (U0(i) >= U0(i+1)){
                    f_interior(i) = std::max(f(i*h), f((i+1)*h));
                }
                else if (U0(i) <= 0 && 0 <= U0(i+1)){
                    f_interior(i) = 0;
                }
                else{
                    f_interior(i) = std::min(f(i*h), f((i+1)*h));
                }
            }
            for (unsigned int i = 1; i < M+1; i++){
                result(i) = U0(i) + t/h*(f_interior(i-1) - f_interior(i));
            }
            result(0) = 0;
            result(M+1) = 0;
        }

        else if (time != 0){
           for (unsigned int i = 0; i < M+1; i++){
                if (previous(i) >= previous(i+1)){
                    f_interior(i) = std::max(f(i*h), f((i+1)*h));
                }
                else if (previous(i) <= 0 && 0 <= previous(i+1)){
                    f_interior(i) = 0;
                }
                else{
                    f_interior(i) = std::min(f(i*h), f((i+1)*h));
                }
            }
            for (unsigned int i = 1; i < M+1; i++){
                result(i) = previous(i) + t/h*(f_interior(i-1) - f_interior(i));
            }
            result(0) = 0;
            result(M+1) = 0; 
        }
        previous = result;
        time += t;
        for (unsigned int i = 0; i < result.size(); i++){
            if (i == result.size()-1){Solution << result(i) << std::endl;}
            else{Solution << result(i) << " ";}         
        }
    }
}
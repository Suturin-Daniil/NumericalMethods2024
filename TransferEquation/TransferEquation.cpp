#include <iostream>
#include <functional>
#include <type_traits>
#include <eigen-3.4.0/Eigen/Dense> 
#include <fstream>

template <typename Xtype> // initial condition
Xtype u0(Xtype x, Xtype L){
    return sin(4*M_PI*x/L);
}

template <typename Xtype, typename Ttype> // free function
Xtype F(Xtype x, Ttype t){
    return 0;
}

template <typename Xtype, typename Ttype>
struct BoundaryCondition{
    Xtype CFL = 1.00;

    Xtype a = 1;

    Xtype Xmax = 20; // size of computational domain
    Xtype Xmin = 0;
    Ttype Tmax = 20;
    Ttype Tmin = 0;

    std::function<Xtype(Xtype, Xtype)> U0 = u0<Xtype>; // function of initial condition

    // std::function<Xtype(Xtype, Ttype)> South = S<Xtype, Ttype>; // 4 functions of boundary condition
    // std::function<Xtype(Xtype, Ttype)> North = N<Xtype, Ttype>; // for rectangular computational domain
    // std::function<Xtype(Xtype, Ttype)> West = W<Xtype, Ttype>;
    // std::function<Xtype(Xtype, Ttype)> East = E<Xtype, Ttype>;

    std::function<Xtype(Xtype, Ttype)> f = F<Xtype, Ttype>; // free function
};

template <typename Xtype, typename Ttype>
void CornerSolver(
Xtype & stepx, 
Ttype & stept, 
struct BoundaryCondition<Xtype, Ttype> & BC){
    unsigned int M = (unsigned int)((BC.Xmax - BC.Xmin)/stepx) + 1;
    Xtype h = (BC.Xmax - BC.Xmin)/M;
    //N = (unsigned int)((BC.Tmax - BC.Tmin)/stept) + 1;
    //t = (BC.Tmax - BC.Tmin)/N;
    Xtype t = BC.CFL*h/BC.a;

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> result;
    result.resize(M+1, 1);
    result.setZero();

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> previous;
    previous.resize(M+1, 1);
    previous.setZero();

    Ttype time = BC.Tmin;

    std::ofstream solution("CornerSolution0.1h.txt");

    solution << "StepX " << h << " Xmax " << BC.Xmax << " StepT " << t << std::endl; 
    
    Xtype L = BC.Xmax - BC.Xmin;

    //time loop
    while (time <= BC.Tmax)
    {
        if (time == BC.Tmin){
        for (unsigned int i = 1; i <= M; i++){
            result(i) = BC.f(BC.Xmin+i*h, time) + BC.U0(BC.Xmin + i*h, L)*(1 - BC.CFL) + BC.CFL*BC.U0(BC.Xmin + (i-1)*h, L);
            if (i == M) {result(0) = result(i);}
        }
        }

        if (time != BC.Tmin){
            for (unsigned int i = 1; i <= M; i++){
            result(i) = BC.f(BC.Xmin+i*h, time) + previous(i)*(1 - BC.CFL) + BC.CFL*previous(i-1);
            if (i == M) {result(0) = result(i);}
        }
        }

        time += t;
        for (unsigned int i = 0; i <= M; i++){
            if (i != M) {solution << result(i) << " ";}
            if (i == M) {solution << result(i) << std::endl;}
        }
        previous = result;
        
    }
    solution.close();
}

template <typename Xtype, typename Ttype>
void LaxWendroff(
Xtype & stepx, 
Ttype & stept, 
struct BoundaryCondition<Xtype, Ttype> & BC){
    unsigned int M = (unsigned int)((BC.Xmax - BC.Xmin)/stepx) + 1;
    Xtype h = (BC.Xmax - BC.Xmin)/M;
    //N = (unsigned int)((BC.Tmax - BC.Tmin)/stept) + 1;
    //t = (BC.Tmax - BC.Tmin)/N;
    Xtype t = BC.CFL*h/BC.a;

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> result;
    result.resize(M+1, 1);
    result.setZero();

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> previous;
    previous.resize(M+1, 1);
    previous.setZero();

    Ttype time = BC.Tmin;

    Xtype L = BC.Xmax - BC.Xmin;

    std::ofstream solution("LaxSolution0.1h.txt");

    solution << "StepX " << h << " Xmax " << BC.Xmax << " StepT " << t << std::endl; 

    //time loop
    while(time <= BC.Tmax){
        if (time == BC.Tmin){
            for (unsigned int i = 0; i <= M; i++){
                if (i == 0){result(i) = BC.f(BC.Xmin + i*h, time) + BC.U0(BC.Xmin, L) - BC.CFL/2*(BC.U0(BC.Xmin + h, L) - BC.U0(BC.Xmin + M*h, L)) + pow(BC.CFL, 2)/2*(BC.U0(BC.Xmin + h, L) - 2*BC.U0(BC.Xmin, L) + BC.U0(BC.Xmin + M*h, L));}
                else if (i == M){result(i) = BC.f(BC.Xmin + i*h, time) + BC.U0(BC.Xmin + M*h, L) - BC.CFL/2*(BC.U0(BC.Xmin, L) - BC.U0(BC.Xmin + (M-1)*h, L)) + pow(BC.CFL, 2)/2*(BC.U0(BC.Xmin, L) - 2*BC.U0(BC.Xmin + M*h, L) + BC.U0(BC.Xmin + (M-1)*h, L));}
                else{result(i) = BC.f(BC.Xmin + i*h, time) + BC.U0(BC.Xmin + i*h, L) - BC.CFL/2*(BC.U0(BC.Xmin + (i+1)*h, L) - BC.U0(BC.Xmin + (i-1)*h, L)) + pow(BC.CFL, 2)/2*(BC.U0(BC.Xmin + (i+1)*h, L) - 2*BC.U0(BC.Xmin + i*h, L) + BC.U0(BC.Xmin + (i-1)*h, L));}
            }
        }

        if (time != BC.Tmin){
            for (unsigned int i = 0; i <= M; i++){
                if (i == 0){result(i) = BC.f(BC.Xmin + i*h, time) + previous(i) - BC.CFL/2*(previous(i+1) - previous(M)) + pow(BC.CFL, 2)/2*(previous(1) - 2*previous(i) + previous(M));}
                else if (i == M){result(i) = BC.f(BC.Xmin + i*h, time) + previous(M) - BC.CFL/2*(previous(0) - previous(M-1)) + pow(BC.CFL, 2)/2*(previous(0) - 2*previous(M) + previous(M-1));}
                else{result(i) = BC.f(BC.Xmin + i*h, time) + previous(i) - BC.CFL/2*(previous(i+1) - previous(i-1)) + pow(BC.CFL, 2)/2*(previous(i+1) - 2*previous(i) + previous(i-1));}
            }
        }

        time += t;
        for (unsigned int i = 0; i <= M; i++){
        if (i != M) {solution << result(i) << " ";}
        if (i == M) {solution << result(i) << std::endl;}
        }
        previous = result;
    }
    solution.close();

}

int main(){
    BoundaryCondition<float, float> BC;

    float stepx = 0.11;
    float stept = 0.11;

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> res;
    res.setZero();

    LaxWendroff<float, float>(stepx, stept, BC);
}


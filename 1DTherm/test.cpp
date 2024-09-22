#include <vector>
#include <iostream>
#include <functional>
#include <type_traits>
#include <eigen-3.4.0/Eigen/Dense> 
#include <fstream>

template<typename vtype>
void printVec(std::vector<vtype> &vec){
    for (vtype el: vec) {
        std::cout << el << " ";
    }
    std::cout << std::endl;
}

template <typename Ttype> // function of left boundary condition
Ttype g0(Ttype t){
    return sqrt(0.1/t)*exp(-25/(4*t));
}

template <typename Ttype> // function of right boundary condition
Ttype g1(Ttype  t){
    return sqrt(0.1/t)*exp(-25/(4*t));
}

template <typename Xtype> // initial condition
Xtype u0(Xtype x){
    return exp(-pow(x,2)/0.4);
}

template <typename Xtype, typename Ttype> // free function
Xtype F(Xtype  x, Ttype  t){
    return 0;
}

template <typename Xtype, typename Ttype>
struct BoundaryCon{
    float a;

    float Xmax; // size of computational domain
    float Xmin;
    float Tmax;
    float Tmin;

    float a0;
    float b0;
    std::function<Ttype(Ttype)> J0 = g0<Ttype>; // left boundary condition

    float a1;
    float b1;
    std::function<Ttype(Ttype )> J1 = g1<Ttype>; // right boundary condition

    std::function<Xtype(Xtype )> U0 = u0<Xtype>; // initial condition

    std::function<Xtype(Xtype , Ttype )> f = F<Xtype, Ttype>; // free function
};

// Three diagonal matrix class. Construct matrix
template<typename Xtype, typename Ttype>
class ThreeDiagonalMatrix {
    public:
    unsigned int M; // 
    unsigned int N; // 
    Xtype h;
    Ttype t;
    Xtype sigma;

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> A;

    ThreeDiagonalMatrix(Xtype & steph, Ttype & stept, struct BoundaryCon<Xtype, Ttype> & BounCon){
        M = (unsigned int)((BounCon.Xmax - BounCon.Xmin)/steph) + 1;
        h = (BounCon.Xmax - BounCon.Xmin)/M;
        N = (unsigned int)((BounCon.Tmax - BounCon.Tmin)/stept) + 1;
        t = (BounCon.Tmax - BounCon.Tmin)/N;

        sigma = (float)(BounCon.a*t)/(2*pow(h, 2)); // Courant number

        A.resize(M+1, M+1);
        A.setIdentity();
        A *= (1+2*sigma);
        A(0, 0) = BounCon.a0 - (3*BounCon.b0)/(2*this->h);
        A(M, M) = BounCon.a1 + (3*BounCon.b1)/(2*this->h);
        for (unsigned int i = 0; i <= M; i++){
            for (unsigned int j = 0; j <= M; j++){
                if (j-i == 1 || i-j == 1) {A(i,j) = -sigma;}
            }
        }
        A(0,1) = (2*BounCon.b0)/h;
        A(0,2) = -(BounCon.b0)/(2*h);
        A(M, M-1) = -(2*BounCon.b1)/h;
        A(M, M-2) = -(BounCon.b1)/(2*h);

        // A.row(0) = A.row(0) + A.row(1)*(A(0,2)/sigma);
        // A.row(M) = A.row(M) + A.row(M-1)*(A(M, M-2)/sigma);
        // A(0,2) = 0;
        // A(M, M-2) = 0;
    }; 

    void PrintMatrix(void) {
        std::cout << A << std::endl;
    }

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> GetMatrix(void){
        return A;
    }
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

// Функция для решения методм прогонки 
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(ThreeDiagonalMatrix<mType, mType>& matrix, const std::vector<cType>& column){
    unsigned int size = column.size();
    std::vector<DivisType<cType, mType>> alpha((size - 1), 0);
    std::vector<DivisType<cType, mType>> beta(size, 0);
    std::vector<DivisType<cType, mType>> y(size, 0);
    std::vector<DivisType<cType, mType>> m(size, 0);
    Eigen::Matrix<mType, Eigen::Dynamic, Eigen::Dynamic> mat = matrix.GetMatrix();

    y[0] = mat(0,0);
    alpha[0] = (-1)*(mat(0,1) / y[0]);
    beta[0] = column[0] / y[0];

    for (int j = 1; j < size - 1; j++) {
        y[j] = mat(j,j) + alpha[j - 1] * mat(j,j - 1);
        alpha[j] = - mat(j,j + 1) / y[j];
        beta[j] = (column[j] - mat(j,j - 1) * beta[j - 1]) / y[j];
    }

    y[size - 1] = mat(size - 1,size - 1) + mat(size - 1,size - 2) * alpha[size - 2];
    beta[size - 1] = (column[size - 1] - mat(size - 1,size - 2) * beta[size - 2]) / y[size - 1];

    m[size - 1] = beta[size - 1];

    for (int i = size - 2; i >= 0; i--){
        m[i] = alpha[i] * m[i + 1] + beta[i];
    }
    return m;
}

int main(){
    struct BoundaryCon<float, float> BC;
    BC.a = 1;
    BC.a0 = 1;
    BC.a1 = 1;
    BC.b0 = 0;
    BC.b1 = 0;
    BC.Tmax = 5;
    BC.Tmin = 0.1;
    BC.Xmax = 5;
    BC.Xmin = -5;
    float t = 0.01;
    float h = 0.05;
    ThreeDiagonalMatrix<float, float> TriDigMat(h, t, BC);

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> A = TriDigMat.GetMatrix();

    unsigned int sizeOfMat = A.rows(); // system matrix dimension

    //Create free column of the system
    std::vector<float> d(sizeOfMat, 0);
    printVec(d);

    std::vector<DivisType<float, float>> result(sizeOfMat,0);

    std::ofstream Solution("1DThermSolutionH0.05.txt");

    std::ofstream SolutionT("1DThermT5.txt");
    std::ofstream SolutionX("1DThermXH0.05.txt");

    for (unsigned int k = 0; k <= TriDigMat.M; k++){ // write x coordinate
        if (k != TriDigMat.M) {SolutionX << BC.Xmin + k*TriDigMat.h << " ";}
        if (k == TriDigMat.M) {SolutionX << BC.Xmin + k*TriDigMat.h;}
    }

    for (unsigned int k = 0; k <= TriDigMat.N; k++){ // write t coordinate
        if (k != TriDigMat.N) {SolutionT << BC.Tmin + k*TriDigMat.t << " ";}
        if (k == TriDigMat.N) {SolutionT << BC.Tmin + k*TriDigMat.t;}
    }

    for (unsigned int k = 0; k <= TriDigMat.M; k++){
       if (k != TriDigMat.M) {Solution << BC.U0((float)(BC.Xmin + k*TriDigMat.h)) << " ";}
       if (k == TriDigMat.M) {Solution << BC.U0((float)(BC.Xmin + k*TriDigMat.h)) << std::endl;}
    }
    
    for (double time = BC.Tmin; time <= BC.Tmax; time += TriDigMat.t){ // time step loop

        for (unsigned int i = 0; i < sizeOfMat; i++){
            if (i == 0){d[0] = BC.J0((float)(time + TriDigMat.t));}
            if (i == sizeOfMat-1){d[i] = BC.J1((float)(time + TriDigMat.t));}
            else if (i != 0 & i != sizeOfMat-1){
                if (time == BC.Tmin){
                    d[i] = TriDigMat.sigma*(BC.U0((float)(BC.Xmin + (i - 1)*TriDigMat.h)) - 2*BC.U0((float)(BC.Xmin + i*TriDigMat.h)) + BC.U0((float)(BC.Xmin + (i+1)*TriDigMat.h))) + BC.U0((float)(BC.Xmin + i*TriDigMat.h)) + TriDigMat.t/2*(BC.f(BC.Xmin + i*TriDigMat.h, time) + BC.f(BC.Xmin + i*TriDigMat.h, time + TriDigMat.t));
                }
                if (time != BC.Tmin){
                    d[i] = TriDigMat.sigma*(result[i-1] - 2*result[i] + result[i+1]) + result[i] + TriDigMat.t/2*(BC.f(BC.Xmin + i*TriDigMat.h, time) + BC.f(BC.Xmin + i*TriDigMat.h, time + TriDigMat.t));
                }
            }
        }
    
        A.row(0) = A.row(0) + A.row(1)*(A(0,2)/TriDigMat.sigma); //
        A.row(sizeOfMat - 1) = A.row(sizeOfMat - 1) + A.row(sizeOfMat-2)*(A(sizeOfMat - 1, sizeOfMat - 3)/TriDigMat.sigma);
        A(0,2) = 0;
        A(sizeOfMat - 1, sizeOfMat - 3) = 0;

        d[0] = d[0] + d[1]*(A(0, 2)/TriDigMat.sigma);
        d[TriDigMat.M] = d[TriDigMat.M] + d[TriDigMat.M - 1] * (A(sizeOfMat - 1, sizeOfMat - 3)/TriDigMat.sigma);

        result = solve<float, float>(TriDigMat, d);
        for (unsigned int l = 0; l <= TriDigMat.M; l++){
            if (l == TriDigMat.M) {Solution << result[l] << std::endl;}
            else {Solution << result[l] << " ";}
        }
        
    }

    Solution.close();
    SolutionT.close();
    SolutionX.close();
}
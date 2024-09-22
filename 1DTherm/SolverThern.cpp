#include <vector>
#include <iostream>
#include <functional>
#include <type_traits>
#include <eigen-3.4.0/Eigen/Dense> 


template<typename Xtype, typename Ttype>
Xtype f(Xtype & x, Ttype & t){
    return;
}

template<typename Ttype> // function of boundary condition J0(t) = a0*U(t,A) + b0*dU/dx(t,A)
Ttype J0(Ttype & t){ 
    return;
}

template<typename Ttype> 
Ttype J1(Ttype & t){ // J1(t) = a1*U(t,B) + b1*dU/dx(t,B)
    return;
}


template<typename Xtype> // function of initial condition U(0, x) = U0(x)
Xtype U0(Xtype & x){
    return ;
}


template<typename Type, typename Xtype, typename Ttype>
struct BoundaryCon{
    Type a0;
    Type b0;
    std::function<Ttype(Ttype)> g0 = J0;
    Type a1;
    Type b1;
    std::function<Xtype(Ttype)> g1 = J1;
    Xtype Xmax;
    Xtype Xmin;
    Ttype Tmax;
    Ttype Tmin;
    std::function<Xtype(Xtype)> n = U0;
    std::function<Xtype(Xtype, Ttype)> F = f;
};


template<typename Type, typename Xtype, typename Ttype>
class ThreeDiagonalMatrix {
    unsigned int M;
    unsigned int N;
    Xtype h;
    Ttype t;
    Eigen::Matrix<double, -1,-1> A;

    ThreeDiagonalMatrix(Xtype &steph, Ttype &stept, struct BoundaryCon<Type, Xtype, Ttype> & BounCon){
        M = (unsigned int)((BounCon.Xmax - BounCon.Xmin)/steph) + 1;
        h = (BounCon.Xmax - BounCon.Xmin)/M;
        N = (unsigned int)((BounCon.Tmax - BounCon.Tmin)/stept) + 1;
        t = (BounCon.Tmax - BounCon.Tmin)/N;

        this->A = Eigen::Matrix<double, this->N, this->M> B;
        A.setZero();
    } 

    void PrintMatrix(void) {
        std::cout << A << std::endl;
    }
};


int main(){
// Here I just try to use Eigen library in order to use Eigen matrix operation in code above in ThreeDiagonalMatrix class
Eigen::Matrix<float, 3, 3> mat;
mat.setZero();
std::cout << mat << std::endl;

Eigen::Matrix3f matf;
matf.setOnes();
std::cout << matf << std::endl;

Eigen::MatrixXf A(10,10);
A.setZero();
for (int i = 0; i < 10; i++){
    for (int j = 0; j < 10; j++){
        if (i == j){A(i, j) = 1;}
        if (abs(i-j) == 1){A(i,j) = 5;}
    }
}
std::cout << A << std::endl;
std::cout << std::endl;
A.row(0) = A.row(0) - A.row(1);
std::cout << A << std::endl;

Eigen::MatrixXf B(4,4);
B << 1, 2, 3, 4,
    5, 6, 7, 8,
    8, 9, 10, 11,
    12, 13, 14, 15;
std::cout << B.block(0,0,3,3) << std::endl;

std::cout << B(1,1);

BoundaryCon<double, double, double> BC;
BC.a0=1.7;

ThreeDiagonalMatrix<double, double, double> A(0.3, 0.7, BC);
A.PrintMatrix();

return 0;
}

/**template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

// Функция для решения методм  прогонки 
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(ThreeDiagonalMatrix<mType>& matrix, const std::vector<cType>& column){
    unsigned int size = column.size();
    std::vector<DivisType<cType, mType>> alpha((size - 1), 0);
    std::vector<DivisType<cType, mType>> beta(size, 0);
    std::vector<DivisType<cType, mType>> y(size, 0);
    std::vector<DivisType<cType, mType>> m(size, 0);
    std::vector<std::vector<mType>> mat = matrix.ReturnMatrix();

    y[0] = mat[0][0];
    alpha[0] = (-1)*(mat[0][1] / y[0]);
    beta[0] = column[0] / y[0];

    for (int j = 1; j < size - 1; j++) {
        y[j] = mat[j][j] + alpha[j - 1] * mat[j][j - 1];
        alpha[j] = - mat[j][j + 1] / y[j];
        beta[j] = (column[j] - mat[j][j - 1] * beta[j - 1]) / y[j];
    }

    y[size - 1] = mat[size - 1][size - 1] + mat[size - 1][size - 2] * alpha[size - 2];
    beta[size - 1] = (column[size - 1] - mat[size - 1][size - 2] * beta[size - 2]) / y[size - 1];

    m[size - 1] = beta[size - 1];

    for (int i = size - 2; i >= 0; i--){
        m[i] = alpha[i] * m[i + 1] + beta[i];
    }
    return m;
}**/
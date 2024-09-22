#include <vector>
#include <iostream>
#include <functional>
#include <type_traits>
#include <eigen-3.4.0/Eigen/Dense> 
#include <fstream>

template <typename Xtype> // initial condition
Xtype u0(Xtype x, Xtype y){
    return exp(-pow((x-0.5),2)/2 - pow((y-0.5),2)/2);
}

template <typename Xtype, typename Ttype> // function of south sboundary condition 
Xtype S(Xtype x, Ttype t){
    return 10;
}

template <typename Xtype, typename Ttype> // function of west boundary condition 
Xtype W(Xtype y, Ttype t){
    return 10;
}

template <typename Xtype, typename Ttype> // function of east sboundary condition 
Xtype E(Xtype y, Ttype t){
    return 10;
}

template <typename Xtype, typename Ttype> // function of north sboundary condition 
Xtype N(Xtype x, Ttype t){
    return 10;
}

template <typename Xtype, typename Ttype> // free function
Xtype F(Xtype x, Xtype y, Ttype t){
    return 0;
}

template <typename Xtype, typename Ttype>
struct BoundaryCondition{
    Xtype a = 1;

    Xtype Xmax = 1; // size of computational domain
    Xtype Xmin = 0;
    Xtype Ymax = 1;
    Xtype Ymin = 0;
    Ttype Tmax = 1;
    Ttype Tmin = 0;

    std::function<Xtype(Xtype, Xtype)> U0 = u0<Xtype>; // function of initial condition

    std::function<Xtype(Xtype, Ttype)> South = S<Xtype, Ttype>; // 4 functions of boundary condition
    std::function<Xtype(Xtype, Ttype)> North = N<Xtype, Ttype>; // for rectangular computational domain
    std::function<Xtype(Xtype, Ttype)> West = W<Xtype, Ttype>;
    std::function<Xtype(Xtype, Ttype)> East = E<Xtype, Ttype>;

    std::function<Xtype(Xtype, Xtype, Ttype)> f = F<Xtype, Ttype>; // free function

};

template <typename Xtype, typename Ttype>
class FiveDiagonalMatrix{
    private:
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> Upper;
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> Lower;
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> Diagonal;

    public:
    unsigned int M; // number of nodes in x direction
    unsigned int N; // number of time step
    unsigned int K; // number of nodes in y direction
    Xtype hx, hy;
    Ttype t;
    Xtype sigmax, sigmay;

    FiveDiagonalMatrix(Xtype & stepx, Xtype & stepy, Ttype & stept, struct BoundaryCondition<Xtype, Ttype> & BC){
        M = (unsigned int)((BC.Xmax - BC.Xmin)/stepx) + 1;
        hx = (BC.Xmax - BC.Xmin)/M;
        N = (unsigned int)((BC.Tmax - BC.Tmin)/stept) + 1;
        t = (BC.Tmax - BC.Tmin)/N;
        K = (unsigned int)((BC.Ymax - BC.Ymin)/stepy) + 1;
        hy = (BC.Ymax - BC.Ymin)/K;

        sigmax = (Xtype)( (BC.a*t) / (2*pow(hx,2)) );
        sigmay = (Xtype)( (BC.a*t) / (2*pow(hy,2)) );

        A.resize((M-1)*(K-1), (M-1)*(K-1));
        Lower.resize((M-1)*(K-1), (M-1)*(K-1));
        Diagonal.resize((M-1)*(K-1), 1);
        A.setZero();
        Lower.setZero();
        Diagonal.setZero();

        A.topRightCorner((M-1)*(K-2), (M-1)*(K-2)).setIdentity();
        A.topRightCorner((M-1)*(K-2), (M-1)*(K-2)) *= -sigmay;
        Upper = A;

        A.bottomLeftCorner((M-1)*(K-2), (M-1)*(K-2)).setIdentity();
        A.bottomLeftCorner((M-1)*(K-2), (M-1)*(K-2)) *= -sigmay;

        Lower.bottomLeftCorner((M-1)*(K-2), (M-1)*(K-2)).setIdentity();
        Lower.bottomLeftCorner((M-1)*(K-2), (M-1)*(K-2)) *= -sigmay;

        for (unsigned int i = 0; i < (M-1)*(K-1); i++){
            A(i, i) = 1 + 2*sigmax + 2*sigmay;
            Diagonal(i) = 1 + 2*sigmax + 2*sigmay;
        }

        int k = 1;
        for (unsigned int i = 0; i < (M-1)*(K-1); i++){
            if (i == 0){A(i, i+1) = -sigmax; Upper(i, i+1) = -sigmax;}
            else if (i == (M-1)*(K-1) - 1){A(i, i-1) = -sigmax; Lower(i, i-1) = -sigmax;}
            else{
                A(i, i+1) = -sigmax;
                Upper(i, i+1) = -sigmax;
                A(i, i-1) = -sigmax;   
                Lower(i, i-1) = -sigmax;
            }
        }
        for (unsigned int i = 0; i < (M-1)*(K-1); i++){
            if (k == M-1 && i != (M-1)*(K-1) - 1){
                A(i, i+1) = 0;
                A(i+1, i) = 0;
                Upper(i, i+1) = 0;
                Lower(i+1, i) = 0;
                k = 0;
            }
            k += 1;
        }
    }

    void PrintMatrix(void){
        std::cout << A << std::endl;
    }

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> GetMatrix(void){
        return A;
    }

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> U(void){return Upper;}

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> L(void){return Lower;}

    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> D(void){return Diagonal;}

};

template <typename Xtype, typename Ttype>
Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> GaussSeidelSolver(
FiveDiagonalMatrix<Xtype, Ttype> & FDM, // matrix of the system
Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> & column, // free vector column
Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> & initialCondition, // initial approximation
float accuracy){
    const unsigned int size = column.size();
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> mat = FDM.GetMatrix();
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> U = FDM.U();
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> L = FDM.L();
    Eigen::Matrix<Xtype, Eigen::Dynamic, Eigen::Dynamic> D = FDM.D(); // vector of diagonal elements
    
     
    Eigen::Matrix<Xtype, Eigen::Dynamic, 1> CurrentX;
    CurrentX.resize(size, 1);
    CurrentX(0) = 1000; // in order to enter in the loop below
    Eigen::Matrix<Xtype, Eigen::Dynamic, 1> PreviousX;
    PreviousX.resize(size, 1);
    PreviousX = initialCondition;

    Eigen::Matrix<Xtype, Eigen::Dynamic, 1> FreeColumn = column;

    int k = 0;
    while ((CurrentX - PreviousX).squaredNorm() >= accuracy){
        if (k != 0){PreviousX = CurrentX;}
        k += 1;
        FreeColumn = column - U*PreviousX;

        for (unsigned int i = 0; i < size; i++){
            Xtype sum = 0;
            for (unsigned int j = 0; j < i; j++){
                sum += L(i, j)*CurrentX(j);
            }
            CurrentX(i) = (FreeColumn(i) - sum)/D(i);
        }
    }
    return CurrentX;
}

int main(){
    BoundaryCondition<float, float> BC;
    float stepx = 0.21;
    float stepy = 0.34;
    float t = 0.11;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> A;
    FiveDiagonalMatrix<float, float> FDM(stepx, stepy, t, BC);
    A = FDM.GetMatrix();

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> column;
    column.resize((FDM.M-1)*(FDM.K-1), 1);

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> CurrentState;
    CurrentState.resize((FDM.M-1)*(FDM.K-1), 1);

    std::cout << A << std::endl;

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> d;
    d.resize((FDM.M-1)*(FDM.K-1), 1);
    d << 1,
         2,
         3,
         4,
         5,
         1,
         6,
         8;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> init;
    init.resize((FDM.M-1)*(FDM.K-1), 1);
    init.setZero();

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> res;
    res = GaussSeidelSolver<float, float>(FDM, d, init, 1e-5);

    std::cout << res << std::endl;
    std::cout << A.inverse()*d << std::endl;

    return 0;
}
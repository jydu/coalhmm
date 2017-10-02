#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
namespace bcs {
bpp::RowMatrix<double> getRateMatrix_C_H( double C_C, double C_H, double R ) {
bpp::RowMatrix<double> Q_C_H;
Q_C_H.resize(15, 15);
bpp::MatrixTools::fill(Q_C_H, 0);

Q_C_H(0,0) = -(R+R);
Q_C_H(0,1) = R;
Q_C_H(0,2) = R;

Q_C_H(1,1) = -(C_C+R);
Q_C_H(1,0) = C_C;
Q_C_H(1,3) = R;

Q_C_H(2,2) = -(C_H+R);
Q_C_H(2,0) = C_H;
Q_C_H(2,3) = R;

Q_C_H(3,3) = -(C_H+C_C);
Q_C_H(3,1) = C_H;
Q_C_H(3,2) = C_C;

return Q_C_H;
}
bpp::RowMatrix<double> getRateMatrix_HC( double C_HC, double R ) {
bpp::RowMatrix<double> Q_HC;
Q_HC.resize(15, 15);
bpp::MatrixTools::fill(Q_HC, 0);

Q_HC(0,0) = -(C_HC+R+R);
Q_HC(0,4) = C_HC;
Q_HC(0,1) = R;
Q_HC(0,2) = R;

Q_HC(1,1) = -(C_HC+C_HC+C_HC+R);
Q_HC(1,0) = C_HC;
Q_HC(1,13) = C_HC;
Q_HC(1,12) = C_HC;
Q_HC(1,3) = R;

Q_HC(2,2) = -(C_HC+C_HC+C_HC+R);
Q_HC(2,0) = C_HC;
Q_HC(2,5) = C_HC;
Q_HC(2,6) = C_HC;
Q_HC(2,3) = R;

Q_HC(3,3) = -(C_HC+C_HC+C_HC+C_HC+C_HC+C_HC);
Q_HC(3,1) = C_HC;
Q_HC(3,7) = C_HC;
Q_HC(3,8) = C_HC;
Q_HC(3,9) = C_HC;
Q_HC(3,10) = C_HC;
Q_HC(3,2) = C_HC;

Q_HC(4,4) = -(R);
Q_HC(4,11) = R;

Q_HC(5,5) = -(C_HC+R);
Q_HC(5,4) = C_HC;
Q_HC(5,7) = R;

Q_HC(6,6) = -(C_HC+R);
Q_HC(6,4) = C_HC;
Q_HC(6,10) = R;

Q_HC(7,7) = -(C_HC+C_HC+C_HC);
Q_HC(7,11) = C_HC;
Q_HC(7,13) = C_HC;
Q_HC(7,5) = C_HC;

Q_HC(8,8) = -(C_HC+C_HC+R+C_HC);
Q_HC(8,12) = C_HC;
Q_HC(8,5) = C_HC;
Q_HC(8,3) = R;
Q_HC(8,14) = C_HC;

Q_HC(9,9) = -(C_HC+C_HC+C_HC+R);
Q_HC(9,13) = C_HC;
Q_HC(9,14) = C_HC;
Q_HC(9,6) = C_HC;
Q_HC(9,3) = R;

Q_HC(10,10) = -(C_HC+C_HC+C_HC);
Q_HC(10,11) = C_HC;
Q_HC(10,12) = C_HC;
Q_HC(10,6) = C_HC;

Q_HC(11,11) = -(C_HC);
Q_HC(11,4) = C_HC;

Q_HC(12,12) = -(C_HC+R);
Q_HC(12,4) = C_HC;
Q_HC(12,10) = R;

Q_HC(13,13) = -(C_HC+R);
Q_HC(13,4) = C_HC;
Q_HC(13,7) = R;

Q_HC(14,14) = -(C_HC+R+R);
Q_HC(14,4) = C_HC;
Q_HC(14,9) = R;
Q_HC(14,8) = R;

return Q_HC;
}
}

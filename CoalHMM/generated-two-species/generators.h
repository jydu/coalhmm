#ifndef GENERATORS_H
#define GENERATORS_H
namespace bcs {
bpp::RowMatrix<double> getRateMatrix_C_H( double C_C, double C_H, double R );
bpp::RowMatrix<double> getRateMatrix_HC( double C_HC, double R );
}
#endif

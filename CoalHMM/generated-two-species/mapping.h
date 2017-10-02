
#ifndef MAPPING_H
#define MAPPING_H
#include <vector>
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace gf {
void project_HC(bpp::Matrix<double>& m);
void expand_HC(const std::vector<double>& sep_system, std::vector<double> &mi_system);
}
#endif

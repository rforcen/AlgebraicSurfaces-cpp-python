#include "AlgebraicSurfaces.h"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;

#include "python.h"

static p::list algebraic_surface(int n_func, int resolution) {
    SurfaceMesh sm;
    Mesh mesh = sm.generate_mesh(n_func, resolution);
    return mesh2list(mesh);
}

static p::object func_name(int n_func) {
    SurfaceMesh sm;
    return p::object(sm.get_name(n_func));
}

static p::list func_names() {
    SurfaceMesh sm;
    p::list l;
    for (auto fn:sm.get_names())
        l.append(p::object(fn));
    return l;
}

BOOST_PYTHON_MODULE(AlgebraicSurfaces) {
    def("algebraic_surface", algebraic_surface, (p::arg("n_func")), (p::arg("resolution")));
    def("func_name", func_name, (p::arg("n_func")));
    def("func_names", func_names);
}
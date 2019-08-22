#include "AlgebraicSurfaces.h"
#include "Surfaces.h"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;

#include "python.h"

static p::list algebraic_surface(int n_func, int resolution, int color_map_index) {
    return mesh2list( SurfaceMesh().generate_mesh(n_func, resolution, color_map_index) );
}

static p::object func_name(int n_func) {
    return p::object( SurfaceMesh().get_name(n_func) );
}

static p::list func_names() {
    p::list l;
    for (auto fn:SurfaceMesh().get_names())
        l.append(p::object(fn));
    return l;
}



BOOST_PYTHON_MODULE(AlgebraicSurfaces) {
    def("algebraic_surface", algebraic_surface, (p::arg("n_func")), (p::arg("resolution")), (p::arg("color_map_index"))=1);
    def("func_name", func_name, (p::arg("n_func")));
    def("func_names", func_names);
}
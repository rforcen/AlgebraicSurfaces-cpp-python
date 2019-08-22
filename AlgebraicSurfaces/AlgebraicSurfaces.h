#ifndef AlgebraicSurfacesH
#define AlgebraicSurfacesH

// algebraic surfaces in parametric form
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string>
#include <vector>
#include <simd/simd.h>
#include "Thread.h"
#include "color_map.h"

using std::vector, std::string, std::tuple;

typedef vector<float>Coord;
typedef vector<int>Indexes;
typedef Coord Normal;
typedef Coord Texture;

typedef Coord Colors;
typedef Coord Textures;
typedef Coord Coords;
typedef Coord Normals;

typedef tuple<Indexes, Coords, Textures, Normals, Colors> Mesh;


class AlgebraicSurface {
public:
    int res=100, nCoords=0;

    float fromU, toU, fromV, toV, difU, difV;
    float diff;
    float _max, _min;

    float twoPi=M_PI * 2.0f;
    float fscale=1;
    int n_color_map=1;

    Mesh mesh;

    Indexes indexes;
    Coords coords;
    Textures textures;
    Normals normals;
    Colors colors;

    float max(float x, float y) {return (x > y) ? x : y;	}
    float sqr(float x)          {return x*x;	}
    float sqr3(float x)         {return x*x*x;	}
    float cube(float x)         {return x*x*x;	}
    float sqr4(float x)         {return x*x*x*x;	}
    float sqr5(float x)         {return x*x*x*x*x;	}

    AlgebraicSurface() {}

    /*
     Base Parametric Surface class
     required implementing virtual func: XYZ ParametricFunc(float u, float v)
     p.x=x(u,v)
     p.y=y(u,v)
     p.z=z(u,v)

     default func returns a 0 point
     */

    virtual simd_float3 ParametricFunc(float u, float v) = 0;

    void setResol(int resol) {
        if (resol < Thread::getnthreads()) resol=Thread::getnthreads()*2; // min resolution

        this->res = resol;	nCoords = res*res;
    }
    float rad2Deg(float rad) {	return rad * 180.0f / M_PI;	}
    float scaleU(float val)  {	return val*difU + fromU;	}
    float scaleV(float val)  {  return val*difV + fromV;	}

    simd_float3 evaluate(float u, float v) {
        return ParametricFunc(scaleU(u), scaleV(v));
    }

    void calcDif() {
        diff = fabs(_max - _min);
        fscale = (diff == 0) ? 1 : 1. / diff; // autoscale
    }

    simd_float3  calcNormal(simd_float3 p0, simd_float3 p1, simd_float3 p2) {
        auto cr = simd::cross(p1-p2, p1-p0);

        if (simd::length(cr)!=0) return simd::normalize( cr );
        else                     return simd_float3{1,1,1};
    }


    simd_float3  calcNormal(int i, int j, int ix_coord) { // ix_coord = i*res+j
        simd_float3 v[3];
        int ix=ix_coord*3;

        v[0] = simd_float3 { coords[ix], coords[ix+1], coords[ix+2] };

        ix += 3 * ( (j==res-1) ? -1 : +1 );
        v[1] = simd_float3 { coords[ix], coords[ix+1], coords[ix+2] };

        ix += 3*res * ( (i==res-1) ? -1 : +1 );
        v[2] = simd_float3 { coords[ix], coords[ix+1], coords[ix+2] };

        return calcNormal( v[0], v[1], v[2] );
    }

     void minMaxP(simd_float3 &p) { // update _min/_max in p
        _max = std::max(_max, simd::reduce_max(p));
        _min = std::min(_min, simd::reduce_min(p));
    }

    void addNormals() {
        for (int i=0; i<res; i++)
            for (int j=0; j<res; j++) {
                int ix_coord = i*res + j;

                auto n = calcNormal(i,j, ix_coord);

                for (int k=0; k<3; k++) normals[ix_coord*3+k] = n[k];
            }
    }

    void addIndexes() {
        int ic, i, ir;

        for (i=0, ir=res, ic=0; ic<nCoords && ir<nCoords; i++)
            indexes[i] = !(i&1) ? ic++ : ir++;
        indexes.resize(i);
    }

    void scaleCoords() {
        calcDif();
        for (auto &c:coords) c*=fscale;
    }

    void init_mesh_vectors() { //
        int n_coords = nCoords*3, n_txt=nCoords*2;

        indexes = Indexes(n_coords);
        coords  = Coords(n_coords);
        normals = Normals(n_coords);
        textures= Textures(n_txt);
        colors  = Colors(n_coords);

        _max = -FLT_MAX;
        _min = +FLT_MAX;
    }

    void set_limits(float _fromU, float _toU, float _fromV, float _toV) {
        fromU = _fromU;        toU = _toU;
        fromV = _fromV;        toV = _toV; // define limits
        difU = fabs(fromU - toU);     difV = fabs(fromV - toV);
    }

    void addTextVertex(int i, int j, int ix_coord) { // add textures coords and vertex

        float tx = (float)j/res, ty=(float)i/res;

        auto p = evaluate(ty, tx);

        textures[ix_coord * 2 + 0]  = tx;
        textures[ix_coord * 2 + 1]  = ty;

        for (int k=0; k<3; k++)
            coords[ix_coord * 3 + k] = p[k];

        auto cm=color_map(ty, 0, M_PI*2, n_color_map);

        colors[ix_coord*3+0] = cm.r;
        colors[ix_coord*3+1] = cm.g;
        colors[ix_coord*3+2] = cm.b;


        minMaxP(p);
    }


    // multithreaded generation
    Mesh calcMesh(int resol, float _fromU, float _toU, float _fromV, float _toV, int color_map_index) { // calc & load

        n_color_map = color_map_index;

        set_limits(_fromU,  _toU,  _fromV,	 _toV);
        setResol(resol);

        init_mesh_vectors();

        Thread(res).run([this](int i) {  // use all CPU threads to generate res x res triangle fan
            for (int j=0; j<res; j++) {

                int ix_coord = i*res + j;
                addTextVertex(i, j, ix_coord);

            }
        });
        
        scaleCoords();

        addNormals();
        addIndexes();

        return Mesh(indexes, coords, textures, normals, colors);
    }

};


#endif

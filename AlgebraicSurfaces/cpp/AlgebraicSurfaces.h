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

using std::vector, std::string, std::tuple;

typedef vector<float>Coord;
typedef Coord Normal;
typedef Coord Texture;

typedef vector<vector<float>> Textures;
typedef vector<vector<float>> Coords;
typedef Coords Normals;

typedef tuple<Coords, Textures, Normals> Mesh;


class AlgebraicSurface {
public:
    int res, nCoords;

    float fromU, toU, fromV, toV, difU, difV;
    float mx, my, mz, mix, miy, miz, MaxVal, MinVal, Dif;

    bool scaled;
    float twoPi;
    float fscale;

    Mesh mesh;
    Coords coords;
    Textures textures;
    Normals normals;

    float max(float x, float y) {return (x > y) ? x : y;	}
    float sqr(float x)          {return x*x;	}
    float sqr3(float x)         {return x*x*x;	}
    float cube(float x)         {return x*x*x;	}
    float sqr4(float x)         {return x*x*x*x;	}
    float sqr5(float x)         {return x*x*x*x*x;	}

    AlgebraicSurface() {
        twoPi = M_PI * 2.0f;
        res = 35;
        scaled = false;
        fscale = 1;
    }

    /*
     Base Parametric Surface class
     required implementing virtual func: XYZ ParametricFunc(float u, float v)
     p.x=x(u,v)
     p.y=y(u,v)
     p.z=z(u,v)

     default func returns a 0 point
     */

    virtual simd_float3 ParametricFunc(float u, float v) = 0;
    // assign 'p' with x,y,z

    void setResol(int resol) {	this->res = resol;	}
    float rad2Deg(float rad) {	return rad * 180.0f / M_PI;	}
    float scaleU(float val)  {	return val*difU + fromU;	}
    float scaleV(float val)  {  return val*difV + fromV;	}
    float scale01(float val) {	return val / Dif;	} // keep center
    float scale0to1(float val) {	return (val - MinVal) / Dif;	} // scale 0..1
    void scale(simd_float3&p) {
        p.x = scale01(p.x);
        p.y = scale01(p.y);
        p.z = scale01(p.z);
    }
    simd_float3 Eval(float u, float v) {
        return ParametricFunc(scaleU(u), scaleV(v));
    }

    void minMaxP(simd_float3&p) { // update min/max in p
        if (p.x > mx)
            mx = p.x;
        if (p.y > my)
            my = p.y;
        if (p.z > mz)
            mz = p.z;
        if (p.x < mix)
            mix = p.x;
        if (p.y < miy)
            miy = p.y;
        if (p.z < miz)
            miz = p.z;
    }

    void calcDif() {
        if (mx < my) {
            if (mx < mz)
                MaxVal = mx;
            else
                MaxVal = mz;
        }
        else {
            if (my < mz)
                MaxVal = my;
            else
                MaxVal = mz;
        }
        if (mix < miy) {
            if (mix < miz)
                MinVal = mix;
            else
                MinVal = miz;
        }
        else {
            if (miy < miz)
                MinVal = miy;
            else
                MinVal = miz;
        }

        Dif = fabs(MaxVal - MinVal);
        fscale = (Dif == 0) ? 1 : 1. / Dif; // autoscale
    }

    simd_float3  calcNormal(simd_float3 p0, simd_float3 p1, simd_float3 p2) {
        return simd::normalize( simd::cross(p1-p2, p1-p0) );
    }

    void addTextVertex(int ix_coord, float tx, float ty) { // add textures coords and vertex
        auto p = Eval(tx, ty);
        p *= fscale;

        textures[ix_coord]  = Coord{tx,ty};
        coords[ix_coord]    = Coord{p.x, p.y, p.z};

        minMaxP(p);
    }

    simd_float3 c2f3(Coord &c) { // convert vector<float> to simd_float3
        return simd_float3{c[0], c[1], c[2]};

    }
    void addNormals(int ix_coord) { // from last 3 coords
        auto p0=coords[ix_coord+0], p1=coords[ix_coord+1], p2=coords[ix_coord+3];

        auto n = calcNormal(c2f3(p0), c2f3(p1), c2f3(p2));
        for (int i=0; i<4; i++)
            normals[ix_coord + i] = Coord{n.x, n.y, n.z};
    }

    void scaleCoords() {
        for (auto &c:coords) {
            c[0]*=fscale;
            c[1]*=fscale;
            c[2]*=fscale;
        }
    }

    void initMinMax() {
        mx = my = mz = -FLT_MAX;
        mix = miy = miz = FLT_MAX;
    }

    void init_mesh_vectors() { //
        int n_coords = res*res*4;
        coords  = Coords(n_coords, Coord());
        normals = Normals(n_coords, Normal());
        textures= Textures(n_coords, Texture());
    }

    // multithreaded generation
    Mesh calcMesh(int resol, float _fromU, float _toU, float _fromV,	float _toV) { // calc & load

        fromU = _fromU;
        toU = _toU;
        fromV = _fromV;
        toV = _toV; // define limits
        difU = fabs(fromU - toU);
        difV = fabs(fromV - toV);

        setResol(resol);
        initMinMax();

        init_mesh_vectors();

        float dt = 1.0f / res;

        Thread(res).run([this, dt](int i) {  // use all CPU threads to generate res x res QUADS

            float idt = i * dt, jdt=0;

            for (int j=0; j<res; j++, jdt+=dt) {

                int ix_coord=(i * res+j) * 4;

                addTextVertex(ix_coord+0, idt,          jdt);
                addTextVertex(ix_coord+1, idt + dt,     jdt);
                addTextVertex(ix_coord+2, idt + dt,     jdt + dt);
                addTextVertex(ix_coord+3, idt,          jdt + dt);

                addNormals(ix_coord); // 4 x last 3 coords

            }
        });

        calcDif();
        fscale = (Dif == 0) ? 1 : 1. / Dif; // autoscale
        scaleCoords();

        return Mesh(coords, textures, normals);
    }
};

// surface implementation

class Tanaka : public AlgebraicSurface {
public:
    float paramSets[4][7] = {
        {0, 4, 3, 4, 5, 7, 4}, {0, 4, 3, 0, 5, 7, 4}, {0, 3, 4, 8, 5, 5, 2},
        {14, 3, 1, 8, 5, 5, 2}};
    float a, b1, b2, c, d, w, h;

    Tanaka() {
        setParam(0);
        a = 0, // center hole size of a torus
                b1 = 4, // number of cross
                b2 = 3, // number of cross
                c = 4, // distance from the center of rotation
                d = 5, // number of torus
                w = 7, // gap width
                h = 4; // height

    }

    void setParam(int np) {
        np=np%4;
        a = paramSets[np][0];        b1 = paramSets[np][1];        b2 = paramSets[np][2];
        c = paramSets[np][3];        d = paramSets[np][4];        w = paramSets[np][5];
    }
    float getNTorus() {      return d;    } // number of torus
    float f(float v) {  return sin(2*sin(sin(sin(v))));    }

    simd_float3 ParametricFunc(float s, float t) {
        simd_float3 p;
        p.x = (a - cos(t) + w * sin(b1 * s)) * cos(b2 * s);
        p.y = (a - cos(t) + w * sin(b1 * s)) * f(b2 * s);
        p.z = h * (w * sin(b1 * s) + f(t)) + c;
        return p;
    }
};

class Kuen : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = 2 * cosh(v) * (cos(u) + u * sin(u)) / (cosh(v) * cosh(v) + u * u);
        p.y = 2 * cosh(v) * (-u * cos(u) + sin(u)) /
                (cosh(v) * cosh(v) + u * u);
        p.z = v - (2 * sinh(v) * cosh(v)) / (cosh(v) * cosh(v) + u * u);
        return p;
    }
};

class ButterFly : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float t1 = (exp(cos(u)) - 2 * cos(4 * u) + sqr5(sin(u / 12))) * sin(v);
        p.x = sin(u) * t1;
        p.y = cos(u) * t1;
        p.z = sin(v);
        return p;
    }
};

class RoseSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float a = 1, n = 7;
        p.x = a * sin(n * u) * cos(u) * sin(v);
        p.y = a * sin(n * u) * sin(u) * sin(v);
        p.z = cos(v) / (n * 3);
        return p;
    }
};

class UpDownShellSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        return simd_float3 {
         u * sin(u) * cos(v),
         u * cos(u) * cos(v),
         u * sin(v) }; // -10,10, -10,10
    }
};

class CayleySurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        // plot3d([u*sin(v)-u*cos(v), u^2*sin(v)*cos(v), u^3*sin(v)^2*cos(v)], u=0..0.5, v=0..2*M_PI, numpoints=1000);
        p.x = u * sin(v) - u * cos(v);
        p.y = sqr(u) * sin(v) * cos(v);
        p.z = cube(u) * sqr(sin(v)) * cos(v);
        return p;
    }
};

class PluckerConoidSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u * v;
        p.y = u * sqrt(1 - sqr(v));
        p.z = 1 - sqr(v);
        return p;
    }
};

class AmmoniteSurface : public AlgebraicSurface {
    inline float W(float u) { return pow(u / (2*M_PI), 2.2f);	}

    simd_float3 ParametricFunc(float u, float v) {
        float N = 5.6f; // number of turns
        float F = 120.0f; // wave frequency
        float A = 0.2f; // wave amplitude
        simd_float3 p;
        p.x = W(u) * cos(N * u) * (2 + sin(v + cos(F * u) * A));
        p.y = W(u) * sin(N * u) * (2 + sin(v + cos(F * u) * A));
        p.z = W(u) * cos(v);
        return p;
    }
};

class AppleSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        float R1 = 4, R2 = 3.8f;
        simd_float3 p;

        p.x = cos(u) * (R1 + R2 * cos(v)) + pow((v / M_PI), 100);
        p.y = sin(u) * (R1 + R2 * cos(v)) + 0.25f * cos(5 * u);
        p.z = -2.3f * log(1 - v * 0.3157f) + 6 * sin(v) + 2 * cos(v);
        return p;
    }
};

class AstroidalEllipseSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        float A = 1, B = 1, C = 1;
        simd_float3 p;

        p.x = pow(A * cos(u) * cos(v), 3);
        p.y = pow(B * sin(u) * cos(v), 3);
        p.z = pow(C * sin(v), 3);
        return p;
    }
    // <0,0>,<2*pi,2*pi>
};

class BohemianDomeSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float A = 0.5f, B = 1.5f, C = 1;
        p.x = A * cos(u);
        p.y = B * cos(v) + A * sin(u);
        p.z = C * sin(v);
        return p;
    }
};

class ConicalSpiralSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u * v * sin(15 * v);
        p.y = v;
        p.z = u * v * cos(15 * v);
        return p;
        // <0,-1>,<1,1>
    }
};

class EnneperSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u - u * u * u / 3 + u * v * v;
        p.y = v - v * v * v / 3 + v * u * u;
        p.z = u * u - v * v;
        return p;
    }
};

class Scherk : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float aa = 0.1f;
        v += 0.1;
        p.x = u;
        p.y = v;
        p.z = (log(fabs(cos(aa * v) / cos(aa * u)))) / aa;
        return p;
    }
};

class Dini : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float psi = 0.3f; // aa;
        if (psi < 0.001f)
            psi = 0.001f;
        if (psi > 0.999f)
            psi = 0.999f;
        psi = psi * M_PI;
        float sinpsi = sin(psi);
        float cospsi = cos(psi);
        float g = (u - cospsi * v) / sinpsi;
        float s = exp(g);
        float r = (2 * sinpsi) / (s + 1 / s);
        float t = r * (s - 1 / s) * 0.5f;

        p.x = (u - t);
        p.y = (r * cos(v));
        p.z = (r * sin(v));

        return p;
    }
};

class KleinBottle : public AlgebraicSurface {
public:
    float t;

    KleinBottle() {	t = 4.5f;	}

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float tmp = (4 + 2 * cos(u) * cos(t * v) - sin(2 * u) * sin(t * v));
        p.x = sin(v) * tmp;
        p.y = cos(v) * tmp;
        p.z = 2 * cos(u) * sin(t * v) + sin(2 * u) * cos(t * v);
        return p;
    }
};

class KleinBottle0 : public AlgebraicSurface {
public:
    float t;

    KleinBottle0() { t = 4.5f;	}

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = (0 <= u && u < M_PI) ? 6 * cos(u) * (1 + sin(u)) +
                                     4 * (1 - 0.5f * cos(u)) * cos(u) * cos(v) :
                                     6 * cos(u) * (1 + sin(u)) + 4 * (1 - 0.5f * cos(u)) * cos(v + M_PI);
        p.y = (0 <= u && u < M_PI) ? 16 * sin(u) + 4 * (1 - 0.5f * cos(u)) * sin
                                     (u) * cos(v) : 16 * sin(u);
        p.z = 4 * (1 - 0.5f * cos(u)) * sin(v);
        return p;
    }
};

class Bour : public AlgebraicSurface {
public:
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u * cos(v) - 0.5f * u * u * cos(2 * v);
        p.y = -u * sin(v) - 0.5f * u * u * sin(2 * v);
        p.z = 4 / 3 * pow(u, 1.5f) * cos(1.5f * v);
        return p;
    }
};

class BreatherSurface : public AlgebraicSurface {
public:
    float aa, w1, w;

    BreatherSurface() {
        aa = 0.45f; // Values from 0.4 to 0.6 produce sensible results
        w1 = 1 - aa * aa;
        w = sqrt(w1);
    }

    float d(float u, float v) {
        return aa * (pow((w * cosh(aa * u)), 2) + pow((aa * sin(w * v)), 2));
    }

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = -u + (2 * w1 * cosh(aa * u) * sinh(aa * u) / d(u, v));
        p.y = 2 * w * cosh(aa * u) * (-(w * cos(v) * cos(w * v)) -
                                      (sin(v) * sin(w * v))) / d(u, v);
        p.z = 2 * w * cosh(aa * u) * (-(w * sin(v) * cos(w * v)) +
                                      (cos(v) * sin(w * v))) / d(u, v);
        return p;
    }

};

/*
 Seashell surfaces
 These are parametric isosurfaces using variations of functions like:

 #declare N=5.6;  // number of turns
 #declare H=3.5;  // height
 #declare P=2;    // power
 #declare L=4;    // Controls spike length
 #declare K=9;    // Controls spike sharpness
 */

class Seashell : public AlgebraicSurface {
public:
    float N, H, P, L, K;

    Seashell() {
        N = 5.6f; // number of turns
        H = 3.5f; // height
        P = 2; // power
        L = 4; // Controls spike length
        K = 9;
    }

    inline float W(float u) {	return pow(u / (2*M_PI), P);	}

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = W(u) * cos(N * u) * (1 + cos(v));
        p.y = W(u) * sin(N * u) * (1 + cos(v));
        p.z = W(u) * (sin(v) + pow(sin(v / 2), K) * L) +
                H * pow(u / (2 * M_PI), P + 1);
        return p;
    }
};


class TudorRoseSurface : public AlgebraicSurface {

    float R(float u, float v) {
        return cos(v) * cos(v) * max(fabs(sin(4 * u)),
                                     0.9f - 0.2f * fabs(cos(8 * u)));
    }

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = R(u, v) * cos(u) * cos(v);
        p.y = R(u, v) * sin(u) * cos(v);
        p.z = R(u, v) * sin(v) * 0.5f;
        return p;
    }
};

// Cap's surface
class CapSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = 0.5f * cos(u) * sin(2 * v);
        p.y = 0.5f * sin(u) * sin(2 * v);
        p.z = 0.5f * (sqr(cos(v)) - sqr(cos(u)) * sqr(sin(v)));
        return p;
    }
};

// Boy's surface
class BoySurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float dv = (2 - sqrt(2) * sin(3 * u) * sin(2 * v)), d1 =
                (cos(u) * sin(2 * v)), d2 = sqrt(2) * sqr(cos(v));

        p.x = (d2 * cos(2 * u) + d1) / dv;
        p.y = (d2 * sin(2 * u) + d1) / dv;
        p.z = (3 * sqr(cos(v))) / (2 - sqrt(2) * sin(3 * u) * sin(2 * v));

        /* x,y,z deformations
         float a=2.1,b=1.3,c=1;
         p.x*=a; p.y*=b; p.z*=c;
         */
        return p;
    }
};


class RomanSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float r, float t) // 0 <= r <= 1   |  0 <= t <= 2*PI
    {
        simd_float3 p;
        float r2 = r * r, rq = sqrt((1 - r2)), st = sin(t), ct = cos(t);
        p.x = r2 * st * ct;
        p.y = r * st * rq;
        p.z = r * ct * rq;
        return p;
    }
};

class SurfaceMesh { // TODO: add here new funcs

    RoseSurface RoseSrf;    RomanSurface Rsrf;    CapSurface Csrf;
    BoySurface Bsrf;        Seashell Ssrf;        TudorRoseSurface Tsrf;
    BreatherSurface Brsrf;  KleinBottle KBsrf;    KleinBottle0 KBsrf0;
    Bour bourS;             Dini biniS;           EnneperSurface EnneperS;
    Scherk scherk;          ConicalSpiralSurface ConSpi;
    BohemianDomeSurface BDsurf;
    AstroidalEllipseSurface AESurf;
    AppleSurface SApple;    AmmoniteSurface AmmSrf; PluckerConoidSurface PCSrf;
    CayleySurface CYSrf;    UpDownShellSurface UDSSrf;
    ButterFly bfy;          Kuen KuenSrf;       Tanaka tanaka;


    float twoPi;

public:
    static const int n_funcs=27;
    vector<string>srfNames = {
            "cap","boy", "roman", "sea shell", "tudor rose", "breather",
            "klein bottle", "klein bottle 0", "bour", "dini", "enneper",
            "scherk", "conical spiral", "bohemian dome", "astrodial ellipse",
            "apple", "ammonite", "plucker comoid", "cayley", "up down shell",
            "butterfly", "rose", "kuen", "tanaka-0", "tanaka-1", "tanaka-2",
            "tanaka-3"};

    SurfaceMesh() : twoPi(M_PI * 2) {}

    Mesh generate_mesh(int ns=0, int resol=100) {

        Mesh mesh;
        switch (ns % n_funcs) {
            case 0:	mesh=Csrf.calcMesh(resol, 0, M_PI, 0, M_PI);		break;
            case 1:	mesh=Bsrf.calcMesh(resol, 0, M_PI, 0, M_PI);		break;
            case 2:	mesh=Rsrf.calcMesh(resol, 0, 1, 0, twoPi);			break;
            case 3:	mesh=Ssrf.calcMesh(resol, 0, twoPi, 0, twoPi);		break;
            case 4:	mesh=Tsrf.calcMesh(resol, 0, M_PI, 0, M_PI);		break;
            case 5:	mesh=Brsrf.calcMesh(resol, -20, 20, 20, 80);		break;
            case 6:	mesh=KBsrf.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 7:	mesh=KBsrf0.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 8:	mesh=bourS.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 9:	mesh=biniS.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 10:mesh=EnneperS.calcMesh(resol, -1, 1, -1, 1);		break;
            case 11:mesh=scherk.calcMesh(resol, 1, 30, 1, 30);			break;
            case 12:mesh=ConSpi.calcMesh(resol, 0, 1, -1, 1);			break;
            case 13:mesh=BDsurf.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 14:mesh=AESurf.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 15:mesh=SApple.calcMesh(resol, 0, twoPi, -M_PI, M_PI);break;
            case 16:mesh=AmmSrf.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 17:mesh=PCSrf.calcMesh(resol, -2, 2, -1, 1);			break;
            case 18:mesh=CYSrf.calcMesh(resol, 0, 3, 0, twoPi);		break;
            case 19:mesh=UDSSrf.calcMesh(resol, -10, 10, -10, 10);		break;
            case 20:mesh=bfy.calcMesh(resol, 0, twoPi, 0, twoPi);		break;
            case 21:mesh=RoseSrf.calcMesh(resol, 0, twoPi, 0, twoPi);	break;
            case 22:mesh=KuenSrf.calcMesh(resol, -4, 4, -3.75f, +3.75f);break;
            case 23:
            case 24:
            case 25:
            case 26:tanaka.setParam(ns - 23);
                    mesh=tanaka.calcMesh(resol, 0, twoPi, 0, twoPi);
                    break;
        }

        return mesh;
    }

    string get_name(int nf)         { return srfNames[nf % n_funcs];   }
    vector<string> get_names()      { return srfNames; }
};

#endif

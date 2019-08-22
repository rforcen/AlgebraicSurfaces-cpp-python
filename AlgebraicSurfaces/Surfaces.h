// samples surfaces implementation

#include "AlgebraicSurfaces.h"

class Tanaka : public AlgebraicSurface {
public:
    float paramSets[4][7] = {
        {0, 4, 3, 4, 5, 7, 4}, {0, 4, 3, 0, 5, 7, 4}, {0, 3, 4, 8, 5, 5, 2},
        {14, 3, 1, 8, 5, 5, 2}};
    float a, b1, b2, c, d, w, h;

    Tanaka() {
        setParam(0);
        a = 0; // center hole size of a torus
                b1 = 4; // number of cross
                b2 = 3; // number of cross
                c = 4; // distance from the center of rotation
                d = 5; // number of torus
                w = 7; // gap width
                h = 4; // height

    }

    void setParam(int np) {
        np=np%4;
        a = paramSets[np][0];        b1 = paramSets[np][1];        b2 = paramSets[np][2];
        c = paramSets[np][3];        d = paramSets[np][4];        w = paramSets[np][5];
    }
    float getNTorus() {      return d;    } // number of torus
    float f(float v) {  return sinf(2*sinf(sinf(sinf(v))));    }

    simd_float3 ParametricFunc(float s, float t) {
        return simd_float3 {
         (a - cosf(t) + w * sinf(b1 * s)) * cosf(b2 * s),
         (a - cosf(t) + w * sinf(b1 * s)) * f(b2 * s),
         h * (w * sinf(b1 * s) + f(t)) + c};
    }
};

class Kuen : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = 2 * cosh(v) * (cosf(u) + u * sinf(u)) / (cosh(v) * cosh(v) + u * u);
        p.y = 2 * cosh(v) * (-u * cosf(u) + sinf(u)) /
                (cosh(v) * cosh(v) + u * u);
        p.z = v - (2 * sinh(v) * cosh(v)) / (cosh(v) * cosh(v) + u * u);
        return p;
    }
};

class ButterFly : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float t1 = (exp(cosf(u)) - 2 * cosf(4 * u) + sqr5(sinf(u / 12))) * sinf(v);
        p.x = sinf(u) * t1;
        p.y = cosf(u) * t1;
        p.z = sinf(v);
        return p;
    }
};

class RoseSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float a = 1, n = 7;
        p.x = a * sinf(n * u) * cosf(u) * sinf(v);
        p.y = a * sinf(n * u) * sinf(u) * sinf(v);
        p.z = cosf(v) / (n * 3);
        return p;
    }
};

class UpDownShellSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        return simd_float3 {
         u * sinf(u) * cosf(v),
         u * cosf(u) * cosf(v),
         u * sinf(v) }; // -10,10, -10,10
    }
};

class CayleySurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        // plot3d([u*sinf(v)-u*cosf(v), u^2*sinf(v)*cosf(v), u^3*sinf(v)^2*cosf(v)], u=0..0.5, v=0..2*M_PI, numpoints=1000);
        p.x = u * sinf(v) - u * cosf(v);
        p.y = sqr(u) * sinf(v) * cosf(v);
        p.z = cube(u) * sqr(sinf(v)) * cosf(v);
        return p;
    }
};

class PluckerConoidSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u * v;
        p.y = u * sqrtf(1 - sqr(v));
        p.z = 1 - sqr(v);
        return p;
    }
};

class AmmoniteSurface : public AlgebraicSurface {
    inline float W(float u) { return powf(u / (2*M_PI), 2.2f);	}

    simd_float3 ParametricFunc(float u, float v) {
        float N = 5.6f; // number of turns
        float F = 120.0f; // wave frequency
        float A = 0.2f; // wave amplitude
        simd_float3 p;
        p.x = W(u) * cosf(N * u) * (2 + sinf(v + cosf(F * u) * A));
        p.y = W(u) * sinf(N * u) * (2 + sinf(v + cosf(F * u) * A));
        p.z = W(u) * cosf(v);
        return p;
    }
};

class AppleSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        float R1 = 4, R2 = 3.8f;
        simd_float3 p;

        p.x = cosf(u) * (R1 + R2 * cosf(v)) + powf((v / M_PI), 100);
        p.y = sinf(u) * (R1 + R2 * cosf(v)) + 0.25f * cosf(5 * u);
        p.z = -2.3f * logf(1 - v * 0.3157f) + 6 * sinf(v) + 2 * cosf(v);
        return p;
    }
};

class AstroidalEllipseSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        float A = 1, B = 1, C = 1;
        simd_float3 p;

        p.x = powf(A * cosf(u) * cosf(v), 3);
        p.y = powf(B * sinf(u) * cosf(v), 3);
        p.z = powf(C * sinf(v), 3);
        return p;
    }
    // <0,0>,<2*pi,2*pi>
};

class BohemianDomeSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float A = 0.5f, B = 1.5f, C = 1;
        p.x = A * cosf(u);
        p.y = B * cosf(v) + A * sinf(u);
        p.z = C * sinf(v);
        return p;
    }
};

class ConicalSpiralSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u * v * sinf(15 * v);
        p.y = v;
        p.z = u * v * cosf(15 * v);
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
        p.z = (logf(fabs(cosf(aa * v) / cosf(aa * u)))) / aa;
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
        float sinpsi = sinf(psi);
        float cospsi = cosf(psi);
        float g = (u - cospsi * v) / sinpsi;
        float s = exp(g);
        float r = (2 * sinpsi) / (s + 1 / s);
        float t = r * (s - 1 / s) * 0.5f;

        p.x = (u - t);
        p.y = (r * cosf(v));
        p.z = (r * sinf(v));

        return p;
    }
};

class KleinBottle : public AlgebraicSurface {
public:
    float t;

    KleinBottle() {	t = 4.5f;	}

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float tmp = (4 + 2 * cosf(u) * cosf(t * v) - sinf(2 * u) * sinf(t * v));
        p.x = sinf(v) * tmp;
        p.y = cosf(v) * tmp;
        p.z = 2 * cosf(u) * sinf(t * v) + sinf(2 * u) * cosf(t * v);
        return p;
    }
};

class KleinBottle0 : public AlgebraicSurface {
public:
    float t;

    KleinBottle0() { t = 4.5f;	}

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = (0 <= u && u < M_PI) ? 6 * cosf(u) * (1 + sinf(u)) +
                                     4 * (1 - 0.5f * cosf(u)) * cosf(u) * cosf(v) :
                                     6 * cosf(u) * (1 + sinf(u)) + 4 * (1 - 0.5f * cosf(u)) * cosf(v + M_PI);
        p.y = (0 <= u && u < M_PI) ? 16 * sinf(u) + 4 * (1 - 0.5f * cosf(u)) * sinf
                                     (u) * cosf(v) : 16 * sinf(u);
        p.z = 4 * (1 - 0.5f * cosf(u)) * sinf(v);
        return p;
    }
};

class Bour : public AlgebraicSurface {
public:
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = u * cosf(v) - 0.5f * u * u * cosf(2 * v);
        p.y = -u * sinf(v) - 0.5f * u * u * sinf(2 * v);
        p.z = 4 / 3 * powf(u, 1.5f) * cosf(1.5f * v);
        return p;
    }
};

class BreatherSurface : public AlgebraicSurface {
public:
    float aa, w1, w;

    BreatherSurface() {
        aa = 0.45f; // Values from 0.4 to 0.6 produce sensible results
        w1 = 1 - aa * aa;
        w = sqrtf(w1);
    }

    float d(float u, float v) {
        return aa * (powf((w * cosh(aa * u)), 2) + powf((aa * sinf(w * v)), 2));
    }

    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = -u + (2 * w1 * cosh(aa * u) * sinh(aa * u) / d(u, v));
        p.y = 2 * w * cosh(aa * u) * (-(w * cosf(v) * cosf(w * v)) -
                                      (sinf(v) * sinf(w * v))) / d(u, v);
        p.z = 2 * w * cosh(aa * u) * (-(w * sinf(v) * cosf(w * v)) +
                                      (cosf(v) * sinf(w * v))) / d(u, v);
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

    inline float W(float u) {	return powf(u / (2*M_PI), P);	}

    simd_float3 ParametricFunc(float u, float v) {
        return simd_float3 {
         W(u) * cosf(N * u) * (1 + cosf(v)),
         W(u) * sinf(N * u) * (1 + cosf(v)),
         W(u) * (sinf(v) + powf(sinf(v / 2), K) * L) +  H * powf(u / (2 * M_PI), P + 1) };
    }
};


class TudorRoseSurface : public AlgebraicSurface {

    float R(float u, float v) {
        return cosf(v) * cosf(v) * max(fabs(sinf(4 * u)),
                                     0.9f - 0.2f * fabs(cosf(8 * u)));
    }

    simd_float3 ParametricFunc(float u, float v) {
        return simd_float3 {
         R(u, v) * cosf(u) * cosf(v),
         R(u, v) * sinf(u) * cosf(v),
         R(u, v) * sinf(v) * 0.5f };
    }
};

// Cap's surface
class CapSurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        p.x = 0.5f * cosf(u) * sinf(2 * v);
        p.y = 0.5f * sinf(u) * sinf(2 * v);
        p.z = 0.5f * (sqr(cosf(v)) - sqr(cosf(u)) * sqr(sinf(v)));
        return p;
    }
};

// Boy's surface
class BoySurface : public AlgebraicSurface {
    simd_float3 ParametricFunc(float u, float v) {
        simd_float3 p;
        float dv = (2 - sqrtf(2) * sinf(3 * u) * sinf(2 * v)), d1 =
                (cosf(u) * sinf(2 * v)), d2 = sqrtf(2) * sqr(cosf(v));

        p.x = (d2 * cosf(2 * u) + d1) / dv;
        p.y = (d2 * sinf(2 * u) + d1) / dv;
        p.z = (3 * sqr(cosf(v))) / (2 - sqrtf(2) * sinf(3 * u) * sinf(2 * v));

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
        float r2 = r * r, rq = sqrtf((1 - r2)), st = sinf(t), ct = cosf(t);
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

    Mesh generate_mesh(int ns=0, int resol=100, int color_map_index=1) {

        Mesh mesh;
        switch (ns % n_funcs) {
            case 0:	mesh=Csrf.calcMesh(resol, 0, M_PI, 0, M_PI, color_map_index);		break;
            case 1:	mesh=Bsrf.calcMesh(resol, 0, M_PI, 0, M_PI, color_map_index);		break;
            case 2:	mesh=Rsrf.calcMesh(resol, 0, 1, 0, twoPi, color_map_index);			break;
            case 3:	mesh=Ssrf.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);		break;
            case 4:	mesh=Tsrf.calcMesh(resol, 0, M_PI, 0, M_PI, color_map_index);		break;
            case 5:	mesh=Brsrf.calcMesh(resol, -20, 20, 20, 80, color_map_index);		break;
            case 6:	mesh=KBsrf.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 7:	mesh=KBsrf0.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 8:	mesh=bourS.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 9:	mesh=biniS.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 10:mesh=EnneperS.calcMesh(resol, -1, 1, -1, 1, color_map_index);		break;
            case 11:mesh=scherk.calcMesh(resol, 1, 30, 1, 30, color_map_index);			break;
            case 12:mesh=ConSpi.calcMesh(resol, 0, 1, -1, 1, color_map_index);			break;
            case 13:mesh=BDsurf.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 14:mesh=AESurf.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 15:mesh=SApple.calcMesh(resol, 0, twoPi, -M_PI, M_PI, color_map_index);break;
            case 16:mesh=AmmSrf.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 17:mesh=PCSrf.calcMesh(resol, -2, 2, -1, 1, color_map_index);			break;
            case 18:mesh=CYSrf.calcMesh(resol, 0, 3, 0, twoPi, color_map_index);		break;
            case 19:mesh=UDSSrf.calcMesh(resol, -10, 10, -10, 10, color_map_index);		break;
            case 20:mesh=bfy.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);		break;
            case 21:mesh=RoseSrf.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);	break;
            case 22:mesh=KuenSrf.calcMesh(resol, -4, 4, -3.75f, +3.75f, color_map_index);break;
            case 23:
            case 24:
            case 25:
            case 26:tanaka.setParam(ns - 23);
                    mesh=tanaka.calcMesh(resol, 0, twoPi, 0, twoPi, color_map_index);
                    break;
        }

        return mesh;
    }

    string get_name(int nf)         { return srfNames[nf % n_funcs];   }
    vector<string> get_names()      { return srfNames; }
};

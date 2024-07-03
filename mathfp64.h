#pragma once

#if defined(_WINDOWS)
#include <cmath>
#endif

// https://github.com/MikeLankamp/fpm
//#include "fixed.h"
//#include "math.h"

// https://sourceforge.net/projects/fixedptc/
#include "fixedptc.h"

 typedef fixed1616 fixed_16_16;
 //typedef float fixed_16_16;

struct vec2fp
{
    fixed_16_16 x, y;
    vec2fp() : x(0), y(0) {}
    vec2fp(fixed_16_16 inx, fixed_16_16 iny)
    {
        x = inx;
        y = iny;
    }
    vec2fp operator+(const vec2fp& a) const
    {
        vec2fp out;
        out.x = x + a.x;
        out.y = y + a.y;
        return out;
    }
    vec2fp& operator+=(const vec2fp& a)
    {
        this->x += a.x;
        this->y += a.y;
        return *this;
    }
    vec2fp& operator-=(const vec2fp& a)
    {
        this->x -= a.x;
        this->y -= a.y;
        return *this;
    }
    vec2fp& operator*=(const vec2fp& a)
    {
        this->x *= a.x;
        this->y *= a.y;
        return *this;
    }
    vec2fp& operator/=(const vec2fp& a)
    {
        this->x /= a.x;
        this->y /= a.y;
        return *this;
    }
    vec2fp& operator*=(const fixed_16_16& a)
    {
        this->x *= a;
        this->y *= a;
        return *this;
    }
    vec2fp& operator/=(const fixed_16_16& a)
    {
        this->x /= a;
        this->y /= a;
        return *this;
    }
    vec2fp operator-(const vec2fp& a) const
    {
        vec2fp out;
        out.x = x - a.x;
        out.y = y - a.y;
        return out;
    }
    vec2fp operator*(const vec2fp& a) const
    {
        vec2fp out;
        out.x = x * a.x;
        out.y = y * a.y;
        return out;
    }
    vec2fp operator*(const fixed_16_16& a) const
    {
        vec2fp out;
        out.x = x * a;
        out.y = y * a;
        return out;
    }
    vec2fp operator/(const vec2fp& a) const
    {
        vec2fp out;
        out.x = x / a.x;
        out.y = y / a.y;
        return out;
    }
    vec2fp operator/(const fixed_16_16& a) const
    {
        vec2fp out;
        out.x = x / a;
        out.y = y / a;
        return out;
    }
    fixed_16_16 dot(const vec2fp& a) const
    {
        return x * a.x + y * a.y;
    }
    fixed_16_16 mag() const
    {
        return sqrt(dot(*this));
    }
    vec2fp getnormalized() const
    {
        return *this / mag();
    }
};

struct vec3fp
{
    fixed_16_16 x = (fixed_16_16)0, y = (fixed_16_16)0, z = (fixed_16_16)0;
    vec3fp() {}
    vec3fp(fixed_16_16 inx, fixed_16_16 iny, fixed_16_16 inz)
    {
        x = inx;
        y = iny;
        z = inz;
    }
    vec3fp& operator=(vec3fp other)
    {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }
    vec3fp operator+(const vec3fp& a) const
    {
        vec3fp out;
        out.x = x + a.x;
        out.y = y + a.y;
        out.z = z + a.z;
        return out;
    }
    vec3fp operator-(const vec3fp& a) const
    {
        vec3fp out;
        out.x = x - a.x;
        out.y = y - a.y;
        out.z = z - a.z;
        return out;
    }
    vec3fp operator-() const
    {
        vec3fp out(x, y, z);
        out = out * (fixed_16_16 )-1.0f;
        return out;
    }
    vec3fp operator*(const vec3fp& a) const
    {
        vec3fp out;
        out.x = x * a.x;
        out.y = y * a.y;
        out.z = z * a.z;
        return out;
    }
    vec3fp operator*(const fixed_16_16& a) const
    {
        vec3fp out;
        out.x = x * a;
        out.y = y * a;
        out.z = z * a;
        return out;
    }
    vec3fp operator/(const vec3fp& a) const
    {
        vec3fp out;
        out.x = x / a.x;
        out.y = y / a.y;
        out.z = z / a.z;
        return out;
    }
    vec3fp operator/(const fixed_16_16& a) const
    {
        vec3fp out;
        out.x = x / a;
        out.y = y / a;
        out.z = z / a;
        return out;
    }
    vec3fp cross(const vec3fp& a) const
    {
        vec3fp out;
        out.x = y * a.z - z * a.y;
        out.y = z * a.x - x * a.z;
        out.z = x * a.y - y * a.x;
        return out;
    }
    fixed_16_16 mag() const
    {
        return sqrt(dot(*this));
    }
    fixed_16_16 dot(const vec3fp& a) const
    {
        return x * a.x + y * a.y + z * a.z;
    }
    vec3fp getnormalized() const
    {
        fixed_16_16 vectormag = mag();
        if (vectormag == (fixed_16_16)0) {
            return vec3fp((fixed_16_16)0, (fixed_16_16)0, (fixed_16_16)0);
        }
        else {
            return *this / vectormag;
        }
    }
    vec2fp xy() const
    {
        return vec2fp(x, y);
    }
};

struct vec4fp
{
    fixed_16_16 x, y, z, w;
    vec4fp() {}
    vec4fp(fixed_16_16 inx, fixed_16_16 iny, fixed_16_16 inz, fixed_16_16 inw)
    {
        x = inx;
        y = iny;
        z = inz;
        w = inw;
    }
    vec4fp(vec3fp invec, fixed_16_16 inw)
    {
        x = invec.x;
        y = invec.y;
        z = invec.z;
        w = inw;
    }
    vec4fp& operator=(vec4fp a)
    {
        x = a.x;
        y = a.y;
        z = a.z;
        w = a.w;
        return *this;
    }
    vec4fp operator*(const vec4fp& a) const
    {
        vec4fp out;
        out.x = x * a.x;
        out.y = y * a.y;
        out.z = z * a.z;
        out.w = w * a.w;
        return out;
    }
    vec3fp xyz() const
    {
        return vec3fp(x, y, z);
    }
    fixed_16_16 dot(const vec4fp& a) const
    {
        return x * a.x + y * a.y + z * a.z + w * a.w;
    }
};

struct mat44fp
{
    vec4fp rows[4];
    mat44fp() {}
    mat44fp(vec4fp r0, vec4fp r1, vec4fp r2, vec4fp r3)
    {
        rows[0] = r0;
        rows[1] = r1;
        rows[2] = r2;
        rows[3] = r3;
    }
    mat44fp& operator=(mat44fp other)
    {
        rows[0] = other.rows[0];
        rows[1] = other.rows[1];
        rows[2] = other.rows[2];
        rows[3] = other.rows[3];
        return *this;
    }
    inline
        mat44fp operator*(const mat44fp& mat)
    {
        mat44fp matout;
        vec4fp colvec0 = vec4fp(mat.rows[0].x, mat.rows[1].x, mat.rows[2].x, mat.rows[3].x);
        vec4fp colvec1 = vec4fp(mat.rows[0].y, mat.rows[1].y, mat.rows[2].y, mat.rows[3].y);
        vec4fp colvec2 = vec4fp(mat.rows[0].z, mat.rows[1].z, mat.rows[2].z, mat.rows[3].z);
        vec4fp colvec3 = vec4fp(mat.rows[0].w, mat.rows[1].w, mat.rows[2].w, mat.rows[3].w);
        for (int row = 0; row < 4; ++row)
        {
            matout.rows[row].x = rows[row].dot(colvec0);
            matout.rows[row].y = rows[row].dot(colvec1);
            matout.rows[row].z = rows[row].dot(colvec2);
            matout.rows[row].w = rows[row].dot(colvec3);
        }
        return matout;
    }
    inline
        vec4fp operator*(const vec4fp& vec)
    {
        vec4fp vecout;
        vecout.x = rows[0].dot(vec);
        vecout.y = rows[1].dot(vec);
        vecout.z = rows[2].dot(vec);
        vecout.w = rows[3].dot(vec);
        return vecout;
    }
};

template<class T>
void swapfp(T& a, T& b) {
    T temp = a;
    a = b;
    b = temp;
}

struct TriangleFP {
    vec3fp a, b, c;
    TriangleFP() {}
    TriangleFP(vec3fp ain, vec3fp bin, vec3fp cin) {
        a = ain; b = bin; c = cin;
    }
    vec3fp facenormal()
    {
        vec3fp cross = (a - b).cross(c - a);
        return cross.getnormalized();
    }
};

struct TriangleUVFP : public TriangleFP {
    vec2fp auv, buv, cuv;
    fixed_16_16 aw, bw, cw;
    TriangleUVFP() {}
    TriangleUVFP(vec3fp ain, vec3fp bin, vec3fp cin,
        vec2fp auvin, vec2fp buvin, vec2fp cuvin,
        fixed_16_16 awin, fixed_16_16 bwin, fixed_16_16 cwin) {
        a = ain; b = bin; c = cin;
        auv = auvin; buv = buvin; cuv = cuvin;
        aw = awin; bw = bwin; cw = cwin;
    }
    void scale(fixed_16_16 scale) {
        vec3fp center = ( a + b + c ) / (fixed_16_16)3.0f;
        a.x = (a.x - center.x) * scale + center.x;
        a.y = (a.y - center.y) * scale + center.y;
        b.x = (b.x - center.x) * scale + center.x;
        b.y = (b.y - center.y) * scale + center.y;
        c.x = (c.x - center.x) * scale + center.x;
        c.y = (c.y - center.y) * scale + center.y;
    }
    void sortascending() {
        if (a.y > b.y) {
            swapfp(a, b);
            swapfp(auv, buv);
            swapfp(aw, bw);
        }
        if (b.y > c.y) {
            swapfp(b, c);
            swapfp(buv, cuv);
            swapfp(bw, cw);
        }
        if (a.y > b.y) {
            swapfp(a, b);
            swapfp(auv, buv);
            swapfp(aw, bw);
        }
    }
};
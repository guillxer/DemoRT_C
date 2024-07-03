#pragma once

#if defined(_WINDOWS)
#include <cmath>
#endif

struct vec2
{
    float x, y;
    vec2() : x(0), y(0) {}
    vec2(float inx, float iny)
    {
        x = inx;
        y = iny;
    }
    vec2 operator+(const vec2& a) const
    {
        vec2 out;
        out.x = x + a.x;
        out.y = y + a.y;
        return out;
    }
    vec2& operator+=(const vec2& a)
    {
        this->x += a.x;
        this->y += a.y;
        return *this;
    }
    vec2& operator-=(const vec2& a)
    {
        this->x -= a.x;
        this->y -= a.y;
        return *this;
    }
    vec2& operator*=(const vec2& a)
    {
        this->x *= a.x;
        this->y *= a.y;
        return *this;
    }
    vec2& operator/=(const vec2& a)
    {
        this->x /= a.x;
        this->y /= a.y;
        return *this;
    }
    vec2& operator*=(const float& a)
    {
        this->x *= a;
        this->y *= a;
        return *this;
    }
    vec2& operator/=(const float& a)
    {
        this->x /= a;
        this->y /= a;
        return *this;
    }
    vec2 operator-(const vec2& a) const
    {
        vec2 out;
        out.x = x - a.x;
        out.y = y - a.y;
        return out;
    }
    vec2 operator*(const vec2& a) const
    {
        vec2 out;
        out.x = x * a.x;
        out.y = y * a.y;
        return out;
    }
    vec2 operator*(const float& a) const
    {
        vec2 out;
        out.x = x * a;
        out.y = y * a;
        return out;
    }
    vec2 operator/(const vec2& a) const
    {
        vec2 out;
        out.x = x / a.x;
        out.y = y / a.y;
        return out;
    }
    vec2 operator/(const float& a) const
    {
        vec2 out;
        out.x = x / a;
        out.y = y / a;
        return out;
    }
    float dot(const vec2& a) const
    {
        return x * a.x + y * a.y;
    }
    float mag() const
    {
        return sqrt(dot(*this));
    }
    vec2 getnormalized() const
    {
        return *this / mag();
    }
};

struct ivec2
{
    int x = 0, y = 0;
    ivec2() {}
    ivec2(int inx, int iny)
    {
        x = inx;
        y = iny;
    }
    ivec2 operator+(const ivec2& a) const
    {
        ivec2 out;
        out.x = x + a.x;
        out.y = y + a.y;
        return out;
    }
    ivec2 operator-(const ivec2& a) const
    {
        ivec2 out;
        out.x = x - a.x;
        out.y = y - a.y;
        return out;
    }
};

struct ivec3
{
    int x, y, z;
    ivec3() {}
    ivec3(int inx, int iny, int inz)
    {
        x = inx;
        y = iny;
        z = inz;
    }
    ivec3 operator+(const ivec3& a) const
    {
        ivec3 out;
        out.x = x + a.x;
        out.y = y + a.y;
        out.z = z + a.z;
        return out;
    }
    ivec3 operator-(const ivec3& a) const
    {
        ivec3 out;
        out.x = x - a.x;
        out.y = y - a.y;
        out.z = z - a.z;
        return out;
    }
    ivec3 operator*(const int& a) const
    {
        ivec3 out;
        out.x = x * a;
        out.y = y * a;
        out.z = z * a;
        return out;
    }
    ivec3 operator/(const int& a) const
    {
        ivec3 out;
        out.x = x / a;
        out.y = y / a;
        out.z = z / a;
        return out;
    }
    bool operator!=(const ivec3& a) const
    {
        return a.x != x ||
            a.y != y ||
            a.z != z;
    }
    bool operator==(const ivec3& a) const
    {
        return a.x == x ||
            a.y == y ||
            a.z == z;
    }
};

struct vec3
{
    float x = 0, y = 0, z = 0;
    vec3() {}
    vec3(float inx, float iny, float inz)
    {
        x = inx;
        y = iny;
        z = inz;
    }
    vec3& operator=(vec3 other)
    {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }
    vec3 operator+(const vec3& a) const
    {
        vec3 out;
        out.x = x + a.x;
        out.y = y + a.y;
        out.z = z + a.z;
        return out;
    }
    vec3 operator-(const vec3& a) const
    {
        vec3 out;
        out.x = x - a.x;
        out.y = y - a.y;
        out.z = z - a.z;
        return out;
    }
    vec3 operator-() const
    {
        vec3 out(x, y, z);
        out = out * -1.0f;
        return out;
    }
    vec3 operator*(const vec3& a) const
    {
        vec3 out;
        out.x = x * a.x;
        out.y = y * a.y;
        out.z = z * a.z;
        return out;
    }
    vec3 operator*(const float& a) const
    {
        vec3 out;
        out.x = x * a;
        out.y = y * a;
        out.z = z * a;
        return out;
    }
    vec3 operator/(const vec3& a) const
    {
        vec3 out;
        out.x = x / a.x;
        out.y = y / a.y;
        out.z = z / a.z;
        return out;
    }
    vec3 operator/(const float& a) const
    {
        vec3 out;
        out.x = x / a;
        out.y = y / a;
        out.z = z / a;
        return out;
    }
    vec3& operator+=(const vec3& a) {
        this->x += a.x;
        this->y += a.y;
        this->z += a.z;
        return *this;
    }
    vec3 cross(const vec3& a) const
    {
        vec3 out;
        out.x = y * a.z - z * a.y;
        out.y = z * a.x - x * a.z;
        out.z = x * a.y - y * a.x;
        return out;
    }
    float mag() const
    {
        return sqrt(dot(*this));
    }
    float dot(const vec3& a) const
    {
        return x * a.x + y * a.y + z * a.z;
    }
    vec3 getnormalized() const
    {
        return *this / mag();
    }
};

struct vec4
{
    float x, y, z, w;
    vec4() {}
    vec4(float inx, float iny, float inz, float inw)
    {
        x = inx;
        y = iny;
        z = inz;
        w = inw;
    }
    vec4(vec3 invec, float inw)
    {
        x = invec.x;
        y = invec.y;
        z = invec.z;
        w = inw;
    }
    vec4& operator=(vec4 a)
    {
        x = a.x;
        y = a.y;
        z = a.z;
        w = a.w;
        return *this;
    }
    vec4 operator*(const vec4& a) const
    {
        vec4 out;
        out.x = x * a.x;
        out.y = y * a.y;
        out.z = z * a.z;
        out.w = w * a.w;
        return out;
    }
    vec3 xyz() const
    {
        return vec3(x, y, z);
    }
    float dot(const vec4& a) const
    {
        return x * a.x + y * a.y + z * a.z + w * a.w;
    }
};

struct mat44
{
    vec4 rows[4];
    mat44() {}
    mat44(vec4 r0, vec4 r1, vec4 r2, vec4 r3)
    {
        rows[0] = r0;
        rows[1] = r1;
        rows[2] = r2;
        rows[3] = r3;
    }
    mat44& operator=(mat44 other)
    {
        rows[0] = other.rows[0];
        rows[1] = other.rows[1];
        rows[2] = other.rows[2];
        rows[3] = other.rows[3];
        return *this;
    }
    inline
        mat44 operator*(const mat44& mat)
    {
        mat44 matout;
        vec4 colvec0 = vec4(mat.rows[0].x, mat.rows[1].x, mat.rows[2].x, mat.rows[3].x);
        vec4 colvec1 = vec4(mat.rows[0].y, mat.rows[1].y, mat.rows[2].y, mat.rows[3].y);
        vec4 colvec2 = vec4(mat.rows[0].z, mat.rows[1].z, mat.rows[2].z, mat.rows[3].z);
        vec4 colvec3 = vec4(mat.rows[0].w, mat.rows[1].w, mat.rows[2].w, mat.rows[3].w);
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
        vec4 operator*(const vec4& vec)
    {
        vec4 vecout;
        vecout.x = rows[0].dot(vec);
        vecout.y = rows[1].dot(vec);
        vecout.z = rows[2].dot(vec);
        vecout.w = rows[3].dot(vec);
        return vecout;
    }
};

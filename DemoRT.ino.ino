#if defined(_WINDOWS)
#define RP2 0
#define MTDRAW 0
#else
#define RP2 1
#define MTDRAW 1
#endif

#define STREAM_DISPLAY 0
#define FIXED_POINT_RENDER 1
#define FIXED_POINT_CLIPPING 1

#include "math64.h"
#include "mathfp64.h"

#if RP2
#include "pico/mutex.h"
#include "pico/divider.h"
#else
#include <windows.h>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <SDL.h>
#include <winsock.h>
#pragma comment(lib, "Ws2_32.lib")
#define PI 3.14159265359f
#endif

#define MAX_ENTITIES 20

#if RP2
// The modified thumby library that includes grayscale support for OLED2
#include "ThumbySystem.h"
Thumby thumby = Thumby();
#else
SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;
#endif
// console state
#if RP2
const float consolescale = 1.0f;
#else
const float consolescale = 1.0f;
#endif
float resscale = 1.0f;
const float resscalestep = 0.01f;
const float minresscale = 1.0f;
const int hardwarewidth = 72;
const int hardwareheight = 40;
const int maximumwidth = hardwarewidth * consolescale;
const int maximumheight = hardwareheight * consolescale;
//
int width = hardwarewidth * consolescale;
int height = hardwareheight * consolescale;
int screenwidth = width;
int screenheight = height;
fixed_16_16 depthbuffer[maximumwidth * maximumheight] = {};
unsigned char colorbuffer[maximumwidth * maximumheight * 3] = {};
unsigned char compressionbuffer[maximumwidth * maximumheight * 3] = {};
float maxdrawdistance = 2.0f;
float spritedrawdistance = 3;
float entitydrawdistance = 3;
int backgroundx = 0;
int backgroundy = 0;
//
const int blockwidth = 8; // in pixels
const int blockheight = 8; // in pixels
const int bytesperendpoint = 3; // rgb
const int numendpoints = 2;
const int bytespercolorvector = blockwidth * blockheight / 8;
const int bytesperblock = bytesperendpoint * numendpoints + bytespercolorvector;
const int compressedbytesperimage = ((maximumwidth * maximumheight) / (blockwidth * blockheight)) * bytesperblock;


void rgb888to565(unsigned char r, unsigned char g, unsigned char b, unsigned char& r0, unsigned char& r1) {
    unsigned int rgb = ((r & 0b11111000) << 8) | ((g & 0b11111100) << 3) | (b >> 3);
    r0 = rgb & 0x00ff;
    r1 = (rgb & 0xff00) >> 8;
}

//https://forum.arduino.cc/t/help-converting-rgb565-to-rgb888/275681/2
void rgb565to888(unsigned char& r, unsigned char& g, unsigned char& b, unsigned char r0, unsigned char r1) {
    unsigned int color = (unsigned int)r0 | ((unsigned int)r1 >> 8);
    r = ((((color >> 11) & 0x1F) * 527) + 23) >> 6;
    g = ((((color >> 5) & 0x3F) * 259) + 33) >> 6;
    b = (((color & 0x1F) * 527) + 23) >> 6;
}

vec3fp vec3tovec3fp(vec3 a) {
    vec3fp ret;
    ret.x = (fixed_16_16)a.x;
    ret.y = (fixed_16_16)a.y;
    ret.z = (fixed_16_16)a.z;
    return ret;
}
vec3 vec3fptovec3(vec3fp a) {
    vec3 ret;
    ret.x = (float)a.x;
    ret.y = (float)a.y;
    ret.z = (float)a.z;
    return ret;
}

struct Ray;
struct Sphere;
struct Plane;

bool visiblecheck(const Ray& ray, Sphere* spheres, int numspheres, float minlength);

struct Ray {
    vec3 origin, direction;

    Ray(const vec3& origin, const vec3& direction) : origin(origin), direction(direction.getnormalized()) {}

    vec3 point_at_parameter(double t) const {
        return origin + direction * t;
    }
};
struct Sphere {
    float radius;
    vec3 center;
    Sphere(vec3 acenter, float aradius) {
        center = acenter;
        radius = aradius;
    }
    bool intersect(const Ray& ray, float& t) const {
        vec3 oc = ray.origin - center;
        double a = ray.direction.dot(ray.direction);
        double b = 2.0 * oc.dot(ray.direction);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            return false;
        }
        else {
            t = (-b - std::sqrt(discriminant)) / (2.0 * a);
            return true;
        }
    }
};
struct Plane {
    vec3 position;
    vec3 normal;
    Plane(vec3 aposition, vec3 anormal) {
        position = aposition;
        normal = anormal;
    }
    bool intersect(const Ray& ray, float& t)
    {
        float denom = normal.dot(ray.direction);
        if (denom > 0.0001) {
            vec3 p0l0 = position - ray.origin;
            t = p0l0.dot(normal) / denom;
            return (t >= 0);
        }
        return false;
    }
};

bool visiblecheck(const Ray& ray, Sphere* spheres, int numspheres, float minlength) {
    bool bvis = true;
    for (int i = 0; i < numspheres; ++i) {
        float t = 0;
        if (spheres[i].intersect(ray, t)) {
            if (t <= minlength) {
                bvis = false;
            }
        }
    }
    return bvis;
}

vec3 lightpos = vec3(0, -0.4f, -1);

void raytracetest(int coreid) {
    Sphere lightsphere = Sphere(lightpos, 0.1);
    Sphere sphere0 = Sphere(vec3(0.5, 0, -1.0), 0.5);
    Sphere sphere1 = Sphere(vec3(-0.8, 0, -1.5), 0.5);
    Sphere spheres[] = { lightsphere, sphere0, sphere1 };
    Sphere visspheres0[] = { sphere1 };
    Sphere visspheres1[] = { sphere0 };
    Sphere visspheres2[] = { sphere0, sphere1 };
    Sphere* visspheres[] = { nullptr, visspheres0, visspheres1, visspheres2 };
    int visspheresnum[] = {
        0,
        sizeof(visspheres0) / sizeof(Sphere),
        sizeof(visspheres1) / sizeof(Sphere),
        sizeof(visspheres2) / sizeof(Sphere) };
    const int numspheres = sizeof(spheres) / sizeof(Sphere);
    Plane plane0 = Plane(vec3(0, 1, 0), vec3(0, 1, 0));
    Plane planes[] = { plane0 };
    const int numplanes = sizeof(planes) / sizeof(Plane);

    float aspectratio = (float)width / height;

    int starty = coreid == 0 ? height - 1 : height / 2 - 1;
    int endy = coreid == 0 ? height / 2 : 0;

    float visbias = 0.001f;

    for (int j = starty; j >= endy; --j) {
        for (int i = 0; i < width; ++i) {
            float u = (float)i / width;
            float v = (float)j / height;
            Ray ray(vec3(0, 0, 0), vec3((2.0f * u - 1.0f) * aspectratio, 2.0f * v - 1.0f, -1));

            int r = 0;
            int g = 0;
            int b = 0;
            float t = 0;
            float mint = 100000;
            for (int spherei = 0; spherei < numspheres; ++spherei) {
                if (spheres[spherei].intersect(ray, t)) {
                    if (spherei == 0) {
                        if (t < mint) {
                            mint = t;
                            r = 255;
                            g = 255;
                            b = 255;
                        }
                    }
                    else {
                        if (t < mint) {
                            vec3 point = ray.point_at_parameter(t);
                            vec3 normal = (point - spheres[spherei].center).getnormalized();
                            vec3 lightvec = -(point - lightpos).getnormalized();
                            float length = sqrt((point - lightpos).dot(point - lightpos));
                            Ray visray(point + normal * visbias, lightvec);
                            mint = t;
                            if (visiblecheck(visray, visspheres[spherei], visspheresnum[spherei], length)) {
                                float color = (lightvec.dot(normal) + 1.0f) / 2.0f;
                                r = int(255.99 * color);
                                g = int(255.99 * color);
                                b = int(255.99 * color);
                            }
                        }
                    }
                }
            }
            for (int planei = 0; planei < numplanes; ++planei) {
                if (planes[planei].intersect(ray, t)) {
                    if (t < mint) {
                        vec3 point = ray.point_at_parameter(t);
                        vec3 normal = planes[planei].normal;
                        vec3 lightvec = (point - lightpos).getnormalized();
                        float length = sqrt((point - lightpos).dot(point - lightpos));
                        Ray visray(point + normal * visbias, lightvec);
                        mint = t;
                        if (visiblecheck(visray, visspheres[3], visspheresnum[3], length)) {
                            float color = (lightvec.dot(normal) + 1.0f) / 2.0f;
                            r = int(255.99 * color);
                            g = int(255.99 * color);
                            b = int(255.99 * color);
                        }
                    }
                }
            }
            colorbuffer[i * height * 3 + j * 3 + 0] = r;
            colorbuffer[i * height * 3 + j * 3 + 1] = g;
            colorbuffer[i * height * 3 + j * 3 + 2] = b;
        }
    }
}

int lumantsc(vec3 rgb) {
    return 0.3f * rgb.x + 0.59f * rgb.y + 0.11f * rgb.z;
}

void getcccvectoraddress(int linearblockaddress, int& byteoffset, int& bitoffset) {
    byteoffset = linearblockaddress / 8;
    bitoffset = linearblockaddress % 8;
}

void compresscolorbuffer() {
    int cccaddr = 0;
    unsigned char* cccbuffer = compressionbuffer;
    for (int x = 0; x < maximumwidth; x += blockwidth) {
        for (int y = 0; y < maximumheight; y += blockheight) {
            // for each block compute the mean luminance
            float meanluma = 0;
            float pixelluma[blockwidth * blockheight] = {};
            int lumaaddr = 0;
            for (int bx = x; bx < x + blockwidth; ++bx) {
                for (int by = y; by < y + blockheight; ++by) {
                    unsigned char r = colorbuffer[bx * maximumheight * 3 + by * 3 + 0];
                    unsigned char g = colorbuffer[bx * maximumheight * 3 + by * 3 + 1];
                    unsigned char b = colorbuffer[bx * maximumheight * 3 + by * 3 + 2];
                    pixelluma[lumaaddr] = lumantsc(vec3(r, g, b));
                    meanluma += pixelluma[lumaaddr];
                    lumaaddr++;
                }
            }
            meanluma /= blockwidth * blockheight;
            // for each pixel determine if the pixel is smaller or larger than the median and mark each pixel in group 0 or 1 respectively
            lumaaddr = 0;
            char cccvector[bytespercolorvector] = {};
            vec3 endpoint0 = {};
            vec3 endpoint1 = {};
            int num0pixels = 0;
            int num1pixels = 0;
            for (int bx = x; bx < x + blockwidth; ++bx) {
                for (int by = y; by < y + blockheight; ++by) {
                    float luma = pixelluma[lumaaddr];
                    // for each group 0 or 1 compute the mean color, these are endpoints 0 and 1
                    int byteoffset = 0;
                    int bitoffset = 0;
                    int vectorvalue = 0;
                    getcccvectoraddress(lumaaddr, byteoffset, bitoffset);
                    if (luma > meanluma) {
                        endpoint1.x += colorbuffer[bx * maximumheight * 3 + by * 3 + 0];
                        endpoint1.y += colorbuffer[bx * maximumheight * 3 + by * 3 + 1];
                        endpoint1.z += colorbuffer[bx * maximumheight * 3 + by * 3 + 2];
                        vectorvalue = 1;
                        num1pixels++;
                    }
                    else {
                        endpoint0.x += colorbuffer[bx * maximumheight * 3 + by * 3 + 0];
                        endpoint0.y += colorbuffer[bx * maximumheight * 3 + by * 3 + 1];
                        endpoint0.z += colorbuffer[bx * maximumheight * 3 + by * 3 + 2];
                        vectorvalue = 0;
                        num0pixels++;
                    }
                    if (vectorvalue > 0) {
                        cccvector[byteoffset] |= 1 << bitoffset;
                    }
                    lumaaddr++;
                }
            }
            if (num0pixels > 0) {
                endpoint0 = endpoint0 / num0pixels;
            }
            if (num1pixels > 0) {
                endpoint1 = endpoint1 / num1pixels;
            }
            // pack the data as follows: endpoint0.x, endpoint0.y, endpoint0.z, endpoint1.x, endpoint1.y, endpoint1.z, cccvector
            cccbuffer[cccaddr++] = endpoint0.x;
            cccbuffer[cccaddr++] = endpoint0.y;
            cccbuffer[cccaddr++] = endpoint0.z;
            cccbuffer[cccaddr++] = endpoint1.x;
            cccbuffer[cccaddr++] = endpoint1.y;
            cccbuffer[cccaddr++] = endpoint1.z;
            for (int i = 0; i < bytespercolorvector; ++i) {
                cccbuffer[cccaddr++] = cccvector[i];
            }
        }
    }
}

void decompresscolorbuffer() {
    int cccaddr = 0;
    unsigned char* cccbuffer = compressionbuffer;
    // for each block
    for (int x = 0; x < maximumwidth; x += blockwidth) {
        for (int y = 0; y < maximumheight; y += blockheight) {
            // store endpoints 0 and 1
            vec3 endpoint0(cccbuffer[cccaddr + 0], cccbuffer[cccaddr + 1], cccbuffer[cccaddr + 2]);
            vec3 endpoint1(cccbuffer[cccaddr + 3], cccbuffer[cccaddr + 4], cccbuffer[cccaddr + 5]);
            cccaddr += 6;
            char cccvector[bytespercolorvector] = {};
            for (int i = 0; i < bytespercolorvector; ++i) {
                cccvector[i] = cccbuffer[cccaddr++];
            }
            int linearblockaddress = 0;
            for (int bx = x; bx < x + blockwidth; ++bx) {
                for (int by = y; by < y + blockheight; ++by) {
                    // for each pixel in the block grab the correct endpoint by reading from the cccvector
                    int byteoffset = 0;
                    int bitoffset = 0;
                    getcccvectoraddress(linearblockaddress, byteoffset, bitoffset);
                    int vectorvalue = cccvector[byteoffset] & 1 << bitoffset;
                    vec3 selectedendpoint = endpoint0;
                    if (vectorvalue > 0) {
                        selectedendpoint = endpoint1;
                    }
                    colorbuffer[bx * maximumheight * 3 + by * 3 + 0] = selectedendpoint.x;
                    colorbuffer[bx * maximumheight * 3 + by * 3 + 1] = selectedendpoint.y;
                    colorbuffer[bx * maximumheight * 3 + by * 3 + 2] = selectedendpoint.z;
                    linearblockaddress++;
                }
            }
        }
    }
}


volatile float playerangle = 90.0f;
float lastplayerangle = playerangle;
float playerspeed = 1.0f;
float playerturnspeed = 180.0f;
float playerturnspeedair = playerturnspeed / 2.0f;
vec3 playerpos = vec3(0, 0.0, 0.0);
vec3 lastplayerpos = playerpos;


vec3 angletovector(float angle) {
    float rads = PI * angle / 180.0f;
    return vec3(sin(rads), 0.0f, -cos(rads));
}

// proj mat calc
constexpr float ratio = (float)maximumheight / (float)maximumwidth;
float focalLength = 2.0f;
#undef near
#undef far
#define near 0.01f


const int dither88[8][8] = {
    { 0, 32, 8, 40, 2, 34, 10, 42}, /* 8x8 Bayer ordered dithering */
    {48, 16, 56, 24, 50, 18, 58, 26}, /* pattern. Each input pixel */
    {12, 44, 4, 36, 14, 46, 6, 38}, /* is scaled to the 0..63 range */
    {60, 28, 52, 20, 62, 30, 54, 22}, /* before looking in this table */
    { 3, 35, 11, 43, 1, 33, 9, 41}, /* to determine the action. */
    {51, 19, 59, 27, 49, 17, 57, 25},
    {15, 47, 7, 39, 13, 45, 5, 37},
    {63, 31, 55, 23, 61, 29, 53, 21} };

const float dither44[4][4] = {
    { 00.0 / 16.0, 12.0 / 16.0, 03.0 / 16.0, 15.0 / 16.0 },
    { 08.0 / 16.0, 04.0 / 16.0, 11.0 / 16.0, 07.0 / 16.0 },
    { 02.0 / 16.0, 14.0 / 16.0, 01.0 / 16.0, 13.0 / 16.0 },
    { 10.0 / 16.0, 06.0 / 16.0, 09.0 / 16.0, 05.0 / 16.0 } };

#define DITHER8X8 0
#define MONO 1
// 8x8 dither taken from lib retro projects
float find_closest(int x, int y, float c0)
{
    float limit = 0.0;
#if DITHER8X8
    if (x < 8)
#else
    if (x < 4)
#endif
    {
#if DITHER8X8
        limit = (dither88[x][y] + 1) / 64.0;
#else
        limit = (dither44[x][y]);
#endif
    }
#if MONO
    if (c0 <= limit)
        return 0.0;
    return 1.0;
#else
    float a = 0.33f;
    float b = 0.0f;
    if (c0 >= 0.66f) {
        a = 1.0f;
        b = 0.66f;
    }
    else if (c0 >= 0.33f) {
        a = 0.66f;
        b = 0.33f;
    }
    if (c0 >= limit) {
        return a;
    }
    else {
        return b;
    }
#endif
}

void dithershadeintodisplay() {
    for (int x = 0; x < maximumwidth; ++x) {
        for (int y = 0; y < maximumheight; ++y) {
            float shadeestimate =
                (colorbuffer[x * maximumheight * 3 + y * 3 + 0] +
                    colorbuffer[x * maximumheight * 3 + y * 3 + 1] +
                    colorbuffer[x * maximumheight * 3 + y * 3 + 2]) / 3.0f;
#if DITHER8X8
            float color = find_closest(x % 8, y % 8, shadeestimate / 255.0f / 2.0f);
#else
            float color = find_closest(x % 4, y % 4, shadeestimate / 255.0f / 2.0f);
#endif
            color *= 255;
#if RP2
            if (color > 0) {
                thumby.drawPixel(x / consolescale, y / consolescale, 1);
            }
#else
#if 0
            SDL_SetRenderDrawColor(
                renderer,
                colorbuffer[x * maximumheight * 3 + y * 3 + 0],
                colorbuffer[x * maximumheight * 3 + y * 3 + 1],
                colorbuffer[x * maximumheight * 3 + y * 3 + 2],
                255);
#else
            color = (color > 0) ? 255 : 0;
            SDL_SetRenderDrawColor(
                renderer,
                color,
                color,
                color,
                255);
#endif
            SDL_RenderDrawPoint(renderer, x, y);
#endif
        }
    }
}

// Projects a 3D vector from world space to screen space
void project(
    vec3 vertex,
    mat44 viewprojmat,
    vec2& screenvecout,
    bool& isvisible,
    float& distance) {
    vec4 vertexw = vec4(vertex.x, vertex.y, vertex.z, 1.0f);
    vec4 screenvec = viewprojmat * vertexw;
    float w = screenvec.w;
    bool visible = true;
    if (screenvec.z > maxdrawdistance) {
        visible = false;
    }
    distance = screenvec.z;
    vec2 projvec = { screenvec.x / w, screenvec.y / w };
    float x = (projvec.x + 1.0f) * 0.5f * width;
    float y = (1.0f - (projvec.y + 1.0f) * 0.5f) * height;
    screenvecout = vec2(x, y);
    isvisible = visible;
}

// Creates a Euler rotation matrix from yaw, pitch, and roll
mat44 makeeuler(float yaw, float pitch, float roll) {
    float alpha = yaw * static_cast<float>(PI) / 180.0f;
    float beta = pitch * static_cast<float>(PI) / 180.0f;
    float gamma = roll * static_cast<float>(PI) / 180.0f;
    mat44 eulermat = mat44(
        vec4(std::cos(beta) * std::cos(gamma),
            std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha) * std::sin(gamma),
            std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha) * std::sin(gamma),
            0.0f),
        vec4(std::cos(beta) * std::sin(gamma),
            std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha) * std::cos(gamma),
            std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha) * std::cos(gamma),
            0.0f),
        vec4(-std::sin(beta),
            std::sin(alpha) * std::cos(beta),
            std::cos(alpha) * std::cos(beta),
            0.0f),
        vec4(0.0f,
            0.0f,
            0.0f,
            1.0f));
    return eulermat;
}

mat44fp makeeulerfp(float yaw, float pitch, float roll) {
    float alpha = yaw * static_cast<float>(PI) / 180.0f;
    float beta = pitch * static_cast<float>(PI) / 180.0f;
    float gamma = roll * static_cast<float>(PI) / 180.0f;
    mat44fp eulermat = mat44fp(
        vec4fp((fixed_16_16)std::cos(beta) * (fixed_16_16)std::cos(gamma),
            (fixed_16_16)std::sin(alpha) * (fixed_16_16)std::sin(beta) * (fixed_16_16)std::cos(gamma) - (fixed_16_16)std::cos(alpha) * (fixed_16_16)std::sin(gamma),
            (fixed_16_16)std::cos(alpha) * (fixed_16_16)std::sin(beta) * (fixed_16_16)std::cos(gamma) + (fixed_16_16)std::sin(alpha) * (fixed_16_16)std::sin(gamma),
            (fixed_16_16)0.0f),
        vec4fp((fixed_16_16)std::cos(beta) * (fixed_16_16)std::sin(gamma),
            (fixed_16_16)std::sin(alpha) * (fixed_16_16)std::sin(beta) * (fixed_16_16)std::sin(gamma) + (fixed_16_16)std::cos(alpha) * (fixed_16_16)std::cos(gamma),
            (fixed_16_16)std::cos(alpha) * (fixed_16_16)std::sin(beta) * (fixed_16_16)std::sin(gamma) - (fixed_16_16)std::sin(alpha) * (fixed_16_16)std::cos(gamma),
            (fixed_16_16)0.0f),
        vec4fp((fixed_16_16)-std::sin(beta),
            (fixed_16_16)std::sin(alpha) * (fixed_16_16)std::cos(beta),
            (fixed_16_16)std::cos(alpha) * (fixed_16_16)std::cos(beta),
            (fixed_16_16)0.0f),
        vec4fp((fixed_16_16)0.0f,
            (fixed_16_16)0.0f,
            (fixed_16_16)0.0f,
            (fixed_16_16)1.0f));
    return eulermat;
}

template<class T>
void swap(T& a, T& b) {
    T temp = a;
    a = b;
    b = temp;
}
template<class T>
void swap(T*& a, T*& b) {
    T* temp = a;
    a = b;
    b = temp;
}

struct Triangle {
    vec3 a, b, c;
    Triangle() {}
    Triangle(vec3 ain, vec3 bin, vec3 cin) {
        a = ain; b = bin; c = cin;
    }
    vec3 facenormal()
    {
        vec3 cross = (a - b).cross(c - a);
        return cross.getnormalized();
    }
    // Quick sphere approximation check, not an optimal sphere
    // return true on possible overlap
    bool spherecheck(vec3 spherepos, float sphereradius)
    {
        float distancesquared = max((a - b).dot(a - b), (a - c).dot((a - c)));
        return ((a - spherepos).dot((a - spherepos)) - distancesquared * distancesquared) < sphereradius * sphereradius;
    }
    // Real Time Collision Detection(The Morgan Kaufmann Series in Interactive 3 - D Technology)
    // optimized version from blog, removed divides kept early outs since rp doen't support intrinsics
    bool spherecollision(vec3 position, float radius, float& distance) {
        distance = 1000.0f;
        bool separated = false;
        vec3 A = a;
        vec3 B = b;
        vec3 C = c;
        vec3 P = position;
        float r = radius;
        // Translate problem so sphere is centered at origin
        A = A - P;
        B = B - P;
        C = C - P;
        //###########################################
        // test triangle plane
        //###########################################
        // compute a vector normal to triangle plane(V)
        vec3 V = (B - A).cross(C - A);
        // get distance from the plane using the first vertex
        float d = A.dot(V);
        float e = V.dot(V);
        separated = d * d > r * r * e;
        if (separated) return true;
        //###########################################
        // test traingle verts
        //##########################################
        // vert is not inside sphere and sphere does not intersect plane between each vertex
        float aa = A.dot(A);
        float ab = A.dot(B);
        float ac = A.dot(C);
        float bb = B.dot(B);
        float bc = B.dot(C);
        float cc = C.dot(C);
        float rr = r * r;
        bool separateda = (aa > rr) && (ab > aa) && (ac > aa);
        bool separatedb = (bb > rr) && (ab > bb) && (bc > bb);
        bool separatedc = (cc > rr) && (ac > cc) && (bc > cc);
        separated = separated || separateda || separatedb || separatedc;
        if (separated) return true;
        //###########################################
        // test triangle edges
        //###########################################
        // project to closest edge on point and check for any separations
        vec3 AB = (B - A);
        vec3 BC = (C - B);
        vec3 CA = (A - C);
        float d1 = A.dot(AB);
        float e1 = AB.dot(AB);
        float d2 = B.dot(BC);
        float e2 = BC.dot(BC);
        float d3 = C.dot(CA);
        float e3 = CA.dot(CA);
        vec3 Q1 = A * e1 - AB * d1;
        vec3 QC = C * e1 - Q1;
        vec3 Q2 = B * e2 - BC * d2;
        vec3 QA = A * e2 - Q2;
        vec3 Q3 = C * e3 - CA * d3;
        vec3 QB = B * e3 - Q3;
        bool separated1 = (Q1.dot(Q1) > r * r * e1 * e1) && (Q1.dot(QC) > 0);
        bool separated2 = (Q2.dot(Q2) > r * r * e2 * e2) && (Q2.dot(QA) > 0);
        bool separated3 = (Q3.dot(Q3) > r * r * e3 * e3) && (Q3.dot(QB) > 0);
        separated = separated || separated1 || separated2 || separated3;
        vec3 N = V.getnormalized();
        distance = std::abs(A.dot(N));
        return separated;
    }
};

struct TriangleUV : public Triangle {
    vec2 auv, buv, cuv;
    float aw, bw, cw;
    TriangleUV() {}
    TriangleUV(vec3 ain, vec3 bin, vec3 cin,
        vec2 auvin, vec2 buvin, vec2 cuvin,
        float awin, float bwin, float cwin) {
        a = ain; b = bin; c = cin;
        auv = auvin; buv = buvin; cuv = cuvin;
        aw = awin; bw = bwin; cw = cwin;
    }
    void sortascending() {
        if (a.y > b.y) {
            swap(a, b);
            swap(auv, buv);
            swap(aw, bw);
        }
        if (b.y > c.y) {
            swap(b, c);
            swap(buv, cuv);
            swap(bw, cw);
        }
        if (a.y > b.y) {
            swap(a, b);
            swap(auv, buv);
            swap(aw, bw);
        }
    }
};

struct Texture {
    const unsigned char* data;
    int width;
    int height;
    Texture(const unsigned char* datain, int widthin, int heightin) {
        data = datain;
        width = widthin;
        height = heightin;
    }
};
Texture boundtexture0 = Texture(nullptr, 0, 0);
Texture boundtexture1 = Texture(nullptr, 0, 0);


void clearscenedepth() {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            colorbuffer[x * height * 3 + y * 3 + 0] = 0;
            colorbuffer[x * height * 3 + y * 3 + 1] = 0;
            colorbuffer[x * height * 3 + y * 3 + 2] = 0;
        }
    }
}

template<class T>
T clamp(T a, T mina, T maxa) {
    T r = max(mina, a);
    r = min(maxa, r);
    return r;
}

struct Edge {
    vec2 a;
    vec2 b;
    Edge() {}
    Edge(vec2 ain, vec2 bin) {
        a = ain;
        b = bin;
    }
};

struct Edge3 {
    vec3 a;
    vec3 b;
    Edge3() {}
    Edge3(vec3 ain, vec3 bin) {
        a = ain;
        b = bin;
    }
};
struct Edge3UV {
    vec3 a;
    vec3 b;
    vec2 uva;
    vec2 uvb;
    float wa;
    float wb;
    Edge3UV() {}
    Edge3UV(vec3 ain, vec3 bin, vec2 auvin, vec2 buvin, float wain, float wbin) {
        a = ain;
        b = bin;
        uva = auvin;
        uvb = buvin;
        wa = wain;
        wb = wbin;
    }
};
struct Edge3UVFP {
    vec3fp a;
    vec3fp b;
    vec2fp uva;
    vec2fp uvb;
    fixed_16_16 wa;
    fixed_16_16 wb;
    Edge3UVFP() {}
    Edge3UVFP(vec3fp ain, vec3fp bin, vec2fp auvin, vec2fp buvin, fixed_16_16 wain, fixed_16_16 wbin) {
        a = ain;
        b = bin;
        uva = auvin;
        uvb = buvin;
        wa = wain;
        wb = wbin;
    }
};

// Return signed shortest distance from point to plane, plane normal must be normalised
inline float vertexplanedistance(vec3 vertex, vec3 planepoint, vec3 planenormal) {
    return planenormal.dot(vertex) - planenormal.dot(planepoint);
};

// Return signed shortest distance from point to plane, plane normal must be normalised
inline fixed_16_16 vertexplanedistancefp(vec3fp vertex, vec3fp planepoint, vec3fp planenormal) {
    return planenormal.dot(vertex) - planenormal.dot(planepoint);
};

float min3(float a, float b, float c) {
    return min(min(a, b), c);
}

float max3(float a, float b, float c) {
    return max(max(a, b), c);
}

#if !RP2
SOCKET m_socket;
#endif

#define STREAMCOMPRESSED 1
#if STREAMCOMPRESSED
void sendframebuffer() {
    // Send and receive data.
#if !RP2
    int bytesSent;
    int bytesRecv = SOCKET_ERROR;
#endif
    // Be careful with the array bound, provide some checking mechanism
    // colorbuffer
    char recvbuf[200] = "";
    int sendchunksize = compressedbytesperimage;
    char sendbuf[compressedbytesperimage];
    for (int i = 0; i < compressedbytesperimage; i += sendchunksize) {
        unsigned char* compressedchunk = &compressionbuffer[i];
#if RP2
        // Wait to begin transmission
        char donebuf[4];
        Serial.readBytes(donebuf, sizeof(donebuf));
        while (donebuf[0] != 77) {
            Serial.readBytes(donebuf, sizeof(donebuf));
        }
        Serial.write(compressedchunk, sendchunksize);
        Serial.flush();
#else
        bytesSent = send(m_socket, (char*)(compressedchunk), sendchunksize, 0);
        bytesRecv = recv(m_socket, recvbuf, sizeof(recvbuf), 0);
#endif
    }
}
#else
void sendframebuffer() {
    // Send and receive data.
#if !RP2
    int bytesSent;
    int bytesRecv = SOCKET_ERROR;
#endif
    // Be careful with the array bound, provide some checking mechanism
    // colorbuffer
    char recvbuf[200] = "";
    char sendbuf[maximumwidth * maximumheight * 3];
    for (int x = 0; x < width; ++x) {
        volatile unsigned char* framebufferrow = &colorbuffer[x * height * 3 + 0 * 3 + 0];
#if RP2
        // Wait to begin transmission
        char donebuf[4];
        Serial.readBytes(donebuf, sizeof(donebuf));
        while (donebuf[0] != 77) {
            Serial.readBytes(donebuf, sizeof(donebuf));
        }
        //Serial.write(framebufferrow, height * 3);
        Serial.flush();
#else
        bytesSent = send(m_socket, (char*)(framebufferrow), height * 3, 0);
        bytesRecv = recv(m_socket, recvbuf, sizeof(recvbuf), 0);
#endif
    }
}
#endif

void initnetwork() {
#if RP2
    Serial.begin(115200);
    Serial.setTimeout(0);
#else
    // Initialize Winsock.
    WSADATA wsaData;
    int iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
    if (iResult != NO_ERROR)
        printf("Client: Error at WSAStartup().\n");
    else
        printf("Client: WSAStartup() is OK.\n");
    // Create a socket.
    m_socket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (m_socket == INVALID_SOCKET) {
        printf("Client: socket() - Error at socket(): %ld\n", WSAGetLastError());
        WSACleanup();
        return;
    }
    else
        printf("Client: socket() is OK.\n");
    // Connect to a server.
    sockaddr_in clientService;
    clientService.sin_family = AF_INET;
    clientService.sin_addr.s_addr = inet_addr("127.0.0.1");
    clientService.sin_port = htons(65432);
    if (connect(m_socket, (SOCKADDR*)&clientService, sizeof(clientService)) == SOCKET_ERROR) {
        printf("Client: connect() - Failed to connect.\n");
        WSACleanup();
        return;
    }
#endif
}

void initrender() {
#if RP2
#else
    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
    }
    // Create a window
    window = SDL_CreateWindow("SDL Draw Line", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
    if (window == NULL) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
    }
    // Create a renderer
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == NULL) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
    }
#endif
}

float lastUpdateTime = 0.0f;
float deltatime = 0.0f;
#if RP2
void displayFPS() {
    thumby.setCursor(0, 0);
    thumby.print(1000.0f / (millis() - lastUpdateTime));
    deltatime = (millis() - lastUpdateTime) / 1000.0f;
    lastUpdateTime = millis();
}
#endif

void clearscreen(int value)
{
#if RP2
    thumby.clear();
#else
    SDL_PumpEvents();
    // clear rt
    SDL_SetRenderDrawColor(renderer, value, value, value, 255);
    // Clear the screen
    SDL_RenderClear(renderer);
#endif
}

void updateframe()
{
#if RP2
    displayFPS();
    thumby.writeBuffer(thumby.getBuffer(), thumby.getBufferSize());
#else
    SDL_RenderPresent(renderer);
#endif
}

static bool renderwork = false;
bool is_running = true;
volatile bool startstreamingcore1 = false;
volatile bool startrendercore1 = false;
bool renderinitialized = false;


void loopdraw(int corenum);
#if !RP2
void loop();
void loop1();
Uint32 mTicksCount = SDL_GetTicks();
int main(int argc, char** argv) {
    initnetwork();
    initrender();
    while (is_running) {
        loop();
    }
    SDL_Quit();
    exit(0);
    return 0;
}
#else
void setup() {
    auto_init_mutex(rendermtx);
    // Sets up buttons, audio, link pins, and screen
    thumby.begin();
    // Init duplex UART for Thumby to PC comms
    //Serial.begin(115200);
    //Serial.setTimeout(0);
    Serial.print("create mtx");
    initrender();
    renderinitialized = true;
}
#endif

void loop() {
    clearscreen(0);

#if RP2
    float dt = deltatime;
#else
    float dt = (SDL_GetTicks() - mTicksCount) / 1000.0f;
    mTicksCount = SDL_GetTicks();
#endif
    // dynamic resolution calculation
    if (dt > 0.035f) {
        resscale -= resscalestep;
    }
    else if (dt < 0.030f) {
        resscale += resscalestep;
    }
    resscale = clamp(resscale, minresscale, 1.0f);
    width = maximumwidth * resscale;
    height = maximumheight * resscale;
    dt = min(dt, 0.1f);

    // update player position
    bool wpress = false;
    bool spress = false;
    bool apress = false;
    bool dpress = false;
    bool abutton = false;
    bool bbutton = false;

#if RP2
    if (thumby.isPressed(BUTTON_U))
    {
        backgroundy++;
        wpress = true;
    }
    if (thumby.isPressed(BUTTON_D))
    {
        backgroundy--;
        spress = true;
    }
    if (thumby.isPressed(BUTTON_L))
    {
        backgroundx++;
        apress = true;
    }
    if (thumby.isPressed(BUTTON_R))
    {
        backgroundx--;
        dpress = true;
    }
    if (thumby.isPressed(BUTTON_A))
    {
        abutton = true;
    }
    if (thumby.isPressed(BUTTON_B))
    {
        bbutton = true;
    }
#else
    int numkeys = 0;
    SDL_PumpEvents();
    const Uint8* keystate = SDL_GetKeyboardState(&numkeys);
    if (keystate[SDL_SCANCODE_W])
    {
        wpress = true;
    }
    if (keystate[SDL_SCANCODE_S])
    {
        spress = true;
    }
    if (keystate[SDL_SCANCODE_A])
    {
        apress = true;
    }
    if (keystate[SDL_SCANCODE_D])
    {
        dpress = true;
    }
    if (keystate[SDL_SCANCODE_PERIOD])
    {
        abutton = true;
    }
    if (keystate[SDL_SCANCODE_COMMA])
    {
        bbutton = true;
    }
#endif
    float movementscale = 2.0f;
    if (wpress) {
        lightpos.z -= movementscale * dt;
    }
    if (spress) {
        lightpos.z += movementscale * dt;
    }
    if (apress) {
        lightpos.x -= movementscale * dt;
    }
    if (dpress) {
        lightpos.x += movementscale * dt;
    }
    if (abutton) {
        lightpos.y += movementscale * dt;
    }
    if (bbutton) {
        lightpos.y -= movementscale * dt;
    }


    lightpos.x = clamp(lightpos.x, -2.0f, 2.0f);
    lightpos.z = clamp(lightpos.z, -2.0f, -0.1f);
    lightpos.y = clamp(lightpos.y, -1.0f, 1.0f);


#if STREAM_DISPLAY
    // process physics and distplay streaming at the same time
    while (startstreamingcore1) {}
#endif

    // start rendering after this clear
    clearscenedepth();
#if RP2 && MTDRAW
    startrendercore1 = true;
    loopdraw(0 /*corenum*/);
    while (startrendercore1) {}
#else
    loopdraw(0 /*corenum*/);
    loopdraw(1 /*corenum*/);
#endif

#if RP2 && MTDRAW
    // 
#endif
    //upscaletovirtualresolution();
#if STREAM_DISPLAY
    // begin streaming work into the start of the next frame
    startstreamingcore1 = true;
#endif
    dithershadeintodisplay();
    // update frame
    updateframe();
}

void loopdraw(int corenum) {
    raytracetest(corenum);
}

void setup1() {}
void loop1() {
    while (!startrendercore1) {}
    loopdraw(1);
    startrendercore1 = false;
}


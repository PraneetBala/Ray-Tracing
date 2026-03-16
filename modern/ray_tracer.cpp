// =============================================================================
// Ray Tracer — Modernized C++17
// Based on the original implementation (2019), refactored for clarity and
// idiomatic C++17: Vec3 math type, constexpr constants, std::filesystem,
// std::clamp, range-based loops, and named structs throughout.
//
// Build:  make
// Run:    ./ray_tracer <scene.txt>
// Output: <scene>.ppm
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
constexpr float PI      = 3.14159265f;
constexpr float EPSILON = 9e-5f;
constexpr int   MAX_DEPTH = 10;
constexpr float INF     = 1e5f;

// -----------------------------------------------------------------------------
// Vec3 — 3-component float vector with operator overloads
// -----------------------------------------------------------------------------
struct Vec3 {
    float x = 0, y = 0, z = 0;

    Vec3() = default;
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3  operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3  operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3  operator*(float s)       const { return {x*s,   y*s,   z*s};   }
    Vec3  operator/(float s)       const { return {x/s,   y/s,   z/s};   }
    Vec3  operator-()              const { return {-x,    -y,    -z};     }
    Vec3& operator+=(const Vec3& o) { x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vec3& operator*=(float s)       { x*=s;   y*=s;   z*=s;   return *this; }

    float dot(const Vec3& o)  const { return x*o.x + y*o.y + z*o.z; }
    float lengthSq()          const { return dot(*this); }
    float length()            const { return std::sqrt(lengthSq()); }

    Vec3 cross(const Vec3& o) const {
        return { y*o.z - z*o.y,
                 z*o.x - x*o.z,
                 x*o.y - y*o.x };
    }

    Vec3 normalized() const {
        float l = length();
        return l > 0.f ? *this / l : *this;
    }
};

inline Vec3 operator*(float s, const Vec3& v) { return v * s; }

// -----------------------------------------------------------------------------
// Scene data structures
// -----------------------------------------------------------------------------
struct Ray {
    Vec3 origin, dir;
};

struct Material {
    Vec3  diffuse;          // od rgb
    Vec3  specular;         // os rgb
    float ka, kd, ks, n;   // ambient, diffuse, specular coefficients + shininess
    float alpha;            // opacity  (1 = opaque)
    float eta;              // refractive index
};

struct Sphere {
    Vec3  center;
    float radius;
    int   materialIdx;
    int   textureIdx;       // -1 = no texture
};

struct Light {
    Vec3  pos;
    float w;                // 0 = directional, 1 = point
    Vec3  color;
};

struct Face {
    int vx, vy, vz;        // vertex indices
    int nx, ny, nz;        // normal indices  (-1 = flat)
    int tx, ty, tz;        // uv indices      (-1 = none)
    bool smoothShading;
    bool hasTexture;
    int  materialIdx;
    int  textureIdx;
};

struct Texture {
    int w = 0, h = 0;
    std::vector<float> r, g, b;
};

// -----------------------------------------------------------------------------
// Global scene state (mirrors original global layout)
// -----------------------------------------------------------------------------
struct Scene {
    Vec3  eye, viewDir, upDir, bkgColor;
    float vfov   = 60.f;
    int   width  = 512;
    int   height = 512;

    std::vector<Sphere>   spheres;
    std::vector<Material> materials;
    std::vector<Light>    lights;
    std::vector<Vec3>     vertices;
    std::vector<Vec3>     normals;
    std::vector<Vec3>     uvCoords;   // stored as (u, v, 0)
    std::vector<Face>     faces;
    std::vector<Texture>  textures;
};

// -----------------------------------------------------------------------------
// Texture lookup helpers
// -----------------------------------------------------------------------------
static Vec3 texelAt(const Texture& tex, int i, int j) {
    int idx = i + tex.w * j;
    idx = std::clamp(idx, 0, (int)tex.r.size() - 1);
    return { tex.r[idx] / 255.f, tex.g[idx] / 255.f, tex.b[idx] / 255.f };
}

// -----------------------------------------------------------------------------
// Scene file parser
// -----------------------------------------------------------------------------
static std::vector<std::string> split(const std::string& line) {
    std::istringstream ss(line);
    std::vector<std::string> tokens;
    std::string tok;
    while (ss >> tok) tokens.push_back(tok);
    return tokens;
}

static Texture loadTexturePPM(const std::string& path) {
    Texture tex;
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open texture: " << path << "\n"; std::exit(1); }

    std::string line;
    std::getline(f, line);              // header line: "format w h maxval"
    auto hdr = split(line);
    tex.w = std::stoi(hdr[1]);
    tex.h = std::stoi(hdr[2]);

    std::string val;
    while (std::getline(f, val)) {
        if (val.empty()) continue;
        tex.r.push_back(std::stof(val));
        std::getline(f, val); tex.g.push_back(std::stof(val));
        std::getline(f, val); tex.b.push_back(std::stof(val));
    }
    return tex;
}

// Parse "v1/t1/n1" face token into (vertIdx, texIdx, normIdx) — all 0-based
static void parseFaceToken(const std::string& tok,
                           int& vi, int& ti, int& ni,
                           bool& hasTex, bool& hasNorm) {
    auto s1 = tok.find('/');
    vi = std::stoi(tok.substr(0, s1)) - 1;
    ti = ni = -1;
    hasTex = hasNorm = false;
    if (s1 == std::string::npos) return;

    auto s2 = tok.find('/', s1 + 1);
    if (s2 == std::string::npos) {
        // v/t
        ti = std::stoi(tok.substr(s1 + 1)) - 1;
        hasTex = true;
    } else if (s2 == s1 + 1) {
        // v//n
        ni = std::stoi(tok.substr(s2 + 1)) - 1;
        hasNorm = true;
    } else {
        // v/t/n
        ti = std::stoi(tok.substr(s1 + 1, s2 - s1 - 1)) - 1;
        ni = std::stoi(tok.substr(s2 + 1)) - 1;
        hasTex = hasNorm = true;
    }
}

static Scene parseScene(const std::string& filename) {
    std::ifstream f(filename);
    if (!f) { std::cerr << "Cannot open scene file: " << filename << "\n"; std::exit(1); }

    Scene sc;
    int curMat = -1, curTex = -1;
    std::string line;

    while (std::getline(f, line)) {
        auto tok = split(line);
        if (tok.empty()) continue;
        const auto& cmd = tok[0];

        if (cmd == "eye") {
            sc.eye = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
        } else if (cmd == "viewdir") {
            sc.viewDir = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
        } else if (cmd == "updir") {
            sc.upDir = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
        } else if (cmd == "vfov") {
            sc.vfov = std::stof(tok[1]);
        } else if (cmd == "imsize") {
            sc.width  = std::stoi(tok[1]);
            sc.height = std::stoi(tok[2]);
        } else if (cmd == "bkgcolor") {
            sc.bkgColor = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
        } else if (cmd == "mtlcolor") {
            Material m;
            m.diffuse  = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
            m.specular = { std::stof(tok[4]), std::stof(tok[5]), std::stof(tok[6]) };
            m.ka    = std::stof(tok[7]);
            m.kd    = std::stof(tok[8]);
            m.ks    = std::stof(tok[9]);
            m.n     = std::stof(tok[10]);
            m.alpha = std::stof(tok[11]);
            m.eta   = std::stof(tok[12]);
            sc.materials.push_back(m);
            curMat = (int)sc.materials.size() - 1;
        } else if (cmd == "texture") {
            sc.textures.push_back(loadTexturePPM(tok[1]));
            curTex = (int)sc.textures.size() - 1;
        } else if (cmd == "sphere") {
            Sphere sp;
            sp.center      = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
            sp.radius      = std::stof(tok[4]);
            sp.materialIdx = curMat;
            sp.textureIdx  = curTex;
            sc.spheres.push_back(sp);
        } else if (cmd == "light") {
            Light lt;
            lt.pos   = { std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) };
            lt.w     = std::stof(tok[4]);
            lt.color = { std::stof(tok[5]), std::stof(tok[6]), std::stof(tok[7]) };
            sc.lights.push_back(lt);
        } else if (cmd == "v") {
            sc.vertices.push_back({ std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) });
        } else if (cmd == "vn") {
            sc.normals.push_back({ std::stof(tok[1]), std::stof(tok[2]), std::stof(tok[3]) });
        } else if (cmd == "vt") {
            sc.uvCoords.push_back({ std::stof(tok[1]), std::stof(tok[2]), 0 });
        } else if (cmd == "f") {
            Face fc{};
            fc.materialIdx = curMat;
            fc.textureIdx  = curTex;
            bool ht1, ht2, ht3, hn1, hn2, hn3;
            parseFaceToken(tok[1], fc.vx, fc.tx, fc.nx, ht1, hn1);
            parseFaceToken(tok[2], fc.vy, fc.ty, fc.ny, ht2, hn2);
            parseFaceToken(tok[3], fc.vz, fc.tz, fc.nz, ht3, hn3);
            fc.hasTexture   = ht1;
            fc.smoothShading = hn1;
            sc.faces.push_back(fc);
        }
    }
    return sc;
}

// -----------------------------------------------------------------------------
// Intersection helpers
// -----------------------------------------------------------------------------
struct HitRecord {
    float  t;
    Vec3   point;
    Vec3   normal;
    int    materialIdx;
    int    textureIdx;
    Vec3   uv;          // for texture lookup
    bool   isFace;
    int    primitiveIdx;
};

static std::optional<float> raySphereIntersect(const Ray& ray, const Sphere& sp) {
    Vec3  oc = ray.origin - sp.center;
    float B  = 2.f * ray.dir.dot(oc);
    float C  = oc.lengthSq() - sp.radius * sp.radius;
    float D  = B*B - 4.f*C;
    if (D < 0) return std::nullopt;

    float sqrtD = std::sqrt(D);
    float t1 = (-B - sqrtD) / 2.f;
    float t2 = (-B + sqrtD) / 2.f;

    if (t1 > EPSILON) return t1;
    if (t2 > EPSILON) return t2;
    return std::nullopt;
}

// Möller–Trumbore — returns t or nullopt
static std::optional<std::pair<float,Vec3>> rayTriangleIntersect(
    const Ray& ray,
    const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    Vec3 e1 = p2 - p1, e2 = p3 - p1;
    Vec3 n  = e1.cross(e2);
    float d   = -n.dot(p1);
    float den = n.dot(ray.dir);
    if (std::abs(den) < 1e-9f) return std::nullopt;

    float t = -(n.dot(ray.origin) + d) / den;
    if (t <= EPSILON) return std::nullopt;

    Vec3  P  = ray.origin + ray.dir * t;
    float A1 = 0.5f * n.length();
    if (A1 < 1e-12f) return std::nullopt;

    Vec3 e3 = P - p2, e4 = P - p3;
    float a = 0.5f * e3.cross(e4).length() / A1;
    float b = 0.5f * e4.cross(e2).length() / A1;
    float g = 0.5f * e1.cross(e3).length() / A1;

    if (a < 0 || b < 0 || g < 0 || (a+b+g - 1.f) > EPSILON) return std::nullopt;

    Vec3 bary = { a, b, g };
    return std::make_pair(t, bary);
}

// -----------------------------------------------------------------------------
// Shadow attenuation factor S in [0,1]
// -----------------------------------------------------------------------------
static float shadowFactor(const Vec3& point, const Vec3& L, float lightDist,
                          const Scene& sc, int excludeSphere, int excludeFace) {
    float S = 1.f;
    Ray shadowRay{ point + L * (10.f * EPSILON), L };

    for (int i = 0; i < (int)sc.spheres.size(); ++i) {
        if (i == excludeSphere) continue;
        auto t = raySphereIntersect(shadowRay, sc.spheres[i]);
        if (t && *t < lightDist)
            S *= (1.f - sc.materials[sc.spheres[i].materialIdx].alpha);
    }

    for (int j = 0; j < (int)sc.faces.size(); ++j) {
        if (j == excludeFace) continue;
        const Face& fc = sc.faces[j];
        auto res = rayTriangleIntersect(shadowRay,
            sc.vertices[fc.vx], sc.vertices[fc.vy], sc.vertices[fc.vz]);
        if (res && res->first < lightDist)
            S *= (1.f - sc.materials[fc.materialIdx].alpha);
    }
    return S;
}

// -----------------------------------------------------------------------------
// Shading — Blinn-Phong with shadow
// -----------------------------------------------------------------------------
static Vec3 shade(const HitRecord& hit, const Scene& sc) {
    const Material& mat = sc.materials[hit.materialIdx];

    // Base diffuse color (material or texture)
    Vec3 diffColor = mat.diffuse;
    if (hit.textureIdx >= 0) {
        const Texture& tex = sc.textures[hit.textureIdx];
        int i = (int)std::round(hit.uv.x * (tex.w - 1));
        int j = (int)std::round(hit.uv.y * (tex.h - 1));
        diffColor = texelAt(tex, i, j);
    }

    Vec3 N = hit.normal.normalized();
    Vec3 V = (-sc.viewDir).normalized();
    Vec3 color = {};

    for (const auto& lt : sc.lights) {
        Vec3  L_raw  = lt.w == 0 ? Vec3{-lt.pos.x, -lt.pos.y, -lt.pos.z}
                                  : lt.pos - hit.point;
        float dist   = lt.w == 0 ? INF : L_raw.length();
        Vec3  L      = L_raw.normalized();
        Vec3  H      = (L + V).normalized();

        int exSphere = hit.isFace ? -1 : hit.primitiveIdx;
        int exFace   = hit.isFace ? hit.primitiveIdx : -1;
        float S = shadowFactor(hit.point, L, dist, sc, exSphere, exFace);

        float NL = std::max(0.f, N.dot(L));
        float NH = std::max(0.f, N.dot(H));

        color += diffColor * mat.ka;
        color += lt.color * (S * (mat.kd * diffColor * NL
                                + mat.ks * mat.specular * std::pow(NH, mat.n)));
    }
    return color;
}

// -----------------------------------------------------------------------------
// Closest hit in scene
// -----------------------------------------------------------------------------
static std::optional<HitRecord> closestHit(const Ray& ray, const Scene& sc) {
    std::optional<HitRecord> best;
    float minT = INF;

    // Spheres
    for (int i = 0; i < (int)sc.spheres.size(); ++i) {
        const Sphere& sp = sc.spheres[i];
        auto t = raySphereIntersect(ray, sp);
        if (!t || *t >= minT) continue;
        minT = *t;

        HitRecord rec;
        rec.t           = *t;
        rec.point       = ray.origin + ray.dir * *t;
        rec.normal      = (rec.point - sp.center) / sp.radius;
        rec.materialIdx = sp.materialIdx;
        rec.textureIdx  = sp.textureIdx;
        rec.isFace      = false;
        rec.primitiveIdx = i;

        // Spherical UV
        if (sp.textureIdx >= 0) {
            Vec3 d = (rec.point - sp.center) / sp.radius;
            float phi   = std::acos(-d.x / sp.radius);
            float theta = std::atan2(-d.z, d.y);
            if (theta < 0) theta += 2 * PI;
            rec.uv = { 0.5f * theta / PI, phi / PI, 0 };
        }
        best = rec;
    }

    // Faces
    for (int j = 0; j < (int)sc.faces.size(); ++j) {
        const Face& fc = sc.faces[j];
        const Vec3& p1 = sc.vertices[fc.vx];
        const Vec3& p2 = sc.vertices[fc.vy];
        const Vec3& p3 = sc.vertices[fc.vz];

        auto res = rayTriangleIntersect(ray, p1, p2, p3);
        if (!res || res->first >= minT) continue;
        minT = res->first;

        Vec3 bary = res->second;
        HitRecord rec;
        rec.t            = res->first;
        rec.point        = ray.origin + ray.dir * rec.t;
        rec.materialIdx  = fc.materialIdx;
        rec.textureIdx   = fc.textureIdx;
        rec.isFace       = true;
        rec.primitiveIdx = j;

        if (fc.smoothShading) {
            Vec3 n = sc.normals[fc.nx] * bary.x
                   + sc.normals[fc.ny] * bary.y
                   + sc.normals[fc.nz] * bary.z;
            rec.normal = n.normalized();
        } else {
            rec.normal = (p2 - p1).cross(p3 - p1).normalized();
        }

        if (fc.hasTexture) {
            Vec3 uv1 = sc.uvCoords[fc.tx];
            Vec3 uv2 = sc.uvCoords[fc.ty];
            Vec3 uv3 = sc.uvCoords[fc.tz];
            rec.uv = uv1 * bary.x + uv2 * bary.y + uv3 * bary.z;
        }
        best = rec;
    }

    return best;
}

// -----------------------------------------------------------------------------
// Recursive ray trace with Fresnel reflection/refraction
// -----------------------------------------------------------------------------
static Vec3 traceRay(const Ray& ray, const Scene& sc,
                     std::stack<float>& etaStack, int depth) {
    auto hit = closestHit(ray, sc);
    if (!hit) return sc.bkgColor;

    Vec3 color = shade(*hit, sc);
    if (depth >= MAX_DEPTH) return color;

    const Material& mat = sc.materials[hit->materialIdx];
    Vec3 N = hit->normal.normalized();
    Vec3 I = (-ray.dir).normalized();
    float NdotI = N.dot(I);

    // Determine which side of the surface we're on
    Vec3  Nf   = (NdotI >= 0) ? N : -N;
    float ni   = etaStack.top();
    float nt   = (NdotI >= 0) ? mat.eta : (etaStack.size() > 1
                                            ? *(&etaStack.top() - 1) : 1.f);
    float cosI = std::abs(NdotI);

    // Schlick Fresnel
    float F0 = std::pow((nt - ni) / (nt + ni), 2.f);
    float Fr = F0 + (1.f - F0) * std::pow(1.f - cosI, 5.f);

    // Reflection ray
    Vec3 rDir = (2.f * cosI * Nf - I).normalized();
    Ray  rRay{ hit->point + Nf * EPSILON, rDir };

    // Refraction (Snell's law)
    float sinT2 = (ni / nt) * (ni / nt) * (1.f - cosI * cosI);
    Vec3 refractColor = {};
    bool totalInternalReflection = sinT2 > 1.f;

    if (!totalInternalReflection) {
        float cosT = std::sqrt(1.f - sinT2);
        Vec3 tDir  = (ni / nt) * (-I) + (ni / nt * cosI - cosT) * Nf;
        Ray  tRay{ hit->point - Nf * EPSILON, tDir.normalized() };

        if (NdotI >= 0) etaStack.push(mat.eta);
        else if (etaStack.size() > 1) etaStack.pop();

        refractColor = traceRay(tRay, sc, etaStack, depth + 1);

        // restore stack
        if (NdotI >= 0) etaStack.pop();
        else etaStack.push(mat.eta);
    }

    Vec3 reflectColor = traceRay(rRay, sc, etaStack, depth + 1);

    color += Fr * reflectColor
           + (1.f - Fr) * (1.f - mat.alpha) * refractColor;

    return color;
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <scene.txt>\n";
        return 1;
    }

    fs::path inputPath(argv[1]);
    Scene sc = parseScene(argv[1]);

    // Camera basis
    Vec3 viewN = sc.viewDir.normalized();
    Vec3 u     = viewN.cross(sc.upDir).normalized();
    if (std::isnan(u.x)) {
        std::cerr << "View and up directions are parallel.\n"; return 1;
    }
    Vec3 v = u.cross(viewN).normalized();

    float W    = (float)sc.width;
    float H    = (float)sc.height;
    float dist = H / (2.f * std::tan(sc.vfov * PI / 360.f));

    Vec3 eye = sc.eye;
    Vec3 UL  = eye + viewN * dist - u * (W/2) + v * (H/2);
    Vec3 dh  = u * (1.f);   // per-pixel horizontal step
    Vec3 dv  = v * (-1.f);  // per-pixel vertical step  (y flipped)

    // Output PPM
    fs::path outPath = inputPath.stem();
    outPath += ".ppm";

    std::ofstream out(outPath);
    out << "P3\n" << sc.width << " " << sc.height << "\n255\n";

    std::cout << "Rendering " << sc.width << "x" << sc.height << "...\n";

    for (int row = 0; row < sc.height; ++row) {
        for (int col = 0; col < sc.width; ++col) {
            Vec3 pixel = UL + dh * (col + 0.5f) + dv * (row + 0.5f);
            Ray  ray{ eye, (pixel - eye).normalized() };

            std::stack<float> etaStack;
            etaStack.push(1.f);  // start in air

            Vec3 color = traceRay(ray, sc, etaStack, 0);

            int r = std::clamp((int)(color.x * 255), 0, 255);
            int g = std::clamp((int)(color.y * 255), 0, 255);
            int b = std::clamp((int)(color.z * 255), 0, 255);
            out << r << " " << g << " " << b << " ";
        }
        out << "\n";
    }

    std::cout << "Done. Output: " << outPath << "\n";
    return 0;
}

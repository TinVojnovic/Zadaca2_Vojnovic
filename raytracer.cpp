#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "algorithm"
#include "geometry.h"

using namespace std;

typedef std::vector<Vec3f> Image;

struct Light
{
    Vec3f pos;
    float intensity;

    Light(const Vec3f &position, const float &intensity) : pos(position), intensity(intensity) {}
};
typedef std::vector<Light> Lights;

struct Material
{
    Vec2f albedo; // difuzni i spekularni koeficijenti refleksije
    Vec3f color;
    float spec_exp;

    Material(const Vec2f &a, const Vec3f &color, const float &coef) : albedo(a), color(color), spec_exp(coef) {}
    Material() : albedo(Vec2f(1, 0)), color(), spec_exp(1.f) {}
};

struct Object
{
    Material mat;
    virtual bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t) const = 0;
    virtual Vec3f normal(const Vec3f &p) const = 0;
};
typedef std::vector<Object *> Objects;

struct Sphere : Object
{
    Vec3f c; // centar
    float r; // radius

    Sphere(const Vec3f &c, const float &r, const Material &m) : c(c), r(r)
    {
        Object::mat = m;
    }

    bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t) const
    {
        Vec3f v = c - p; // vektor izmedju izvora zrake i centra

        if (d * v < 0) // skalarni produkt
        {
            // sfera se nalazi iza zrake
            return false;
        }
        else
        {
            // izracunaj projekciju
            Vec3f pc = p + d * ((d * v) / d.norm());
            if ((c - pc) * (c - pc) > r * r)
            {
                // nema sjeciste
                return false;
            }
            else
            {
                float dist = sqrt(r * r - (c - pc) * (c - pc));

                if (v * v > r * r) // izvor pogleda izvan sfere
                {
                    t = (pc - p).norm() - dist;
                }
                else // izvor pogleda unutar sfere
                {
                    t = (pc - p).norm() + dist;
                }

                return true;
            }
        }
    }

    Vec3f normal(const Vec3f &p) const
    {
        return (p - c).normalize();
    }
};

struct Cuboid : Object {
  Vec3f s, e;

  Cuboid(const Vec3f &s, const Vec3f &e, const Material &m) : s(s), e(e) {
    Object::mat = m;
  }

  bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t) const {
    float ts = numeric_limits<float>::min(), tb = numeric_limits<float>::max();
    float minX = min(s[0], e[0]),
          minY = min(s[1], e[1]),
          minZ = min(s[2], e[2]),
          maxX = max(s[0], e[0]),
          maxY = max(s[1], e[1]),
          maxZ = max(s[2], e[2]);

        if (d[0] == 0 && (p[0] < minX || p[0] > maxX)) return false;
        else {
          float t1 = (minX - p[0]) / d[0], t2 = (maxX - p[0]) / d[0];
          if (t1 > t2) swap(t1, t2); 
          ts = max(ts, t1);
          tb = min(tb, t2);
          if (ts > tb || tb < 0) return false;
        }
        t = ts;
        if (d[1] == 0 && (p[1] < minY || p[1] > maxY)) return false;
        else {
          float t1 = (minY - p[1]) / d[1], t2 = (maxY - p[1]) / d[1];
          if (t1 > t2) swap(t1, t2);
          ts = max(ts, t1);
          tb = min(tb, t2);
          if (ts > tb || tb < 0) return false;
        }
        t = ts;
        if (d[2] == 0 && (p[2] < minZ || p[2] > maxZ)) return false;
        else {
            float t1 = (minZ - p[2]) / d[2], t2 = (maxZ - p[2]) / d[2];
            if (t1 > t2) swap(t1, t2);
            ts = max(ts, t1);
            tb = min(tb, t2);
            if (ts > tb || tb < 0) return false;
        }
        t = ts;
        return true;
  }

  Vec3f normal(const Vec3f &p) const {
    if(abs(p[0] - s[0]) < 0.0001) return Vec3f(-1,0,0);
    else if(abs(p[0] - e[0]) < 0.0001) return Vec3f(1,0,0);
    else if(abs(p[1] - s[1]) < 0.0001) return Vec3f(0,-1,0);
    else if(abs(p[1] - e[1]) < 0.0001) return Vec3f(0,1,0);
    else if(abs(p[2] - s[2]) < 0.0001) return Vec3f(0,0,-1);
    else if(abs(p[2] - e[2]) < 0.0001) return Vec3f(0,0,1);
  }
};

struct Cylinder : Object {
  Vec3f c;
  float r;
  float h;

  Cylinder(const Vec3f &c, const float &r, const float &h, const Material &m) : c(c), r(r), h(h) {
    Object::mat = m;
  }

  bool ray_intersect(const Vec3f &p, const Vec3f &d, float &t) const {
    if((c - p)*d < 0) return false;
    else {
      float A = (d[0]*d[0])+(d[2]*d[2]);
      float B = 2*(d[0]*(p[0]-c[0])+d[2]*(p[2]-c[2]));
      float C = (p[0]-c[0])*(p[0]-c[0])+(p[2]-c[2])*(p[2]-c[2])-(r*r);
      float D = B*B - 4*(A*C);
      if(D < 0) return false;
      else {
        float t1 = min((-B - sqrt(D))/(2*A), (-B + sqrt(D))/(2*A));
        float t2 = max((-B - sqrt(D))/(2*A), (-B + sqrt(D))/(2*A));
        if((p[1] + t1*d[1] >= c[1]) && (p[1] + t1*d[1] <= c[1] + h)) {t = t1; return true;}
        else if ((p[1] + t2*d[1] >= c[1]) && (p[1] + t2*d[1] <= c[1] + h)) {t = t2; return true;}
        else return false;
      }
    }
  }

  Vec3f normal(const Vec3f &p) const {
    Vec3f n = (p-c).normalize();
    n[1] = 0;
    return n;
  }
};

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const Objects &objs, Vec3f &hit, Material &material, Vec3f &N)
{
    // inicijalno, pretp. da je predmet daleko
    float dist = std::numeric_limits<float>::max();
    float obj_dist = dist;

    for (auto obj : objs)
    {
        if (obj->ray_intersect(orig, dir, obj_dist) && obj_dist < dist)
        {
            dist = obj_dist;
            hit = orig + dir * obj_dist;
            N = obj->normal(hit);
            material = obj->mat;
        }
    }

    // provjeri je li sjecište predaleko
    return dist < 1000;
}

Vec3f reflect(const Vec3f &I, const Vec3f &N)
{
    return I - (N * (2 * (I * N)));
}

// funkcija koja vraca udaljenost sjecista pravca i sfere
Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const Objects &objs, const Lights &lights, const float &depth)
{
    int maxDepth = 6;
    Vec3f hit_point;
    Vec3f hit_normal; // normala na povrsinu
    Material hit_material;

    if (!scene_intersect(orig, dir, objs, hit_point, hit_material, hit_normal) || depth > maxDepth)
    {
        return Vec3f(0.5, 0.1, 0.5); // vrati boju pozadine
    }
    else
    {
        float diffuse_light_intensity = 0;
        float specular_light_intensity = 0;

        for (auto light : lights)
        {
            Vec3f light_dir = (light.pos - hit_point).normalize(); // smjer svjetla
            float light_dist = (light.pos - hit_point).norm();     // udaljenost do svjetla

            // sjene
            Vec3f shadow_orig;
            Vec3f shadow_hit_point;
            Vec3f shadow_hit_normal;
            Material shadow_material;

            // epsilon pomak od tocke sjecista
            if (light_dir * hit_normal < 0)
            {
                shadow_orig = hit_point - hit_normal * 0.001;
                hit_normal = -hit_normal;
            }
            else
            {
                shadow_orig = hit_point + hit_normal * 0.001;
            }

            if (scene_intersect(shadow_orig, light_dir, objs, shadow_hit_point, shadow_material, shadow_hit_normal) && (shadow_hit_point - shadow_orig).norm() < light_dist)
            {
                continue;
            }

            // sjencanje
            // lambertov model
            diffuse_light_intensity += light.intensity * std::max(0.f, light_dir * hit_normal);

            // blinn-phongov model
            // smjer pogleda
            Vec3f view_dir = (orig - hit_point).normalize();
            // poluvektor
            Vec3f half_vec = (view_dir + light_dir).normalize();
            specular_light_intensity += light.intensity * powf(std::max(0.f, half_vec * hit_normal), hit_material.spec_exp);
        }

        Vec3f R = reflect(dir, hit_normal);

        Vec3f hitColor = hit_material.color * hit_material.albedo[0] * diffuse_light_intensity // diffuse dio
                         + Vec3f(1, 1, 1) * hit_material.albedo[1] * specular_light_intensity;         // specular dio

        return hitColor = hitColor + cast_ray(hit_point + hit_normal * 0.1, R, objs, lights, depth + 1) * 0.3;
        //return hitColor;
    }
}

// funkcija za crtanje
void render(const Objects &objects, const Lights &lights)
{
    // velicina slike
    const int width = 1024;
    const int height = 768;
    const int fov = 3.141592 / 2.0; // pi / 2

    // spremnik za sliku
    Image buffer(width * height);

    // nacrtaj u sliku
    for (size_t j = 0; j < height; ++j)
    {
        for (size_t i = 0; i < width; ++i)
        {
            // pošalji zraku u svaki piksel
            float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);

            // definiraj smjer
            Vec3f dir = Vec3f(x, y, -1).normalize();

            buffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), dir, objects, lights, 0);
        }
    }

    std::ofstream ofs;
    ofs.open("./scena.ppm", std::ofstream::binary);
    // oblikuj po ppm formatu
    ofs << "P6\n"
        << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            // skaliraj na sa [0, 1] na [0, 255]
            unsigned char color = (unsigned char)(255.f * std::max(0.f, std::min(1.f, buffer[i][j])));
            ofs << color;
        }
    }

    ofs.close();
}

int main()
{
    Material red = Material(Vec2f(0.6, 0.3), Vec3f(1, 0, 0), 100);
    Material green = Material(Vec2f(0.6, 0.3), Vec3f(0, 0.5, 0), 60);
    Material blue = Material(Vec2f(0.9, 0.1), Vec3f(0, 0, 1), 10);
    Material l_blue = Material(Vec2f(0.9, 0.1), Vec3f(0.5, 0.5, 1), 10);
    Material dark_gray = Material(Vec2f(0.9, 0.1), Vec3f(0.3, 0.3, 0.3), 10);

    Cuboid surface(Vec3f(-25, -5.1, -30), Vec3f(25, -5, -9), dark_gray);

    Sphere sphere1(Vec3f(-1.0, -3.5, -15), 1.5, green);
    Sphere sphere2(Vec3f(3, -4.5, -11.5), 0.5, blue);
    Cuboid cuboid1(Vec3f(5, -4, -17), Vec3f(3, -1, -15), red);
    Cylinder cylinder(Vec3f(-6, -2, -13), 1, 3, l_blue);

    Objects objs = {&surface, &sphere1, &sphere2, &cuboid1, &cylinder};

    Light l1 = Light(Vec3f(-20, 20, 20), 1);
    Light l2 = Light(Vec3f(20, 30, 20), 1.5);
    Lights lights = {l1, l2};

    render(objs, lights);

    return 0;
}
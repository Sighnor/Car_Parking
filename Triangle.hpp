#pragma once

#include "BVH.hpp"
#include "Intersection.hpp"
#include "Material.hpp"
#include "OBJ_Loader.hpp"
#include "Object.hpp"
#include "Triangle.hpp"
#include <cassert>
#include <array>

bool rayTriangleIntersect(const Vector3f& v0, const Vector3f& v1,
                          const Vector3f& v2, const Vector3f& orig,
                          const Vector3f& dir, float& tnear, float& u, float& v)
{
    Vector3f edge1 = v1 - v0;
    Vector3f edge2 = v2 - v0;
    Vector3f pvec = crossProduct(dir, edge2);
    float det = dotProduct(edge1, pvec);
    if (det == 0 || det < 0)
        return false;

    Vector3f tvec = orig - v0;
    u = dotProduct(tvec, pvec);
    if (u < 0 || u > det)
        return false;

    Vector3f qvec = crossProduct(tvec, edge1);
    v = dotProduct(dir, qvec);
    if (v < 0 || u + v > det)
        return false;

    float invDet = 1 / det;

    tnear = dotProduct(edge2, qvec) * invDet;
    u *= invDet;
    v *= invDet;

    return true;
}

class Triangle : public Object
{
public:
    Vector3f v0, v1, v2; // vertices A, B ,C , counter-clockwise order
    Vector3f e1, e2;     // 2 edges v1-v0, v2-v0;
    Vector2f t0, t1, t2; // texture coords
    Vector3f normal;
    Vector3f normal0, normal1, normal2;
    Material* m;

    Triangle(Vector3f _v0, Vector3f _v1, Vector3f _v2, Vector3f _n0, Vector3f _n1, Vector3f _n2, Vector2f _t0, Vector2f _t1 , Vector2f _t2 ,Material* _m = nullptr)
        : v0(_v0), v1(_v1), v2(_v2), normal0(_n0), normal1(_n1), normal2(_n2), t0(_t0), t1(_t1), t2(_t2), m(_m)
    {
        e1 = v1 - v0;
        e2 = v2 - v0;
        normal = normalize(crossProduct(e1, e2));
    }

    bool intersect(const Ray& ray) override;
    bool intersect(const Ray& ray, float& tnear,
                   uint32_t& index) const override;
    Intersection getIntersection(Ray ray) override;
    void getSurfaceProperties(const Vector3f& P, const Vector3f& I,
        const uint32_t& index, const Vector2f& uv,
        Vector3f& N, Vector2f& st) const override
    {
        N = normal1;
        //        throw std::runtime_error("triangle::getSurfaceProperties not
        //        implemented.");
    }
    Vector3f evalDiffuseColor(const Vector2f&) const override;
    Bounds3 getBounds() override;
};

class MeshTriangle : public Object
{
public:
    MeshTriangle(const std::string& filename)
    {
        objl::Loader loader;
        loader.LoadFile(filename);

        //assert(loader.LoadedMeshes.size() == 7);
        auto mesh = loader.LoadedMeshes;

        Vector3f min_vert = Vector3f{std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity()};
        Vector3f max_vert = Vector3f{-std::numeric_limits<float>::infinity(),
                                     -std::numeric_limits<float>::infinity(),
                                     -std::numeric_limits<float>::infinity()};

        Vector3f axis(0,0,-1);
    	float rotation_angle = 0.0;

        float c = std::cos(rotation_angle/180*M_PI);
    	float k =1-c;
    	float s = std::sin(rotation_angle/180*M_PI);

    	float kxx = k * axis.x * axis.x +  c;
   	float kxy = k * axis.x * axis.y - axis.z * s;
   	float kxz = k * axis.x * axis.z +  axis.y * s;
    	float kyx = k * axis.x * axis.y +  axis.z * s;
    	float kyy = k * axis.y * axis.y +  c;
    	float kyz = k * axis.y * axis.z - axis.x * s;
    	float kzx = k * axis.x * axis.z - axis.y * s;
    	float kzy = k * axis.y * axis.z +  axis.x * s;
    	float kzz = k * axis.z * axis.z +  c;
        for (int k = 0;k < loader.LoadedMeshes.size(); ++k) {
             std::cout<< mesh[k].MeshName<<'\n';
             auto new_mat =
                    new Material(MaterialType::DIFFUSE_AND_GLOSSY,
                        Vector3f(0.5, 0.5, 0.5), Vector3f(0, 0, 0));
                new_mat->Kd = mesh[k].MeshMaterial.Kd.X;
                new_mat->map_Kd = mesh[k].MeshMaterial.map_Kd;
                new_mat->Ks = mesh[k].MeshMaterial.Ks.X;
                new_mat->specularExponent = mesh[k].MeshMaterial.Ns;
            for (int i = 0; i < mesh[k].Vertices.size(); i += 3) {
                std::array<Vector3f, 3> face_vertices;
                std::array<Vector3f, 3> face_normals;
                std::array<Vector2f, 3> face_textures;
                for (int j = 0; j < 3; j++) {
                    auto vert = Vector3f(mesh[k].Vertices[i + j].Position.X,
                        mesh[k].Vertices[i + j].Position.Y,
                        mesh[k].Vertices[i + j].Position.Z) *
                        60.f;
                    face_vertices[j] = vert;

                    face_normals[j].x = mesh[k].Vertices[i + j].Normal.X;
                    face_normals[j].y = mesh[k].Vertices[i + j].Normal.Y;
                    face_normals[j].z = mesh[k].Vertices[i + j].Normal.Z;

                    min_vert = Vector3f(std::min(min_vert.x, vert.x),
                        std::min(min_vert.y, vert.y),
                        std::min(min_vert.z, vert.z));
                    max_vert = Vector3f(std::max(max_vert.x, vert.x),
                        std::max(max_vert.y, vert.y),
                        std::max(max_vert.z, vert.z));
                    
                    face_textures[j].x = mesh[k].Vertices[i + j].TextureCoordinate.X;
                    face_textures[j].y = mesh[k].Vertices[i + j].TextureCoordinate.Y;
                }

                triangles.emplace_back(face_vertices[0], face_vertices[1],
                    face_vertices[2], face_normals[0], face_normals[1], face_normals[2], face_textures[0], face_textures[1], face_textures[2] ,new_mat);
            }
        }

        bounding_box = Bounds3(min_vert, max_vert);

        std::vector<Object*> ptrs;
        for (auto& tri : triangles)
            ptrs.push_back(&tri);

        bvh = new BVHAccel(ptrs);
    }
    
    MeshTriangle(Vector3f pMin, Vector3f pMax, int num[3], Material* _m, Vector2f tempuv)
    {

        bounding_box = Bounds3(pMin, pMax);

        float perx = (pMax.x-pMin.x)/num[0];
        float pery = (pMax.y-pMin.y)/num[1];
        float perz = (pMax.z-pMin.z)/num[2];
        float perux = 1.0f/num[0];
        float peruy = 1.0f/num[1];
        float peruz = 1.0f/num[2];
        int numx = num[0]+1;
        int numy = num[1]+1;
        int numz = num[2]+1;

        auto m0 =new Material(MaterialType::DIFFUSE_AND_GLOSSY, Vector3f(0.5, 0.5, 0.5), Vector3f(0, 0, 0));//材质信息
        m0->Kd = 0.7;
        m0->map_Kd = "P0.jpg";       
        m0->Ks = 0.15;
        m0->specularExponent = 50;

        //BUILDXY
        Vector3f face_vertices[numx][numy];
        Vector2f face_textures[2][2];
        Vector3f face_normal;
        face_normal.x = 0;
        face_normal.y = 0;
        face_normal.z = 1.f;

        face_textures[0][0].x = tempuv.x;
        face_textures[0][0].y = 1 + tempuv.y;
        face_textures[1][0].x = 1 - tempuv.x;
        face_textures[1][0].y = 1 + tempuv.y;
        face_textures[0][1].x = tempuv.x;
        face_textures[0][1].y = 2 - tempuv.y;
        face_textures[1][1].x = 1 - tempuv.x;
        face_textures[1][1].y = 2 - tempuv.y;
        
       
        for (int k = 0;k < numx; ++k) {
            for (int i = 0; i < numy; ++i) {
                    face_vertices[k][i].x = pMin.x + perx * k;
                    face_vertices[k][i].y = pMin.y + pery * i;
                    face_vertices[k][i].z = pMax.z;
            }
        }
   
       for (int k = 0;k < num[0]; ++k) {
            for (int i = 0; i < num[1]; ++i) {
              triangles.emplace_back(face_vertices[k][i+1], face_vertices[k][i],face_vertices[k+1][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[0][0], face_textures[1][0],
                                     m0);
              triangles.emplace_back(face_vertices[k][i+1], face_vertices[k+1][i],face_vertices[k+1][i+1],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][0], face_textures[1][1],
                                     m0);
            }
        }


      face_normal.z = -1.f;
      for (int k = 0;k < numx; ++k) {
            for (int i = 0; i < numy; ++i) {
                    face_vertices[k][i].z = pMin.z;
            }
        }

      for (int k = 0;k < num[0]; ++k) {
            for (int i = 0; i < num[1]; ++i) {
              triangles.emplace_back(face_vertices[k][i+1], face_vertices[k+1][i],face_vertices[k][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][0], face_textures[0][0],
                                     m0);
              triangles.emplace_back(face_vertices[k][i+1], face_vertices[k+1][i+1],face_vertices[k+1][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][1], face_textures[1][0],
                                     m0);
            }
        }

 
      //BUILDXZ
        Vector3f xzface_vertices[numz][numx];
        face_normal.x = 0;
        face_normal.y = 1.f;
        face_normal.z = 0;

        for (int k = 0;k < numz; ++k) {
            for (int i = 0; i < numx; ++i) {
                    xzface_vertices[k][i].x = pMin.x + perx * i;
                    xzface_vertices[k][i].y = pMax.y;
                    xzface_vertices[k][i].z = pMin.z + perz * k;
            }
        }
   
       for (int k = 0;k < num[2]; ++k) {
            for (int i = 0; i < num[0]; ++i) {
              triangles.emplace_back(xzface_vertices[k][i+1], xzface_vertices[k][i],xzface_vertices[k+1][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[0][0], face_textures[1][0],
                                     _m);
              triangles.emplace_back(xzface_vertices[k][i+1], xzface_vertices[k+1][i],xzface_vertices[k+1][i+1],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][0], face_textures[1][1],
                                     _m);
            }
        }


      face_normal.y = -1.f;
      for (int k = 0;k < numz; ++k) {
            for (int i = 0; i < numx; ++i) {
                    xzface_vertices[k][i].y = pMin.y;
            }
        }

      for (int k = 0;k < num[2]; ++k) {
            for (int i = 0; i < num[0]; ++i) {
              triangles.emplace_back(xzface_vertices[k][i+1], xzface_vertices[k+1][i],xzface_vertices[k][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][0], face_textures[0][0],
                                     _m);
              triangles.emplace_back(xzface_vertices[k][i+1], xzface_vertices[k+1][i+1],xzface_vertices[k+1][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][1], face_textures[1][0],
                                     _m);
            }
        }


        //BUILDYZ
        Vector3f yzface_vertices[numy][numz];
        face_normal.x = 1.f;
        face_normal.y = 0;
        face_normal.z = 0;

        for (int k = 0;k < numy; ++k) {
            for (int i = 0; i < numz; ++i) {
                    yzface_vertices[k][i].x = pMax.x;
                    yzface_vertices[k][i].y = pMin.y + pery * k;
                    yzface_vertices[k][i].z = pMin.z + perz * i;
            }
        }
   
       for (int k = 0;k < num[1]; ++k) {
            for (int i = 0; i < num[2]; ++i) {
              triangles.emplace_back(yzface_vertices[k][i+1], yzface_vertices[k][i],yzface_vertices[k+1][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[0][0], face_textures[1][0],
                                     m0);
              triangles.emplace_back(yzface_vertices[k][i+1], yzface_vertices[k+1][i],yzface_vertices[k+1][i+1],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][0], face_textures[1][1],
                                     m0);
            }
        }


      face_normal.x = -1.f;
      for (int k = 0;k < numy; ++k) {
            for (int i = 0; i < numz; ++i) {
                    yzface_vertices[k][i].x = pMin.x;
            }
        }

      for (int k = 0;k < num[1]; ++k) {
            for (int i = 0; i < num[2]; ++i) {
              triangles.emplace_back(yzface_vertices[k][i+1], yzface_vertices[k+1][i],yzface_vertices[k][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][0], face_textures[0][0],
                                     m0);
              triangles.emplace_back(yzface_vertices[k][i+1], yzface_vertices[k+1][i+1],yzface_vertices[k+1][i],
                                     face_normal, face_normal, face_normal, 
                                     face_textures[0][1], face_textures[1][1], face_textures[1][0],
                                     m0);
            }
        }


        std::vector<Object*> ptrs;
        for (auto& tri : triangles)
            ptrs.push_back(&tri);

        bvh = new BVHAccel(ptrs);
    }


    bool intersect(const Ray& ray) { return true; }

    bool intersect(const Ray& ray, float& tnear, uint32_t& index) const
    {
        bool intersect = false;
        for (uint32_t k = 0; k < numTriangles; ++k) {
            const Vector3f& v0 = vertices[vertexIndex[k * 3]];
            const Vector3f& v1 = vertices[vertexIndex[k * 3 + 1]];
            const Vector3f& v2 = vertices[vertexIndex[k * 3 + 2]];
            float t, u, v;
            if (rayTriangleIntersect(v0, v1, v2, ray.origin, ray.direction, t,
                                     u, v) &&
                t < tnear) {
                tnear = t;
                index = k;
                intersect |= true;
            }
        }

        return intersect;
    }

    Bounds3 getBounds() { return bounding_box; }

    void getSurfaceProperties(const Vector3f& P, const Vector3f& I,
                              const uint32_t& index, const Vector2f& uv,
                              Vector3f& N, Vector2f& st) const
    {
        const Vector3f& v0 = vertices[vertexIndex[index * 3]];
        const Vector3f& v1 = vertices[vertexIndex[index * 3 + 1]];
        const Vector3f& v2 = vertices[vertexIndex[index * 3 + 2]];
        Vector3f e0 = normalize(v1 - v0);
        Vector3f e1 = normalize(v2 - v1);
        N = normalize(crossProduct(e0, e1));
        const Vector2f& st0 = stCoordinates[vertexIndex[index * 3]];
        const Vector2f& st1 = stCoordinates[vertexIndex[index * 3 + 1]];
        const Vector2f& st2 = stCoordinates[vertexIndex[index * 3 + 2]];
        st = st0 * (1 - uv.x - uv.y) + st1 * uv.x + st2 * uv.y;
    }

    Vector3f evalDiffuseColor(const Vector2f& st) const
    {
        float scale = 5;
        float pattern =
            (fmodf(st.x * scale, 1) > 0.5) ^ (fmodf(st.y * scale, 1) > 0.5);
        return lerp(Vector3f(0.815, 0.235, 0.031),
                    Vector3f(0.937, 0.937, 0.231), pattern);
    }

    Intersection getIntersection(Ray ray)
    {
        Intersection intersec;

        if (bvh) {
            intersec = bvh->Intersect(ray);
        }

        return intersec;
    }

    Bounds3 bounding_box;
    std::unique_ptr<Vector3f[]> vertices;
    uint32_t numTriangles;
    std::unique_ptr<uint32_t[]> vertexIndex;
    std::unique_ptr<Vector2f[]> stCoordinates;

    std::vector<Triangle> triangles;

    BVHAccel* bvh;

    Material* m;
};

inline bool Triangle::intersect(const Ray& ray) { return true; }
inline bool Triangle::intersect(const Ray& ray, float& tnear,
                                uint32_t& index) const
{
    return false;
}

inline Bounds3 Triangle::getBounds() { return Union(Bounds3(v0, v1), v2); }

inline Intersection Triangle::getIntersection(Ray ray)
{
    Intersection inter;
    if (dotProduct(ray.direction, normal) > 0)
        return inter;
    double u, v, t_tmp = 0;
    Vector3f pvec = crossProduct(ray.direction, e2);
    double det = dotProduct(e1, pvec);
    if (fabs(det) < EPSILON)
        return inter;

    double det_inv = 1. / det;
    Vector3f tvec = ray.origin - v0;
    u = dotProduct(tvec, pvec) * det_inv;
    if (u < 0 || u > 1)
        return inter;
    Vector3f qvec = crossProduct(tvec, e1);
    v = dotProduct(ray.direction, qvec) * det_inv;
    if (v < 0 || u + v > 1)
        return inter;
    t_tmp = dotProduct(e2, qvec) * det_inv;

    // TODO find ray triangle intersection
    if (t_tmp < 0)
        return inter;
    inter.happened = true;
    inter.obj = this;
    inter.m = m;
    //std::cout<<m->map_Kd<<'\n';
    inter.distance = t_tmp;
    inter.coords = ray(t_tmp);
    inter.normal = normal;
    inter.normal0 = normal0;
    inter.normal1 = normal1;
    inter.normal2 = normal2;
    inter.uv0 = t0;
    inter.uv1 = t1;
    inter.uv2 = t2;
    inter.uv.x = (float)u;
    inter.uv.y = (float)v;

    return inter;
}

inline Vector3f Triangle::evalDiffuseColor(const Vector2f& st) const
{
    float scale = 5;
    float pattern =
        (fmodf(st.x * scale, 1) > 0.5) ^ (fmodf(st.y * scale, 1) > 0.5);
    return lerp(Vector3f(0.815, 0.235, 0.031),
        Vector3f(0.937, 0.937, 0.231), pattern);
    //return Vector3f(0.5, 0.5, 0.5);
}

//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
#include "Vector.hpp"
#include "Material.hpp"
class Object;
class Sphere;

struct Intersection
{
    Intersection(){
        happened=false;
        coords=Vector3f();
        normal=Vector3f();
        normal0=Vector3f();
        normal1=Vector3f();
        normal2=Vector3f();
        uv0=Vector2f();
        uv1=Vector2f();
        uv2=Vector2f();
        uv=Vector2f();
        distance= std::numeric_limits<double>::max();
        obj =nullptr;
        m=nullptr;
    }
    bool happened;
    Vector3f coords;
    Vector3f normal;
    Vector3f normal0;
    Vector3f normal1;
    Vector3f normal2;
    Vector2f uv0;
    Vector2f uv1;
    Vector2f uv2;
    Vector2f uv;
    double distance;
    Object* obj;
    Material* m;
};
#endif //RAYTRACING_INTERSECTION_H

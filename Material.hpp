//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"
#include <string>

enum MaterialType { DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION };

class Material{
private:
float D_GGX_TR(const Vector3f &N, const Vector3f &H, const float &temp_roughness)
{
    float a = temp_roughness * temp_roughness;
    float a2 = a * a;
    float NdotH = std::max(dotProduct(N, H), 0.0f);
    float NdotH2 = NdotH * NdotH;

    float nom = a2;
    float denom  = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = M_PI * denom * denom;

    return nom / denom;
}

float GeometrySchlickGGX(const float &NdotV, const float &temp_roughness)
{
    float r = (temp_roughness + 1.0);
    float k = (r * r) / 8.0; 
    float nom = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}

float GeometrySmith(const Vector3f &N, const Vector3f &wr, const Vector3f &wo, const float &temp_roughness)
{
    float NdotV = std::max(dotProduct(N, wr), 0.0f);
    float NdotL = std::max(dotProduct(N, wo), 0.0f);
    float ggx1 = GeometrySchlickGGX(NdotV, temp_roughness);
    float ggx2 = GeometrySchlickGGX(NdotL, temp_roughness);

    return ggx1 * ggx2;
}

Vector3f fresnelSchlick(const float &cosTheta, const Vector3f &F0)
{
    return F0 + (Vector3f(1.f) - F0) * std::pow(1.0 - cosTheta, 5.0);
}

public:
    MaterialType m_type;
    Vector3f m_color;
    Vector3f m_emission;
    Vector3f F0;
    float ior;
    float Kd, Ks;
    float specularExponent;
    float metalness;
    float roughness;
    bool has_map;
    std::string map_Kd;
    std::string map_reflect;
    std::string map_bump;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE_AND_GLOSSY, Vector3f c=Vector3f(1,1,1), Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline Vector3f eval(const Vector3f &wr, const Vector3f &wo, const Vector3f &N, const Vector3f &my_kd, const float &temp_roughness, const float &temp_Kd);

};

Material::Material(MaterialType t, Vector3f c, Vector3f e){
    m_type = t;
    m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}

Vector3f Material::eval(const Vector3f &wr, const Vector3f &wo, const Vector3f &N, const Vector3f &my_kd, const float &temp_roughness, const float &temp_Kd){
    // calculate the contribution of diffuse   model
    float cosalpha = dotProduct(N, wo);
    if (cosalpha > 0.0f) {
        Vector3f h = normalize(wo - wr);
        float costheta = dotProduct(N, -wr);
        Vector3f my_ks = fresnelSchlick(costheta, F0);
        Vector3f diffuse = (Vector3f(1.f) - my_ks) * my_kd / M_PI;
        Vector3f microfacet = my_ks * GeometrySmith(N, -wr, wo, temp_roughness) * D_GGX_TR(N, h, temp_roughness) / std::max(4.f * cosalpha * costheta, 0.001f);
        return temp_Kd * (diffuse + microfacet);
    }
    else
        return Vector3f(0.0f);
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}
#endif //RAYTRACING_MATERIAL_H

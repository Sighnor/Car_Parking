//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include "Texture.hpp"

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of the Whitted-syle light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected object is either reflective or reflective and refractive,
// then we compute the reflection/refracton direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refractin depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is duffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.
//Vector3f Scene::castRay(const Ray &ray, int depth, Texture tex1, Texture tex2, Texture tex3, Texture tex4, Texture tex5, Texture tex6, Texture tex7 ,Texture tex8, Texture tex9, Texture tex10) const
Vector3f Scene::castRay(const Ray &ray, int depth, Texture tex[]) const
{
    if (depth > this->maxDepth) {
        return Vector3f(0.0,0.0,0.0);
    }
    Intersection intersection = Scene::intersect(ray);
    Material *m = intersection.m;
    Object *hitObject = intersection.obj;
    Vector3f hitColor = this->backgroundColor;
//    float tnear = kInfinity;
    Vector2f uv = intersection.uv;
    uint32_t index = 0;
    if(intersection.happened) {

        Vector3f hitPoint = intersection.coords;
        Vector3f N = intersection.normal; // normal
        Vector2f st; // st coordinates
        hitObject->getSurfaceProperties(hitPoint, ray.direction, index, uv, N, st);
        N =  (1 - uv.x - uv.y) * intersection.normal0 + uv.x * intersection.normal1 + uv.y * intersection.normal2;
//        Vector3f tmp = hitPoint;


                st.x = (1 - uv.x - uv.y) * intersection.uv0.x + uv.x * intersection.uv1.x + uv.y * intersection.uv2.x;
                st.y = (1 - uv.x - uv.y) * intersection.uv0.y + uv.x * intersection.uv1.y + uv.y * intersection.uv2.y - 1;
                std::string texture_path = m->map_Kd;
                Vector3f mycolor;

                //float eye_distance = sqrt(pow(ray.origin.x - hitPoint.x, 2) + pow(ray.origin.z - hitPoint.z, 2));
                float eye_distance = dotProduct(hitPoint - ray.origin, ray.direction);

                float roughness = m->roughness;
                float Kd = m->Kd;

                if(texture_path[1] == '0')
                {
                    int id = 0;
                    float kh = 0.2, kn = 0.1;
                    float u = st.x;
                    float v = st.y;
                    Vector3f t;
                    Vector3f b;

                    if (std::fabs(N.x) > std::fabs(N.y))
                    {
                        float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
                        b = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
                    }
                    else 
                    {
                        float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
                        b = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
                    }
                    t = crossProduct(b, N);

                    //std::cout << eye_distance << std::endl;
                    if(eye_distance <= 400)
                    {
                        mycolor = tex[0].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[2].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[3].getLinearColor(st.x,st.y).x / 255.f);
                        id = 3;
                    }
                    else if(eye_distance <= 550)
                    {
                        mycolor = tex[0].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[2].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[16].getLinearColor(st.x,st.y).x / 255.f);
                        id = 0;
                    }
                    else if(eye_distance <= 700)
                    {
                        mycolor = tex[4].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[10].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[16].getLinearColor(st.x,st.y).x / 255.f);
                        id = 0;
                    }
                    else if(eye_distance <= 800)
                    {
                        mycolor = tex[5].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[11].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[17].getLinearColor(st.x,st.y).x / 255.f);
                        id = 0;
                    }
                    else if(eye_distance <= 950)
                    {
                        mycolor = tex[6].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[12].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[18].getLinearColor(st.x,st.y).x / 255.f);
                        id = 0;
                    }
                    else if(eye_distance <= 1000)
                    {
                        mycolor = tex[7].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[13].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[19].getLinearColor(st.x,st.y).x / 255.f);
                        id = 0;
                    }
                    else
                    {
                        mycolor = tex[8].getLinearColor(st.x,st.y) / 255.f;
                        roughness = 1.0f - tex[14].getLinearColor(st.x,st.y).x / 255.f;
                        //Kd = (tex[20].getLinearColor(st.x,st.y).x / 255.f);
                        id = 0;
                    }
                    //mycolor = tex[9].getLinearColor(st.x,st.y) / 255.f;
                    if(id == 0)
                    {

                    }
                    else
                    {
                        float w = tex[id].width; 
                        float h = tex[id].height;
                        float dU = kh * kn * (norm(tex[id].getLinearColor(u + 1.f / w, v)) - norm(tex[id].getLinearColor(u,v)));
                        float dV = kh * kn * (norm(tex[id].getLinearColor(u, v + 1.f / h)) - norm(tex[id].getLinearColor(u,v)));
                        Vector3f ln = {-dU, -dV, 1.0f};
                        hitPoint = hitPoint + kn * norm(tex[id].getLinearColor(u,v)) * N;
                        N = -dU * t + - dV * b + N;
                        N = normalize(N);
                    }
                }
                else if(texture_path[1] == '1')
                {
                    mycolor = tex[1].getLinearColor(st.x,st.y) / 255.f;
                }
                else
                {mycolor = Vector3f(0.7,0.7,0.7);}

                Vector3f lightAmt = 0, specularColor = 0;
                Vector3f shadowPointOrig = (dotProduct(ray.direction, N) < 0) ?
                                           hitPoint + N * EPSILON:
                                           hitPoint - N * EPSILON;
                for (uint32_t i = 0; i < get_lights().size(); ++i)
                {
                    auto area_ptr = dynamic_cast<AreaLight*>(this->get_lights()[i].get());
                    if (area_ptr)
                    {
                        // Do nothing for this assignment
                    }
                    else if(true)
                    {
                        //lightDir = normalize(lightDir);
                        Vector3f lightDir = normalize(Vector3f(0.0, 1.0, -0.4));

                        Object *shadowHitObject = nullptr;
                        float tNearShadow = kInfinity;
                        //bool inShadow = bvh->Intersect(Ray(shadowPointOrig, lightDir)).happened;
                        Intersection tempint = bvh->Intersect(Ray(shadowPointOrig, lightDir));
                        float inShadow = 0;
                        if(tempint.happened)
                        {
                            if(tempint.m->getType()==DIFFUSE_AND_GLOSSY)
                            {inShadow = 1.0;}
                            else
                            {inShadow = 0.45;}
                        }

                        float LdotN = std::max(0.f, dotProduct(lightDir, N));
                        hitColor = (1 - inShadow) * 4.0 * m->eval(-lightDir, -ray.direction, N, mycolor, roughness, Kd) * LdotN;
                    }
                    else
                    {
                        Vector3f lightDir = get_lights()[i]->position - hitPoint;
                        // square of the distance between hitPoint and the light
                        float lightDistance2 = dotProduct(lightDir, lightDir);
                        lightDir = normalize(lightDir);
                        float LdotN = std::max(0.f, dotProduct(lightDir, N));
                        Object *shadowHitObject = nullptr;
                        float tNearShadow = kInfinity;
                        //bool inShadow = bvh->Intersect(Ray(shadowPointOrig, lightDir)).happened;
                        Intersection tempint = bvh->Intersect(Ray(shadowPointOrig, lightDir));
                        float inShadow = 0;
                        if(tempint.happened)
                        {
                          if(tempint.m->getType()==DIFFUSE_AND_GLOSSY)
                          {inShadow = 1.0;}
                          else
                          {inShadow = 0.45;}
                        }
                        lightAmt += (1 - inShadow) * get_lights()[i]->intensity * LdotN / lightDistance2 * 50000;
                        Vector3f reflectionDirection = reflect(-lightDir, N);
                        specularColor += powf(std::max(0.f, -dotProduct(reflectionDirection, ray.direction)),
                        m->specularExponent) * get_lights()[i]->intensity;
                        hitColor = lightAmt * (mycolor * m->Kd + specularColor * m->Ks) + Vector3f(0.10, 0.10, 0.10);
                    }
                }

        switch (m->getType()) {
            case REFLECTION_AND_REFRACTION:
            {
                break;
            }
            case REFLECTION:
            {
                hitColor = hitColor * 0.3;
                float kr;
                fresnel(ray.direction, N, m->ior, kr);
                Vector3f reflectionDirection = reflect(ray.direction, N);
                Vector3f reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?
                                             hitPoint - N * EPSILON :
                                             hitPoint + N * EPSILON;
                hitColor += castRay(Ray(reflectionRayOrig, reflectionDirection),depth + 1, tex) * kr * 0.5;
                break;
            }
            default:
            {
                break;
            }
        }
    }

    return hitColor;
}
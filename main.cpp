#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{
    //Scene scene(1280, 960);
    Scene scene(700, 700);

    //REFLECTION_AND_REFRACTION; REFLECTION; DIFFUSE_AND_GLOSSY;

    auto m0 =new Material(MaterialType::DIFFUSE_AND_GLOSSY, Vector3f(0.5, 0.5, 0.5), Vector3f(0, 0, 0));//材质信息
    m0->Kd = 0.7;
    m0->map_Kd = "P0.jpg";   
    m0->map_reflect = "P2.jpg";    
    m0->map_bump = "P3.jpg";    
    m0->Ks = 0.15;
    m0->specularExponent = 50;
    m0->F0 = Vector3f(0.25, 0.25, 0.25);
    m0->roughness = 0.8;
    m0->metalness = 0.5;
    m0->has_map = true;

    Vector3f verts1 = {-1000,-1,-1000};//立方体最小的x,y,z坐标
    Vector3f verts2 = {1000,0,1000};//立方体最大的x,y,z坐标
    Vector2f uv0 = {0,0};//立方体uv贴图范围
    int num0[3] = {10,1,10};//立方体三角形分割

    MeshTriangle ground(verts1,verts2,num0,m0,uv0);//建立一个立方体模型

    scene.Add(&ground);

    auto m1 =new Material(MaterialType::DIFFUSE_AND_GLOSSY, Vector3f(0.5, 0.5, 0.5), Vector3f(0, 0, 0));
    m1->Kd = 1.0;
    m1->map_Kd = "P1.jpg";       
    m1->Ks = 0.15;
    m1->specularExponent = 50;
    m1->F0 = Vector3f(0.25, 0.25, 0.25);
    m1->roughness = 0.8;
    m1->metalness = 0.5;
    m1->has_map = false;


    Vector3f verts3 = {-50,0,-80};
    Vector3f verts4 = {50,0.1,20};
    Vector2f uv1 = {0,0};
    int num1[3] = {1,1,1};

    MeshTriangle park(verts3,verts4,num1,m1,uv1);

    scene.Add(&park);

    scene.Add(std::make_unique<Light>(Vector3f(-20, 180, 100), 0.7));//加入灯光，坐标和强度
    scene.Add(std::make_unique<Light>(Vector3f(20, 100, -100), 0.65));
    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);//进入渲染函数
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}
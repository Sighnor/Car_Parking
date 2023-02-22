//
// Created by LEI XU on 4/27/19.
//

#ifndef RAYTRACING_TEXTURE_H
#define RAYTRACING_TEXTURE_H
#include "global.hpp"
#include "Vector.hpp"
#include <opencv2/opencv.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\objdetect\objdetect.hpp>
#include <opencv2\imgproc\types_c.h>

class Texture{
private:
    cv::Mat image_data;

public:
    Texture(){}
    Texture(const std::string& name)
    {
        image_data = cv::imread(name);
        cv::cvtColor(image_data, image_data, cv::COLOR_RGB2BGR);
        width = image_data.cols;
        height = image_data.rows;
    }

    int width, height;

    Vector3f getColor(float u, float v)
    {
        auto u_img = u * width;
        auto v_img = (1 - v) * height;
        auto color = image_data.at<cv::Vec3b>(v_img, u_img);
        return Vector3f(color[0], color[1], color[2]);
    }

    Vector3f getLinearColor(float u, float v)
    {
        auto u_img = u * (width-1);
        auto v_img = (1 - v) * (height-1);
        //std::cout << u_img << ',' << v_img << '\n';
        float x1 = std::floor(u_img);
        float y1 = std::floor(v_img);
        float x2 = std::ceil(u_img);
        float y2 = std::ceil(v_img);
        if(x1 > width -1)
        {
            x1 = width - 1;
        }
        if (x2 > width - 1)
        {
            x2 = width - 1;
        }
        if (y1 > height - 1)
        {
            y1 = height - 1;
        }
        if (y2 > height - 1)
        {
            y2 = height - 1;
        }
        if (x1 < 0)
        {
            x1 = 0;
        }
        if (x2 < 0)
        {
            x2 = 0;
        }
        if (y1 < 0)
        {
            y1 = 0;
        }
        if (y2 < 0)
        {
            y2 = 0;
        }
        auto color1 = image_data.at<cv::Vec3b>(y1, x1);
        auto color2 = image_data.at<cv::Vec3b>(y1, x2);
        auto color3 = image_data.at<cv::Vec3b>(y2, x1);
        auto color4 = image_data.at<cv::Vec3b>(y2, x2);
        auto color = (x2 - u_img)*(y2 - v_img)*color1 + (u_img - x1)*(y2 - v_img)*color2 + (x2 - u_img)*(v_img - y1)*color3 + (u_img -x1)*(v_img -y1)*color4;
        return Vector3f(color[0], color[1], color[2]);
    }
   
};
#endif //RAYTRACING_TEXTURE_H

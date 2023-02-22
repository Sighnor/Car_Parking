#include <iostream>
#include <opencv2\opencv.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\objdetect\objdetect.hpp>
#include <opencv2\imgproc\types_c.h>
std::string name_color[7] = {"red", "origin", "yellow", "green", "cyan", "blue", "purple"};
cv::Scalar scalar_color[7] = {cv::Scalar(0, 0, 255), cv::Scalar(0, 128, 255), cv::Scalar(0, 255, 255), 
cv::Scalar(0, 255, 0), cv::Scalar(255, 255, 0), cv::Scalar(255, 0, 0), cv::Scalar(255, 0, 128)};//存储七种颜色的RGB
int my_H[8] = {10, 25, 35, 77, 99, 130, 160, 180};//存储七种颜色色调的阈值

void CVT_HSV(int data, void* usrdata)//操作函数，data无意义，usrdata指针为经过运动学检测后的HSV图像
{
    cv::Mat edge_img;//二值图
    cv::imshow("edge", edge_img);
    //cv::Rect my_rect = Find_rect(edge_img);
    //cv::rectangle(BGR_img, my_rect, cv::Scalar(0, 0 ,255), 3);
    //寻找二值图像素最小外接矩形(在合适的运动学检测后，图像中有效像素只有运动的物体，用矩形将其包围，进行运动学计算)
    /*if(my_rect.area() > 10000 && my_rect.area() < 250000)//过滤过小或过大的矩形
    {
        if(temp_width == 0)//若还未赋值
        {
            temp_width = my_rect.width;
            temp_rect = my_rect;
        }
        float p_vec = 0;//速度的比例
        float d_vec = 0;//速度的微分
        float p_angle = 320.0 - (my_rect.tl().x + my_rect.br().x) / 2.0;//角度的比例，图像中心点横坐标减去矩形中心点横坐标，即运动物体与相机的相对横坐标
        i_angle += p_angle;//角度积分
        float temp_p_angle = 320.0 - (temp_rect.tl().x + temp_rect.br().x) / 2.0;
        float d_angle = p_angle - temp_p_angle;//角度的微分，即角速度
        //if(my_rect.area() / float(temp_rect.area()) > 0.9 && my_rect.area() / float(temp_rect.area()) < 1.1)//限制矩形变化速度
        //{
            cv::rectangle(BGR_img, my_rect, cv::Scalar(0, 0 ,255), 3);//在分割后BGR图上绘制矩形
            if(abs(d_angle) < 10)//在角速度不大时，进行速度的计算(一定的角速度会使矩形膨胀，带来很大误差)
            {
                float p_vec = 1.0 - (my_rect.width / float(temp_width));//速度的比例，通过两个矩形宽度的比例，判断运动物体是在靠近还是远离相机
                i_vec += p_vec;//速度的积分
                d_vec = p_vec - temp_vec;//速度的微分
                std::cout << "p_vec:"  << p_vec << ',' << "i_vec:"  << i_vec << ',' << "d_vec:"  << d_vec << ',';
                temp_width = my_rect.width;//存储当前矩形宽度
                temp_vec = p_vec;//存储当前速度
            }
            std::cout << "p_angle:" << p_angle << ',' << "i_angle:" << i_angle << ',' << "d_angle:" << d_angle << '\n';
            cv::imshow("BGR_identify",BGR_img); 
            temp_rect = my_rect;//存储当前矩形
        //}
    }*/
    //cv::imshow("BGR_identify",BGR_img); 
}


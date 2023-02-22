//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include "Texture.hpp"
#include "omp.h"

inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00001;

int my_choice = 1;

std::vector<Vector3f> ray_div;
std::vector<int> my_length;
int min_H = 99;//色调最低阈值
int max_H = 130;//色调最高阈值
int min_S = 43;//饱和度最低阈值
int min_V = 46;//强度最低阈值
float p_vec = 0;//速度的比例
float i_vec = 0;//速度的积分
float d_vec = 0;//速度的微分
float p_omega = 0;//角度的比例
float i_omega = 0;//角度的积分
float d_omega = 0;//角度的微分
float vec_rotate = 0.25;
int flag_find = 0;//记录是否找到标志物
int flag_lost = 0;
int time_count = 10;
int range_boundingbox = 20;
int frame_count = 0;
float center_x = 0;
float center_y = 0;
float angle = 0;
int turn_omega = 0;
int total_height = 0;
int total_width = 0;
float turn_change = 0;
int turn_flag = 0;

int T1 = 3;
int T2 = 9;

cv::Mat filter(cv::Mat before)//滤波
{
    cv::Mat after = before.clone();
    std::vector<int> my_left(before.rows);//存储每行最左侧像素位置
    std::vector<int> my_right(before.rows);//存储每行最右侧像素位置

        //对每行插值处理断点
        for (int i = 0; i < before.rows; i++)
        {
            int color = 0;//颜色，0代表黑色，1代表有颜色
            int change = 0;//颜色微分
            //对于每个断点，寻找其左侧有效像素和右侧有效像素，利用距离进行插值
            int begin_x = -1;
            int end_x = -1;
            for (int j = 0; j < before.cols; j++)
            {
                if(before.at<uchar>(i, j * 3 + 0) != 0 || before.at<uchar>(i, j * 3 + 1) != 0 || before.at<uchar>(i, j * 3 + 2) != 0)
                {
                    change = 1 - color;
                    color = 1;
                }
                else
                {
                    change = - color;
                    color = 0;
                }
                if(change == -1)//颜色由彩色变成黑色
                {
                    begin_x = j - 1;//该点 - 1即左侧有效像素
                    my_right[i] = j;//该点有可能为最右侧像素
                }
                else if(change == 1)//颜色由黑色变成彩色
                {
                    if(begin_x == -1)//如果这是第一个有效像素，则它是最左侧像素s
                    {
                        my_left[i] = j;
                    }
                    else//它为右侧有效像素
                    {
                        end_x = j;
                        for(int k = begin_x + 1; k < end_x; k++)//对中间断点进行插值
                        {
                            after.at<uchar>(i, k * 3 + 0) = (end_x - k) / float(end_x - begin_x) * before.at<uchar>(i, begin_x * 3 + 0) + (k - begin_x) / float(end_x - begin_x) * before.at<uchar>(i, end_x* 3 + 0);
                            after.at<uchar>(i, k * 3 + 1) = (end_x - k) / float(end_x - begin_x) * before.at<uchar>(i, begin_x * 3 + 1) + (k - begin_x) / float(end_x - begin_x) * before.at<uchar>(i, end_x* 3 + 1);
                            after.at<uchar>(i, k * 3 + 2) = (end_x - k) / float(end_x - begin_x) * before.at<uchar>(i, begin_x * 3 + 2) + (k - begin_x) / float(end_x - begin_x) * before.at<uchar>(i, end_x* 3 + 2);
                        }
                    }
                }
            }
            int j = before.cols - 1;//处理最右边像素为终点情况
            if(before.at<uchar>(i, j * 3 + 0) != 0 || before.at<uchar>(i, j * 3 + 1) != 0 || before.at<uchar>(i, j * 3 + 2) != 0)
            {
                my_right[i] = j;
            }
        }
        //对每列进行插值
        cv::Mat my_img = after.clone();
        int j = after.cols / 2;//中间列开始
        int color = 0;
        int change = 0;
        int top_y = -1;
        int below_y = -1;
        //寻找断行进行处理
        for (int i = 0; i < after.rows; i++)
        {
            if(after.at<uchar>(i, j * 3 + 0) != 0 || after.at<uchar>(i, j * 3 + 1) != 0 || after.at<uchar>(i, j * 3 + 2) != 0)
            {
                change = 1 - color;
                color = 1;
            }
            else
            {
                change = - color;
                color = 0;
            }
            if(change == -1)
            {
                top_y = i - 1;
            }
            else if(change == 1 && top_y != -1)
            {
                below_y = i;
                //对每个断行插值
                for(int k = top_y + 1; k < below_y; k++)
                {
                    int long_width = my_right[top_y] - my_left[top_y];//上方有效行长度
                    int short_width = my_right[below_y] - my_left[below_y];//下方有效行长度s
                    float top_height = k - top_y;//到上方的距离
                    float below_height = below_y - k;//到下方的距离s
                    float y1 = below_height / ((top_height + below_height));
                    float y2 = 1 - y1;
                    float left_x = my_left[top_y] * y1 + my_left[below_y] * y2;//插值得到最左侧坐标
                    float right_x = my_right[top_y] * y1 + my_right[below_y] * y2;//插值得到最右侧坐标s
                    float x_width = right_x - left_x;//当前行长度
                    for(int n = int(left_x + 1); n < int(right_x); n++)
                    {
                        float ratio_x = (n - left_x) / x_width;//该点在当前行比例
                        int top_x = my_left[top_y] + long_width * ratio_x;//用比例算出上方行对应点
                        int below_x = my_left[below_y] + short_width * ratio_x;//用比例算出下方行对应点
                        //插值计算像素
                        my_img.at<uchar>(k, n * 3 + 0) = after.at<uchar>(top_y, top_x * 3 + 0) * y1 + after.at<uchar>(below_y, below_x * 3 + 0) * y2;
                        my_img.at<uchar>(k, n * 3 + 1) = after.at<uchar>(top_y, top_x * 3 + 1) * y1 + after.at<uchar>(below_y, below_x * 3 + 1) * y2;
                        my_img.at<uchar>(k, n * 3 + 2) = after.at<uchar>(top_y, top_x * 3 + 2) * y1 + after.at<uchar>(below_y, below_x * 3 + 2) * y2;
                    }
                }
            }
        }
    return my_img;
}

bool Inside_boundingbox(int x, int y, cv::Rect bounding_box)
{
    return (x > bounding_box.tl().x - range_boundingbox && x < bounding_box.br().x + range_boundingbox
    && y > bounding_box.tl().y - range_boundingbox && y < bounding_box.br().y + range_boundingbox) 
    ? true : false;
}

cv::Rect Find_rect(cv::Mat img)//寻找最小外接矩形
{
    int num_boundingbox = 0;
    std::vector<int> num_pixel;
    std::vector<cv::Rect> my_boundingbox;
    int max_pixel = 80;
    int my_i = -1;
    for(int i = 0; i < img.rows; i++)  
    {  
        uchar* temp = img.ptr<uchar>(i);
        for(int j = 0; j < img.cols; j++)   
        {  
            if(temp[j] > 0)
            {
                int flag = 0;
                for(int k = 0;k < num_boundingbox && flag == 0; k++)
                {
                    if(Inside_boundingbox(j, i, my_boundingbox[k]))
                    {
                        flag = 1;
                        int left_x = my_boundingbox[k].tl().x;
                        int right_x = my_boundingbox[k].br().x;
                        int top_y = my_boundingbox[k].tl().y;
                        int below_y = my_boundingbox[k].br().y;
                        if(right_x < j)
                        {
                            right_x = j;
                        }
                        else if(left_x > j)
                        {
                            left_x = j;
                        }
                        if(below_y < i)
                        {
                            below_y = i;
                        }
                        else if(top_y > i)
                        {
                            top_y = i;
                        }
                        my_boundingbox[k] = cv::Rect(cv::Point(left_x, top_y), cv::Point(right_x, below_y));
                        num_pixel[k]++;
                    }
                }
                if(flag == 0)
                {
                    my_boundingbox.push_back(cv::Rect(cv::Point(j, i), cv::Point(j, i)));
                    num_pixel.push_back(1);
                    num_boundingbox++;
                }
            }
        }  
    }
    for(int i  = 0; i < num_boundingbox; i++)
    {
        float ratio_height_width = my_boundingbox[i].height / float(my_boundingbox[i].width);
        if(num_pixel[i] > max_pixel)// && ratio_height_width < 2.0f && ratio_height_width > 0.5f && my_boundingbox[i].height < 75 && my_boundingbox[i].width < 75)
        {
            my_i = i;
            max_pixel = num_pixel[i];
            //cv::rectangle(img, my_boundingbox[i], cv::Scalar(0, 0, 255), 1);
        }
    } 
    if (my_i != -1)
    {
        return my_boundingbox[my_i];
    }
    else
    {
        return cv::Rect(cv::Point(0, 0),cv::Point(720, 720));
    }
}

cv::RotatedRect Find_minrect(cv::Rect bounding_box, cv::Mat img)
{
    std::vector<cv::Point> contours;

        for(int i = bounding_box.tl().y; i < bounding_box.br().y; i++)
        {
            uchar* temp = img.ptr<uchar>(i);
            for(int j = bounding_box.tl().x; j < bounding_box.br().x; j++)
            { 
                if(temp[j] > 0)
                {
                    contours.push_back(cv::Point(j, i));
                }
            }
        }

    return cv::minAreaRect(cv::Mat(contours));
}

void pid_control(cv::Mat BGR_img)
{
    cv::Mat temp_img;
    cvtColor(BGR_img, temp_img, CV_BGR2GRAY);
    threshold(temp_img, temp_img, 20, 255, CV_THRESH_BINARY); //二值化
    cv::Rect bounding_box = Find_rect(temp_img);
    //cv::rectangle(BGR_img, bounding_box, cv::Scalar(0, 255, 0), 1);

    if(bounding_box.area() < 10000 && flag_lost == 0)
    {
        cv::RotatedRect my_bounding_box = Find_minrect(bounding_box, temp_img); //定义最小外接矩形
        cv::Point2f rect[4];

        std::cout << "angle = " << my_bounding_box.angle << '\n';

        flag_find = 1;

        if(frame_count < 10)
        {
            my_bounding_box.points(rect);
            //for(int j=0; j<4; j++)
            //{
                //line(BGR_img, rect[j], rect[(j+1)%4], cv::Scalar(0, 0, 255), 2, 8);  //绘制最小外接矩形每条边
            //}
            total_height += my_bounding_box.size.height;
            total_width += my_bounding_box.size.width;
            frame_count++;
        }
        else if(frame_count == 10)
        {
            total_height = total_height / 10;
            total_width = total_width / 10;
            frame_count++;
        }
        else
        {
            if((my_bounding_box.size.height > 0.50 * total_height && my_bounding_box.size.height < 1.50 * total_height
            && my_bounding_box.size.width > 0.50 * total_width && my_bounding_box.size.width < 1.50 * total_width) ||
	        (my_bounding_box.size.width > 0.50 * total_height && my_bounding_box.size.width < 1.50 * total_height
            && my_bounding_box.size.height > 0.50 * total_width && my_bounding_box.size.height < 1.50 * total_width))
            {
                my_bounding_box.points(rect);
                for(int j=0; j<4; j++)
                {
                    line(BGR_img, rect[j], rect[(j+1)%4], cv::Scalar(255, 255, 255), 2, 8);  //绘制最小外接矩形每条边
                }
                center_x = my_bounding_box.center.x;
                center_y = my_bounding_box.center.y;
                angle = my_bounding_box.angle;
            }

            float distance_tl = abs((BGR_img.rows - center_y) * cos(angle * M_PI / 180.0) - (BGR_img.cols / 2 - center_x) * sin(angle * M_PI / 180.0));
            float distance_tr = abs((BGR_img.rows - center_y) * sin(angle * M_PI / 180.0) + (BGR_img.cols / 2 - center_x) * cos(angle * M_PI / 180.0));
            int temp_y;

            if(distance_tl < distance_tr)
            {
                temp_y = center_y + (BGR_img.cols / 2 - center_x) * tan(angle * M_PI / 180.0);
            }
            else
            {
                temp_y = center_y + (center_x - BGR_img.cols / 2) / tan(angle * M_PI / 180.0);
            }

            cv::line(BGR_img, cv::Point(BGR_img.cols / 2, BGR_img.rows), cv::Point(BGR_img.cols / 2, temp_y), cv::Scalar(0, 0, 255), 1);
            cv::line(BGR_img, cv::Point(center_x, center_y), cv::Point(BGR_img.cols / 2, temp_y), cv::Scalar(255, 0, 0), 1);

            float distance_ot = (BGR_img.rows - temp_y) * (BGR_img.rows - temp_y);
            float distance_tc = ((center_x - BGR_img.cols / 2) * (center_x - BGR_img.cols / 2) + (center_y - temp_y) * (center_y - temp_y));
            int center_flag = center_x > (BGR_img.cols / 2) ? -1 : 1;

            //std::cout << '\n' << center_flag << '\n';

            float temp_turn_change;

            if(turn_flag == 1)
            {
                turn_omega = 0;
            }
            else if(turn_omega != 0)
            {
                temp_turn_change = distance_ot - distance_tc;
            }
            else if(temp_y >= BGR_img.rows || temp_y <= center_y)
            {
                turn_omega = center_flag;
                temp_turn_change = 0;
            }
            else if(distance_ot < distance_tc)
            {
                turn_omega = center_flag;
                temp_turn_change = distance_ot - distance_tc;
            }
            else if(distance_ot >= distance_tc)
            {
                turn_omega = -center_flag;
                temp_turn_change = distance_ot - distance_tc;
            }

            if((angle > 88 || angle < 2) && center_x > BGR_img.cols / 2 - 5 && center_x < BGR_img.cols / 2 + 5)
            {
                std::cout << "move!" << '\n';
                turn_flag = 1;
                turn_omega = 0;
            }
            else if(center_x > BGR_img.cols / 2 - 50 + my_length[center_y] || center_x < BGR_img.cols / 2 + 50 - my_length[center_y])
            {
                std::cout << "must turn!" << '\n';
                //turn_flag = 0;  
                //turn_change = 0;
                turn_omega = center_flag; 
            }
            else if(turn_change * temp_turn_change < 0)
            {
                std::cout << "change!" << '\n';
                turn_flag = 1;
                turn_change = 0;
                turn_omega = 0;
            }
            else
            {
                turn_change = temp_turn_change;
            }

            if(turn_omega != 0)
            {
                p_vec = 0;
                p_omega = 0.25 * turn_omega * vec_rotate;
                std::cout << "turn around!" << '\n';
            }
            else if((angle > 89 || angle < 1) && center_x > BGR_img.cols / 2 - 10 && center_x < BGR_img.cols / 2 + 10)
            {
                p_vec = 20;
                p_omega = 0;
                std::cout << "straight!" << '\n';
            }
            else
            {
                for(float t = 0;t < 1;t += 0.02)
                {
                    float my_x = t * t * BGR_img.cols / 2 + 2 * t * (1 - t) * BGR_img.cols / 2 + (1 - t) * (1 - t) * center_x;
                    float my_y = t * t * BGR_img.rows + 2 * t * (1 - t) * temp_y + (1 - t) * (1 - t) * center_y;
                    cv::circle(BGR_img, cv::Point(my_x, my_y), 1, cv::Scalar(255, 255, 255));
                }
                float t = 0.9;
                float my_x = t * t * BGR_img.cols / 2 + 2 * t * (1 - t) * BGR_img.cols / 2 + (1 - t) * (1 - t) * center_x;
                float my_y = t * t * BGR_img.rows + 2 * t * (1 - t) * temp_y + (1 - t) * (1 - t) * center_y;
                if(temp_y > BGR_img.rows - 30)
                {
                    p_vec = 0;
                    if(my_x < BGR_img.cols/2)
                    {
                        p_omega = 0.25 * vec_rotate;
                    }
                    else
                    {
                        p_omega = -0.25 * vec_rotate;
                    }
                }
                else
                {
                    p_vec = 3 * (BGR_img.rows - my_y);
                    p_omega =  atan((BGR_img.cols/2 - my_x) / float(BGR_img.rows - my_y));
                }
                std::cout << "Bezier!" << '\n';
            }
	    }
    }
    else
    {
        if(frame_count > 10)
        {
            flag_lost = 1;
            if(time_count > 0)
            {
                std::cout << "go!" << '\n';
                time_count--;
            }
            else
            {
                std::cout << "stop!" << '\n';
                p_vec = 0;
                p_omega = 0;
            }
        }
        else
        {
            p_omega = -0.50 * vec_rotate;
        }
    }
    cv::imshow("pid", BGR_img);
}

cv::Mat color_cut(cv::Mat img)
{
        cv::Mat BGR_img = img.clone();
        cv::Mat HSV_img;
        cv::cvtColor(BGR_img, HSV_img, cv::COLOR_BGR2HSV);
        for(int y = 0; y < HSV_img.rows; y++)//颜色阈值滤波
        {
            for(int x = 0; x < HSV_img.cols; x++)
            {
                        int temp_H = HSV_img.at<uchar>(y, x * 3 + 0);
                        int temp_S = HSV_img.at<uchar>(y, x * 3 + 1);
                        int temp_V = HSV_img.at<uchar>(y, x * 3 + 2);
                        if(temp_V < min_V || temp_S < min_S)
                        {
                            BGR_img.at<uchar>(y, x * 3 + 0) = 0; 
                            BGR_img.at<uchar>(y, x * 3 + 1) = 0; 
                            BGR_img.at<uchar>(y, x * 3 + 2) = 0; 
                        }
                        else
                        {
                            if(temp_H > 99 && temp_H <= 130 || (temp_H > 35 && temp_H <= 77))
                            {
                    
                            }
                            else
                            {
                                BGR_img.at<uchar>(y, x * 3 + 0) = 0; 
                                BGR_img.at<uchar>(y, x * 3 + 1) = 0; 
                                BGR_img.at<uchar>(y, x * 3 + 2) = 0; 
                            } 
                        }
            }
        } 
        return BGR_img;
}

void dip(float eye_pos_y, float ratio, int height, int width, cv::Mat img)
{    
        cv::Mat BGR_img(width, width, CV_8UC3, cv::Scalar(0,0,0));//初始化为黑色 
        int flag_y = 0;
        for(int y = img.rows - 1; y > img.rows / 2 && flag_y == 0; y--)//在下半边寻找
        {
            uchar* temp = img.ptr<uchar>(y);
            float distance = - eye_pos_y * ray_div[y * img.cols + img.cols / 2].y;//距离
            int temp_y = (width + ratio * distance * ray_div[y * img.cols + img.cols / 2].z);//重构y坐标
            if(temp_y < 0)
            {
                flag_y = 1;
            }
            else
            {
                int flag_x = 0;
                for(int x = img.cols / 2 - 1; x >= 0 && flag_x == 0; x--)
                {
                        int temp_x = (ratio * distance * ray_div[y * img.cols + x].x + width / 2.0f);//重构x坐标
                        if(temp_x < 0)
                        {
                            flag_x = 1;
                        }
                        else
                        {
                            BGR_img.at<uchar>(temp_y, temp_x * 3 + 0) = temp[3 * x];
                            BGR_img.at<uchar>(temp_y, temp_x * 3 + 1) = temp[3 * x + 1]; 
                            BGR_img.at<uchar>(temp_y, temp_x * 3 + 2) = temp[3 * x + 2]; 
                        }
                }
                flag_x = 0;
                for(int x = img.cols / 2; x < img.cols && flag_x == 0; x++)
                {
                        int temp_x = (ratio * distance * ray_div[y * img.cols + x].x + width / 2.0f);//重构x坐标
                        if(temp_x >= width)
                        {
                            flag_x = 1;
                        }
                        else
                        {
                            BGR_img.at<uchar>(temp_y, temp_x * 3 + 0) = temp[3 * x];
                            BGR_img.at<uchar>(temp_y, temp_x * 3 + 1) = temp[3 * x + 1]; 
                            BGR_img.at<uchar>(temp_y, temp_x * 3 + 2) = temp[3 * x + 2]; 
                        }
                }
            }
        }
        //cv::imshow("beforefilter", BGR_img);
        BGR_img = filter(BGR_img);//特殊均值滤波，处理图像断线
        cv::imshow("afterfilter", BGR_img);   
        BGR_img = color_cut(BGR_img);
        //cv::imshow("final", BGR_img);
        
        pid_control(BGR_img);//pid控制
        //cv::Mat edge_img;
        //cv::Canny(BGR_img, edge_img, T1, T2);
        //cv::imshow("canny", edge_img);
}

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    float scale = tan(deg2rad(scene.fov * 0.5));//FOV对应的tan
    float imageAspectRatio = scene.width / (float)scene.height;//x,y比例
    Texture tex[23];//贴图信息
    tex[0] = Texture("./textures/P0.jpg" );
    tex[1] = Texture("./textures/P1.jpg" );
    tex[2] = Texture("./textures/P2.jpg" );
    tex[3] = Texture("./textures/P3.jpg" );
    tex[4] = Texture("./textures/mipmap1.jpg" );
    tex[5] = Texture("./textures/mipmap2.jpg" );
    tex[6] = Texture("./textures/mipmap3.jpg" );
    tex[7] = Texture("./textures/mipmap4.jpg" );
    tex[8] = Texture("./textures/mipmap5.jpg" );
    tex[9] = Texture("./textures/mipmap6.jpg" );
    tex[10] = Texture("./textures/reflectmipmap1.jpg" );
    tex[11] = Texture("./textures/reflectmipmap2.jpg" );
    tex[12] = Texture("./textures/reflectmipmap3.jpg" );
    tex[13] = Texture("./textures/reflectmipmap4.jpg" );
    tex[14] = Texture("./textures/reflectmipmap5.jpg" );
    tex[15] = Texture("./textures/reflectmipmap6.jpg" );
    tex[16] = Texture("./textures/bumpmipmap1.jpg" );
    tex[17] = Texture("./textures/bumpmipmap2.jpg" );
    tex[18] = Texture("./textures/bumpmipmap3.jpg" );
    tex[19] = Texture("./textures/bumpmipmap4.jpg" );
    tex[20] = Texture("./textures/bumpmipmap5.jpg" );
    tex[21] = Texture("./textures/bumpmipmap6.jpg" );
    tex[22] = Texture("./textures/mipmap0.jpg" );
   
    int key = 0;

    std::vector<Vector3f> framebuffer(scene.width * scene.height);//存储像素

    int m = 0;
    Vector3f axis(0,0,-1);//初始视线方向
    Vector3f view_ray(0,0,-1);//目标视线方向
    Vector3f eye_pos(0, 100, 75);//初始相机位置
    float rotation_angle = 0.0;

    FILE* myfp = fopen("eyepos&angle.txt", "r");
    (void)fscanf(myfp, "%f\n%f\n%f\n%f\n%f\n%f\n", &eye_pos.x, &eye_pos.y, &eye_pos.z, &view_ray.x, &view_ray.y, &view_ray.z);
    fclose(myfp);

    /*std::cout << "eye_pos.x = " << eye_pos.x << " , Input new eye_pos.x : ";
    std::cin >> eye_pos.x;
    std::cout << "eye_pos.y = " << eye_pos.y << " , Input new eye_pos.y : ";
    std::cin >> eye_pos.y;
    std::cout << "eye_pos.z = " << eye_pos.z << " , Input new eye_pos.z : ";
    std::cin >> eye_pos.z;
    std::cout << "view_ray.x = " << view_ray.x << " , Input new view_ray.x : ";
    std::cin >> view_ray.x;
    std::cout << "view_ray.y = " << view_ray.y << " , Input new view_ray.y : ";
    std::cin >> view_ray.y;
    std::cout << "view_ray.z = " << view_ray.z << " , Input new view_ray.z : ";
    std::cin >> view_ray.z;*/

    view_ray = normalize(view_ray);//视线向量化为单位向量s
    rotation_angle = std::acos(dotProduct(axis,view_ray));//点乘得到旋转角
    axis = normalize(crossProduct(axis,view_ray));//叉乘得到旋转轴

    myfp = fopen("eyepos&angle.txt", "w");
    (void)fprintf(myfp, "%f\n%f\n%f\n%f\n%f\n%f\n", eye_pos.x, eye_pos.y, eye_pos.z, view_ray.x, view_ray.y, view_ray.z);
    fclose(myfp);

    std::cout << "use pid or key?(0 or 1)" << '\n';
    std::cin >> my_choice;

    //罗德里格斯旋转矩阵
    float c = std::cos(rotation_angle);
    float k = 1-c;
    float s = std::sin(rotation_angle);

    float kxx = k * axis.x * axis.x + c;
    float kxy = k * axis.x * axis.y - axis.z * s;
    float kxz = k * axis.x * axis.z + axis.y * s;
    float kyx = k * axis.x * axis.y + axis.z * s;
    float kyy = k * axis.y * axis.y + c;
    float kyz = k * axis.y * axis.z - axis.x * s;
    float kzx = k * axis.x * axis.z - axis.y * s;
    float kzy = k * axis.y * axis.z + axis.x * s;
    float kzz = k * axis.z * axis.z + c;

    std::vector<Vector3f> dir(scene.height * scene.width);

    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                      imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;
            // TODO: Find the x and y positions of the current pixel to get the
            // direction
            //  vector that passes through it.
            // Also, don't forget to multiply both of them with the variable
            // *scale*, and x (horizontal) variable with the *imageAspectRatio*
            //Vector3f dir = Vector3f(x, y, -1.0f); // Don't forget to normalize this direction!
            dir[m] = normalize(Vector3f(x, y, -1.0f));//得到各个像素的对应的视线向量
            //ray_div.push_back(Vector3f(dir[m].x, 1./dir[m].y, dir[m].z));
            ray_div.push_back(Vector3f(x, 1./y, -1.0f));//ray_div数组存储视线向量，y存储向量y的倒数，没有规范化方便后续计算
            m++;
            // Don't forget to normalize this direction!
        }
    }

    /*float myscale = tan(deg2rad(90 * 0.5));//FOV对应的tan
    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                      imageAspectRatio * myscale;
            float y = (1 - 2 * (j + 0.5) / (float)scene.height) * myscale;
            ray_div.push_back(Vector3f(x, 1./y, -1.0f));//ray_div数组存储视线向量，y存储向量y的倒数，没有规范化方便后续计算
        }
    }*/

    cv::Mat test_img(scene.height, scene.width, CV_8UC3, cv::Scalar(0,0,0));
    cv::imshow("image", test_img);
    //cv::createTrackbar("min_S", "image", &min_S, 255);//设定滑动模块，可动态调整相关参数
    //cv::createTrackbar("min_V", "image", &min_V, 255);
    cv::createTrackbar("T1", "image", &T1, 255);
    cv::createTrackbar("T2", "image", &T2, 255);

    while (key != 27) {//按ESC退出循环

        std::fill(framebuffer.begin(), framebuffer.end(), Vector3f{0, 0, 0});

        #pragma omp parallel for 
        for (uint32_t j = 0; j < scene.height; ++j) {
            for (uint32_t i = 0; i < scene.width; ++i) {
                Vector3f newdir;
                int temp_m = i + j * scene.width;
                newdir.x = kxx * dir[temp_m].x + kxy * dir[temp_m].y + kxz * dir[temp_m].z;//视线旋转
                newdir.y = kyx * dir[temp_m].x + kyy * dir[temp_m].y + kyz * dir[temp_m].z;
                newdir.z = kzx * dir[temp_m].x + kzy * dir[temp_m].y + kzz * dir[temp_m].z;
                dir[temp_m] = newdir;//存储新视线
                framebuffer[temp_m] = scene.castRay(Ray(eye_pos,dir[temp_m]), 0, tex) * 255.f;//对每个像素光线追踪，返回一个rgb像素值
                // Don't forget to normalize this direction!

            }
        }

        cv::Mat image(scene.height, scene.width, CV_32FC3, framebuffer.data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imshow("image", image);

        //cv::cvtColor(image, image, cv::COLOR_BGR2HSV);
        
        //二维重构
        int temp_height = eye_pos.y / scale;//暂时没用
        int temp_width = scene.width;//重构图像宽度
        float temp_ratio = (temp_width * 0.25f) / (2 * eye_pos.y * imageAspectRatio);//比例，最近点在图中占0.25倍长度


        float temp_scale = scale * imageAspectRatio;
        for(int i = 0;i < temp_width;i++)
        {
            int length = (temp_width - i) * temp_scale;
            if(length > temp_width / 2)
            {
                my_length.push_back(temp_width / 2);
            }
            else
            {
                my_length.push_back(length);
            }
        }

        dip(eye_pos.y, temp_ratio, temp_height, temp_width, image);//重构函数

        key = cv::waitKey(50);

        //pid控制
        float kp_vec = 0.15;
        float ki_vec = 0;
        float kd_vec = 0;
        float kp_omega = 0.25;
        float ki_omega = 0;
        float kd_omega = 0;

        float my_vec = kp_vec * p_vec + ki_vec * i_vec + kd_vec * d_vec;
        float my_omega = kp_omega * p_omega + ki_omega * i_omega + kd_omega * d_omega;

        if(my_vec > 40)//限制线速度
        {
            my_vec = 40;
        }
        if(my_omega > 0.15)//限制角速度
        {
            my_omega = 0.15;
        }

        std::cout << "vec = " << my_vec << ',' << "omega = " << my_omega << '\n' << '\n';

        switch (my_choice)
        {
            case 0:
                //视线绕y轴旋转,下一次循环中使用
                kxx = std::cos(my_omega);
                kxy = 0;
                kxz = std::sin(my_omega);
                kyx = 0;
                kyy = 1;
                kyz = 0;
                kzx = -std::sin(my_omega);
                kzy = 0;
                kzz = std::cos(my_omega);

                //相机中心视线绕y轴旋转
                view_ray.x = kxx * view_ray.x + kxy * view_ray.y + kxz * view_ray.z;
                view_ray.y = kyx * view_ray.x + kyy * view_ray.y + kyz * view_ray.z;
                view_ray.z = kzx * view_ray.x + kzy * view_ray.y + kzz * view_ray.z;

                //相机位置移动
                eye_pos.x += my_vec * view_ray.x; 
                eye_pos.y += my_vec * view_ray.y; 
                eye_pos.z += my_vec * view_ray.z;
            break;
            default:
                kxx = 1;
                kxy = 0;
                kxz = 0;
                kyx = 0;
                kyy = 1;
                kyz = 0;
                kzx = 0;
                kzy = 0;
                kzz = 1;
                if (key == 'a') {
                    kxx = std::cos(M_PI/12);
                    //kxy = 0;
                    kxz = std::sin(M_PI/12);
                    //kyx = 0;
                    //kyy = 1;
                    //kyz = 0;
                    kzx = -std::sin(M_PI/12);
                    //kzy = 0;
                    kzz = std::cos(M_PI/12);

                    view_ray.x = kxx * view_ray.x + kxy * view_ray.y + kxz * view_ray.z;
                    view_ray.y = kyx * view_ray.x + kyy * view_ray.y + kyz * view_ray.z;
                    view_ray.z = kzx * view_ray.x + kzy * view_ray.y + kzz * view_ray.z;
                }
                else if (key == 'd') {
                    kxx = std::cos(M_PI/12);
                    //kxy = 0;
                    kxz = -std::sin(M_PI/12);
                    //kyx = 0;
                    //kyy = 1;
                    //kyz = 0;
                    kzx = std::sin(M_PI/12);
                    //kzy = 0;
                    kzz = std::cos(M_PI/12);

                    view_ray.x = kxx * view_ray.x + kxy * view_ray.y + kxz * view_ray.z;
                    view_ray.y = kyx * view_ray.x + kyy * view_ray.y + kyz * view_ray.z;
                    view_ray.z = kzx * view_ray.x + kzy * view_ray.y + kzz * view_ray.z;
                }
                else if (key == 'w') {
                    eye_pos.x += 25 * view_ray.x; 
                    eye_pos.y += 25 * view_ray.y; 
                    eye_pos.z += 25 * view_ray.z; 
                }
                else if (key == 's') {
                    eye_pos.x -= 25 * view_ray.x; 
                    eye_pos.y -= 25 * view_ray.y; 
                    eye_pos.z -= 25 * view_ray.z; 
                }
            break;
        }

    }
    // save framebuffer to file

    /*FILE* fp = fopen(file_name, "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * clamp(0, 1, framebuffer[i].x));
        color[1] = (unsigned char)(255 * clamp(0, 1, framebuffer[i].y));
        color[2] = (unsigned char)(255 * clamp(0, 1, framebuffer[i].z));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);   
    std::cout << "continue? (y or n) : ";
    std::cin >> choice;
    } */
}

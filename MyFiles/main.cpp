// main.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "MyTiffFile.h"


int main()
{
    std::string tifpath = "D:\\CODE\\DEM_ENHANCEMENT\\Level16\\gagongdem.tif";
    std::string muti_tifpath = "D:\\CODE\\SHP2RASTER\\updatemask\\outfiles\\multi_band_output_mask.tif";
    cv::Mat fit_mat_img;
    MyTiffFile mytif(tifpath);
    // test1: 将该tif转换为cv::Mat
    mytif.toMat(fit_mat_img);
    // test2: 对该tif进行重采样
    std::string dstImg = "D:\\CODE\\MyTiff\\MyFiles\\resample_3m.tif";
    mytif.resample(dstImg, 3, InterpolationType::BILINEAR);
    
    // test3: 对该tif进行重投影
    std::string reProjImg = "D:\\CODE\\MyTiff\\MyFiles\\reProjImg.tif";
    mytif.reProjection(reProjImg.c_str());
    int m = 0;
    return 0;
}


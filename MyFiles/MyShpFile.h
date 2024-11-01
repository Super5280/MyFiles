#pragma once
#include <iostream>
#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>
#include <gdal/gdal_alg.h>
#include "D:\\CODE\\gdal232\\include\\ogrsf_frmts.h"

class MyShpFile
{
public:
	MyShpFile(std::string shpFilePath);
	~MyShpFile();
	void init();
	struct ShpInfo
	{
		char* projection_infos = nullptr; // 投影信息
		OGREnvelope envelope;
	};

	void shp2raster(std::string rasterImg);

public:
	ShpInfo shpinfos;
private:
	std::string shpFilePath;
	GDALDataset* poDS;
	OGRLayer* poLayer;  // 假设该shp文件只有一种类型的要素，面要素
};


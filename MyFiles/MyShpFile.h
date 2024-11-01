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
		char* projection_infos = nullptr; // ͶӰ��Ϣ
		OGREnvelope envelope;
	};

	void shp2raster(std::string rasterImg);

public:
	ShpInfo shpinfos;
private:
	std::string shpFilePath;
	GDALDataset* poDS;
	OGRLayer* poLayer;  // �����shp�ļ�ֻ��һ�����͵�Ҫ�أ���Ҫ��
};


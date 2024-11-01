#include "MyShpFile.h"

MyShpFile::MyShpFile(std::string shpFilePath)
{
	this->shpFilePath = shpFilePath;
	this->init();
}

MyShpFile::~MyShpFile()
{
	if (this->poDS)
		GDALClose(this->poDS);
}

void MyShpFile::init()
{
	// 初始化GDAL
	GDALAllRegister();
	OGRRegisterAll();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");   // 支持中文路径
	CPLSetConfigOption("SHAPE_ENCODING", "");

	this->poDS = (GDALDataset*)GDALOpenEx(this->shpFilePath.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (this->poDS == NULL) {
		std::cout << "Failed to open shapefile.\n";
		return;
	}

	// 获取矢量图层
	this->poLayer = poDS->GetLayer(0);
	if (!poLayer) {
		std::cout << "Failed to get layer.\n";
		return;
	}

	// 获取矢量范围
	this->poLayer->GetExtent(&this->shpinfos.envelope);

	// 检查空间参考并转换为有效的 WKT
	if (poLayer->GetSpatialRef() != nullptr) {
		poLayer->GetSpatialRef()->exportToWkt(&this->shpinfos.projection_infos);
	}
	else {
		std::cerr << "空间参考无效，使用默认投影" << std::endl;
		return;
	}

}

void MyShpFile::shp2raster(std::string rasterImg)
{

}

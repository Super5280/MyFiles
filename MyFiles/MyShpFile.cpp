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
	// ��ʼ��GDAL
	GDALAllRegister();
	OGRRegisterAll();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");   // ֧������·��
	CPLSetConfigOption("SHAPE_ENCODING", "");

	this->poDS = (GDALDataset*)GDALOpenEx(this->shpFilePath.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (this->poDS == NULL) {
		std::cout << "Failed to open shapefile.\n";
		return;
	}

	// ��ȡʸ��ͼ��
	this->poLayer = poDS->GetLayer(0);
	if (!poLayer) {
		std::cout << "Failed to get layer.\n";
		return;
	}

	// ��ȡʸ����Χ
	this->poLayer->GetExtent(&this->shpinfos.envelope);

	// ���ռ�ο���ת��Ϊ��Ч�� WKT
	if (poLayer->GetSpatialRef() != nullptr) {
		poLayer->GetSpatialRef()->exportToWkt(&this->shpinfos.projection_infos);
	}
	else {
		std::cerr << "�ռ�ο���Ч��ʹ��Ĭ��ͶӰ" << std::endl;
		return;
	}

}

void MyShpFile::shp2raster(std::string rasterImg)
{

}

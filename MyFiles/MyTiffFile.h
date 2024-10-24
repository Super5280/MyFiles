#pragma once
#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>
#include <gdal/gdalwarper.h>
#include <opencv2/opencv.hpp>
#include <iostream>

// ����뾶����λ����
const double EARTH_RADIUS = 6378137.0;


enum class InterpolationType
{
	// ��ֵ�㷨������
	NEAREST,       // ����ڲ�ֵ
	BILINEAR,      // ˫���Բ�ֵ
	CUBIC,         // ���β�ֵ
	CUBICSPLINE,   // ����������ֵ
};



class MyTiffFile
{
	// ר�����ڹ���tiff�ļ���ʽ�������˶���tiff�ļ���һЩ��������
public:
	MyTiffFile(std::string tiffPath);
	~MyTiffFile();
	struct TiffInfo
	{
		int nSrcXSize;   // width
		int nSrcYSize;   // height
		int nBands;      // band
		GDALDataType gdalType;
		double adfGeoTransform[6]; // ����ת������
		const char* projection_info;  // ͶӰ��Ϣ��WKT��ʽ��
		double dfNoDataValue;         // ��Чֵ
		int bGotNoDataValue = FALSE;

	};
	void resample(std::string dstPath, double targetResolution, InterpolationType interType = InterpolationType::BILINEAR);
	void toMat(cv::Mat& matImg);
	int reProjection(const char* outputFile);

	// �����......
	void splitRegion(std::string dstPath, std::string shpFile);
	GDALDataset* reprojectShpToRaster(OGRLayer* poLayer);  // shp--->raster
	
public:
	TiffInfo tifinfos;   // �洢����tif��һЩ������Ϣ

private:
	void init();
	GDALResampleAlg InterpType2GDALResampleAlg(InterpolationType type);
	void ReprojectShapefile(const char* inputShapefile, const char* outputShapefile, const char* targetWKT);
	// �������ϵ�1�Ⱦ��루�ף�
	double degreesToMetersAtEquator(double degrees) {
		return degrees * 2 * M_PI * EARTH_RADIUS / 360.0;
	}
	// ����γ�ȼ���ʵ�ʵ�1�ȵľ��Ⱦ��루�ף�
	double degreesToMetersAtLatitude(double degrees, double latitude) {
		return degrees * 2 * M_PI * EARTH_RADIUS * cos(latitude * M_PI / 180.0) / 360.0;
	}
	// ����Ŀ��ֱ�������ķŴ���
	double calculateScaleFactor(GDALDataset* poDataset, double targetResolutionMeters);
	// ��Ӱ���л�ȡ����γ��
	double getCenterLatitude(GDALDataset* poDataset);
	std::vector<std::string> getUTMAndWGS84WKT(double longitude, double latitude);
	double* getImageCenterCoordinates();
	OGRLayer* transformLayerCoordinates(OGRLayer* poLayer, OGRCoordinateTransformation* poCT);
	void imageReprojection(const char* pszDstFile, const char* pszDstWKT, InterpolationType eResampleMethod,
		double dResX = 0.0, double dResY = 0.0, const char* pszFormat = "GTiff");

private:
	std::string tiffImgPath;
	GDALDataset* poDataset;
};


#pragma once
#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>
#include <gdal/gdalwarper.h>
#include <opencv2/opencv.hpp>
#include <iostream>

// 地球半径，单位：米
const double EARTH_RADIUS = 6378137.0;


enum class InterpolationType
{
	// 插值算法的类型
	NEAREST,       // 最近邻插值
	BILINEAR,      // 双线性插值
	CUBIC,         // 三次插值
	CUBICSPLINE,   // 三次样条插值
};



class MyTiffFile
{
	// 专门用于管理tiff文件格式，包含了对于tiff文件的一些常见操作
public:
	MyTiffFile(std::string tiffPath);
	~MyTiffFile();
	struct TiffInfo
	{
		int nSrcXSize;   // width
		int nSrcYSize;   // height
		int nBands;      // band
		GDALDataType gdalType;
		double adfGeoTransform[6]; // 地理转换参数
		const char* projection_info;  // 投影信息，WKT格式的
		double dfNoDataValue;         // 无效值
		int bGotNoDataValue = FALSE;

	};
	void resample(std::string dstPath, double targetResolution, InterpolationType interType = InterpolationType::BILINEAR);
	void toMat(cv::Mat& matImg);
	int reProjection(const char* outputFile);

	// 待完成......
	void splitRegion(std::string dstPath, std::string shpFile);
	GDALDataset* reprojectShpToRaster(OGRLayer* poLayer);  // shp--->raster
	
public:
	TiffInfo tifinfos;   // 存储该张tif的一些基本信息

private:
	void init();
	GDALResampleAlg InterpType2GDALResampleAlg(InterpolationType type);
	void ReprojectShapefile(const char* inputShapefile, const char* outputShapefile, const char* targetWKT);
	// 计算赤道上的1度距离（米）
	double degreesToMetersAtEquator(double degrees) {
		return degrees * 2 * M_PI * EARTH_RADIUS / 360.0;
	}
	// 根据纬度计算实际的1度的经度距离（米）
	double degreesToMetersAtLatitude(double degrees, double latitude) {
		return degrees * 2 * M_PI * EARTH_RADIUS * cos(latitude * M_PI / 180.0) / 360.0;
	}
	// 计算目标分辨率所需的放大倍数
	double calculateScaleFactor(GDALDataset* poDataset, double targetResolutionMeters);
	// 从影像中获取中心纬度
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


#include "MyTiffFile.h"


MyTiffFile::MyTiffFile(std::string tiffPath)
{
	this->tiffImgPath = tiffPath;
	this->init();
}

MyTiffFile::~MyTiffFile()
{
	if (this->poDataset) {
		GDALClose(this->poDataset);
	}
}

void MyTiffFile::init()
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");  //设置支持中文路径和文件名
	this->poDataset = (GDALDataset*)GDALOpen(this->tiffImgPath.c_str(), GA_ReadOnly);
	if (poDataset == NULL)
	{
		std::cerr << "指定的文件不能打开!" << std::endl;
		return;
	}
	this->poDataset->GetGeoTransform(this->tifinfos.adfGeoTransform);   // 获取六参数信息
	this->tifinfos.nSrcXSize = this->poDataset->GetRasterXSize();
	this->tifinfos.nSrcYSize = this->poDataset->GetRasterYSize();
	this->tifinfos.nBands = this->poDataset->GetRasterCount();
	this->tifinfos.gdalType = this->poDataset->GetRasterBand(1)->GetRasterDataType();
	this->tifinfos.projection_info = this->poDataset->GetProjectionRef();
	this->tifinfos.dfNoDataValue = GDALGetRasterNoDataValue(GDALGetRasterBand(this->poDataset, 1), &this->tifinfos.bGotNoDataValue);
}

void MyTiffFile::resample(std::string dstPath, double targetResolution, InterpolationType interType /*= InterpolationType::BILINEAR*/)
{
	/*
	 * 对该tif进行重采样
	 * @param dstPath: 目标文件的路径
	 * @param targetResolution: 目标分辨率, 必须以m为单位给出
	 * @param interType: 重采样算法名称，默认是双线性插值
	 */
	 // 进入算法之前，已经自动获取了该图像的坐标系和一些参数信息，因此可以直接用，重投影的话，坐标系是不用变的
	std::cout << "开始重采样，目标分辨率为" << targetResolution << " m......\n";
	// 计算放大倍数
	double scaleFactor = calculateScaleFactor(poDataset, targetResolution);
	// 计算新的分辨率
	double newResX = this->tifinfos.adfGeoTransform[1] / scaleFactor;
	double newResY = this->tifinfos.adfGeoTransform[5] / scaleFactor; 

	// 新影像的宽度和高度
	int nWidth = static_cast<int>(this->tifinfos.nSrcXSize * scaleFactor);
	int nHeight = static_cast<int>(this->tifinfos.nSrcYSize * scaleFactor);

	// 创建新的TIF影像
	GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* poNewDataset = pDriver->Create(dstPath.c_str(), nWidth, nHeight, 1, this->tifinfos.gdalType, nullptr);

	// 设置新的地理变换信息
	double newGeoTransform[6] = {
		tifinfos.adfGeoTransform[0],    // 左上角X坐标
		newResX,                        // 新的像素宽度
		tifinfos.adfGeoTransform[2],    // 旋转参数
		tifinfos.adfGeoTransform[3],    // 左上角Y坐标
		tifinfos.adfGeoTransform[4],    // 旋转参数
		newResY                         // 新的像素高度
	};
	poNewDataset->SetGeoTransform(newGeoTransform);
	poNewDataset->SetProjection(this->tifinfos.projection_info);

	float* pData = new float[nWidth * nHeight];
	// 设置其新的投影与重采样方法
	GDALReprojectImage(static_cast<GDALDatasetH>(poDataset), this->tifinfos.projection_info,
		static_cast<GDALDatasetH>(poNewDataset), this->tifinfos.projection_info,
		InterpType2GDALResampleAlg(interType), 0, 0, GDALTermProgress, nullptr, nullptr);
	
	// 设置输出图像的 NoData 值
	if (this->tifinfos.bGotNoDataValue) {
		for (int i = 1; i <= this->tifinfos.nBands; i++) {
			GDALSetRasterNoDataValue(GDALGetRasterBand(poNewDataset, i), this->tifinfos.dfNoDataValue);
		}
	}
	GDALClose(poNewDataset);
	std::cout << "重采样完成，文件保存在" << dstPath << "中!";
}


void MyTiffFile::toMat(cv::Mat& matImg)
{
	/*********************读写GDAL数据*****************************
	**  将数据读取到Mat中,根据GDAL数据类型映射到OpenCV数据类型
	**  需要注意的是cv::Mat创建时，是(height,width)的格式，与GDAL的(width,height)刚好相反。
	**  同时GADL的起始通道数是1，不是0.
	**  输出的是Mat类型的图像matImg,它的图像类型以及通道数量完全由此类中的GDALDataset* poDataset决定
	***************************************************************/

	if (this->poDataset == nullptr) {
		std::cerr << "数据集未初始化" << std::endl;
		return;
	}

	int nXSize = this->tifinfos.nSrcXSize;
	int nYSize = this->tifinfos.nSrcYSize;
	int nBands = this->tifinfos.nBands;
	GDALDataType gdalType = this->tifinfos.gdalType;

	int cvType;
	switch (gdalType) {
	case GDT_Byte:
		cvType = CV_8U;
		break;
	case GDT_UInt16:
		cvType = CV_16U;
		break;
	case GDT_Int16:
		cvType = CV_16S;
		break;
	case GDT_UInt32:
		cvType = CV_32S;
		break;
	case GDT_Int32:
		cvType = CV_32S;
		break;
	case GDT_Float32:
		cvType = CV_32F;
		break;
	case GDT_Float64:
		cvType = CV_64F;
		break;
	default:
		std::cerr << "不支持的GDAL数据类型" << std::endl;
		return;
	}

	// 创建适当的cv::Mat
	if (nBands == 1) {
		matImg = cv::Mat(nYSize, nXSize, cvType);
	}
	else {
		matImg = cv::Mat(nYSize, nXSize, CV_MAKETYPE(cvType, nBands));
	}

	// 读取每个波段的数据到Mat中
	for (int i = 0; i < nBands; ++i)
	{
		GDALRasterBand* poBand = poDataset->GetRasterBand(i + 1);
		void* pData = matImg.data + i * matImg.step[1]; // 每个波段的数据偏移
		poBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize, pData, nXSize, nYSize, gdalType, 0, 0);
	}
}


GDALResampleAlg MyTiffFile::InterpType2GDALResampleAlg(InterpolationType type)
{
	// 函数将枚举值转换为GDALResampleAlg类型
	switch (type)
	{
	case InterpolationType::NEAREST: return GDALResampleAlg::GRA_NearestNeighbour;
	case InterpolationType::BILINEAR: return GDALResampleAlg::GRA_Bilinear;
	case InterpolationType::CUBIC: return GDALResampleAlg::GRA_Cubic;
	case InterpolationType::CUBICSPLINE: return GDALResampleAlg::GRA_CubicSpline;
	}
}

void MyTiffFile::ReprojectShapefile(const char* inputShapefile, const char* outputShapefile, const char* targetWKT)
{
}


double MyTiffFile::calculateScaleFactor(GDALDataset* poDataset, double targetResolutionMeters)
{
	double geoTransform[6];
	poDataset->GetGeoTransform(geoTransform);

	// 获取影像的原始分辨率（度）
	double originalResolutionDegrees = geoTransform[1];

	// 获取影像中心的纬度
	double centerLatitude = getCenterLatitude(poDataset);

	// 计算该纬度下1度的距离（米）
	double metersPerDegree = degreesToMetersAtLatitude(1.0, centerLatitude);

	// 计算原始分辨率对应的米
	double originalResolutionMeters = originalResolutionDegrees * metersPerDegree;

	// 计算放大倍数（原始米分辨率除以目标米分辨率）
	double scaleFactor = originalResolutionMeters / targetResolutionMeters;

	return scaleFactor;
}

double MyTiffFile::getCenterLatitude(GDALDataset* poDataset)
{
	double topLatitude = this->tifinfos.adfGeoTransform[3];
	double pixelHeight = this->tifinfos.adfGeoTransform[5]; // 像素高度通常为负数
	int height = this->tifinfos.nSrcYSize;

	// 计算中心纬度
	double centerLatitude = topLatitude + pixelHeight * height / 2;
	return centerLatitude;
}

void MyTiffFile::imageReprojection(const char* pszDstFile, const char* pszDstWKT, InterpolationType eResampleMethod, 
	double dResX/* = 0.0*/, double dResY/* = 0.0*/, const char* pszFormat/* = "GTiff"*/)
{
	/*
	 --- 如果目标投影pszDstWKT是一个地理坐标系（如 WGS84），
	 --- 那么 dResX 和 dResY 表示每个像素在经度和纬度方向上的分辨率，单位是度。
	 --- 如果目标投影pszDstWKT是一个投影坐标系（如 UTM），
	 --- 那么 dResX 和 dResY 表示每个像素的宽度和高度在投影坐标中的实际距离，单位是米。
	*/
	std::cout << "\n开始进行重投影...";
	if (pszDstWKT == NULL)
	{
		std::cout << "目标投影信息为空!\n";
		return;
	}

	// 构造通用转换选项, CSLSetNameValue会根据提供的键值对更新或创建新的键值对并返回更新后的数组
	// void* 是一个通用指针类型，可以指向任何类型的数据。
	// GDALCreateGenImgProjTransformer2 返回 void* 是为了隐藏具体的实现细节，并提供灵活性。
	// 后续使用时，不需要显式地操作 void* ，只需要将它传递给 GDAL 的相关函数，GDAL 会处理内部的数据结构。
	char** papszTO = NULL;
	papszTO = CSLSetNameValue(papszTO, "SRC_SRS", this->tifinfos.projection_info);
	papszTO = CSLSetNameValue(papszTO, "DST_SRS", pszDstWKT);
	void* hTransformArg = GDALCreateGenImgProjTransformer2(this->poDataset, NULL, papszTO);
	if (hTransformArg == NULL)
	{
		std::cout << "构造通用转换选项失败!\n";
	}

	/*使用 SuggestedwarpOutput函数计算输出图像四至范围、大小、六参数等信息
	* adfGeoTransform[6]：这是一个数组，用于存储输出图像的仿射变换参数，通常称为 "六参数"。
	* 仿射变换矩阵可以用于将像素坐标（行、列）转换为地理坐标（经度、纬度），其格式如下：
	* adfGeoTransform[0]：左上角像素的X坐标（通常是经度）。
	* adfGeoTransform[1]：每个像素对应的X方向上的像素大小（像素宽度）。
	* adfGeoTransform[2]：行偏移量（旋转或倾斜，如果不旋转则为0）。
	* adfGeoTransform[3]：左上角像素的Y坐标（通常是纬度）。
	* adfGeoTransform[4]：列偏移量（旋转或倾斜，如果不旋转则为0）。
	* adfGeoTransform[5]：每个像素对应的Y方向上的像素大小（像素高度，通常是负数，因为Y坐标向下增加）。
	* adfExtent[4]：这是一个数组，存储输出图像的地理四至范围，分别表示图像的西、东、南、北边界（通常是经度和纬度）。
	* 即 adfExtent[0] 为西边界，adfExtent[1] 为东边界，adfExtent[2] 为南边界，adfExtent[3] 为北边界。
	* nPixels, nLines：这两个变量用于存储输出图像的像素宽度（列数）和高度（行数），即输出图像的分辨率大小。
	*/
	double outGeoTransform[6];
	double adfExtent[4];  // 输出图像的四至范围
	int nCols, nRows;  // 输出图像的列数和行数
	if (GDALSuggestedWarpOutput2(this->poDataset, GDALGenImgProjTransform, hTransformArg, outGeoTransform, &nCols, &nRows, adfExtent, 0) != CE_None)
	{
		std::cout << "获取图像的四至范围和六参数等信息出错!\n" << std::endl;
		return;
	}
	GDALDestroyGenImgProjTransformer(hTransformArg);

	// 下面开始根据用户指定的分辨率来反算输出图像的大小和六参数等信息
	double dResXSize = dResX;
	double dResYSize = dResY;

	if (dResXSize == 0.0 && dResYSize == 0.0 && this->tifinfos.projection_info[0] == 'P' && pszDstWKT[0] == 'P')
	{
		// 如果源图像和结果图像都为投影坐标系统，则其分辨率与原始影像一致
		dResXSize = ABS(this->tifinfos.adfGeoTransform[1]);
		dResYSize = ABS(this->tifinfos.adfGeoTransform[5]);
	}

	/*===========================这部分是必须的，因为当用户指定了分辨率需要重新计算=============================*/
	//如果用户指定了输出图像的分辨率
	if (dResXSize != 0.0 || dResYSize != 0.0)
	{
		//如果只指定了一个，使用自动计算的结果
		if (dResXSize == 0.0) dResXSize = outGeoTransform[1];
		if (dResYSize == 0.0) dResYSize = outGeoTransform[5];
		//确保分辨率符号正确
		if (dResXSize < 0.0)
			dResXSize = -dResXSize;
		if (dResYSize > 0.0)
			dResYSize = -dResYSize;
		//计算输出图像的范围
		double minX = outGeoTransform[0];
		double maxX = outGeoTransform[0] + outGeoTransform[1] * nCols;
		double maxY = outGeoTransform[3];
		double minY = outGeoTransform[3] + outGeoTransform[5] * nRows;
		//按照用户指定的分辨率来计算图像的输出大小以及范围
		nCols = (int)(((maxX - minX) / dResXSize) + 0.5);
		nRows = (int)(((minY - maxY) / dResYSize) + 0.5);
		outGeoTransform[0] = minX;
		outGeoTransform[3] = maxY;
		outGeoTransform[1] = dResXSize;
		outGeoTransform[5] = dResYSize;
	}
	/*===========================创建输出图像并进行重投影=============================*/
	// 创建输出图像所需要的信息有:
	/* 图像的行数nLines, 列数nPixels, 波段数，输出图像的数据类型， 输出图像的投影信息，输出图像的六参数信息
	 * 所以最关键的是获取输出图像的 ***行数和列数，输出图像的六参数信息***,上面的代码已经获取到
	 */
	GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* poNewDataset = pDriver->Create(pszDstFile, nCols, nRows, this->tifinfos.nBands, this->tifinfos.gdalType, nullptr);
	poNewDataset->SetProjection(pszDstWKT);
	poNewDataset->SetGeoTransform(outGeoTransform);

	// 设置输出图像的 NoData 值
	if (this->tifinfos.bGotNoDataValue) {
		for (int i = 1; i <= this->tifinfos.nBands; i++) {
			GDALSetRasterNoDataValue(GDALGetRasterBand(poNewDataset, i), this->tifinfos.dfNoDataValue);
		}
	}
	/*====================================Part2 : 该函数用于重投影===============================================*/
	// 函数GDALReprojectImage用于重投影，但是该函数不会创建输出图像，只会对结果图像的信息进行修改，所以在使用的时候必须要先将
	// 输出图像创建好，输出图像里原来的投影信息，元数据，颜色表信息等在执行完该函数之后都会被覆盖
	// 所以在使用这个函数的时候需要先创建好输出图像，因此上面的代码都是用于获取输出图像的大小，分辨率等信息的过程，用于创建hDstDS
	if (GDALReprojectImage(this->poDataset, this->tifinfos.projection_info, poNewDataset, pszDstWKT,
		InterpType2GDALResampleAlg(eResampleMethod), 128, 0.0, GDALTermProgress, NULL, NULL) != CE_None)
	{
		GDALClose(poNewDataset);
		GDALDeleteDataset(pDriver, pszDstFile);//删除结果数据
	}
	std::cout << "重投影完成!\n" << std::endl;
	GDALClose(poNewDataset);
}


double* MyTiffFile::getImageCenterCoordinates()
{
	// 计算图像中心像素的坐标
	int centerX = this->tifinfos.nSrcXSize / 2;
	int centerY = this->tifinfos.nSrcYSize / 2;

	// 通过仿射变换将像素坐标转换为地理坐标
	// 地理坐标X = GT(0) + Pixel * GT(1) + Line * GT(2)
	// 地理坐标Y = GT(3) + Pixel * GT(4) + Line * GT(5)
	double geoX = this->tifinfos.adfGeoTransform[0] + centerX * this->tifinfos.adfGeoTransform[1] + 
		centerY * this->tifinfos.adfGeoTransform[2];
	double geoY = this->tifinfos.adfGeoTransform[3] + centerX * this->tifinfos.adfGeoTransform[4] + 
		centerY * this->tifinfos.adfGeoTransform[5];

	double lon_lat[2] = { geoX, geoY };
	// 输出中心的地理坐标
	// std::cout << "图像中心的地理坐标 (经度, 纬度): (" << geoX << ", " << geoY << ")" << std::endl;
	return lon_lat;
}

std::vector<std::string> MyTiffFile::getUTMAndWGS84WKT(double longitude, double latitude)
{
	// 计算UTM带号
	int zone = static_cast<int>(std::floor((longitude + 180.0) / 6.0)) + 1;

	// 判断是否是北半球
	int isNorth = (latitude >= 0) ? TRUE : FALSE;

	// 创建WGS84地理坐标系的空间参考
	OGRSpatialReference oWGS84SRS;
	oWGS84SRS.SetWellKnownGeogCS("WGS84");

	// 导出WGS84的WKT格式
	char* pszWGS84WKT = nullptr;
	oWGS84SRS.exportToPrettyWkt(&pszWGS84WKT, FALSE);

	// 创建UTM投影坐标系的空间参考
	OGRSpatialReference oUTMSRS;
	oUTMSRS.SetUTM(zone, isNorth); // 设置UTM带
	oUTMSRS.SetWellKnownGeogCS("WGS84"); // UTM基于WGS84地理坐标系

	// 导出UTM的WKT格式
	char* pszUTMWKT = nullptr;
	oUTMSRS.exportToPrettyWkt(&pszUTMWKT, FALSE);

	std::vector<std::string> wkts;
	wkts.push_back(std::string(pszWGS84WKT));
	wkts.push_back(std::string(pszUTMWKT));

	// 输出结果
// 	std::cout << "UTM Zone " << zone << (isNorth ? "N" : "S") << " in WKT format:" << std::endl;
// 	std::cout << pszUTMWKT << std::endl;
// 
// 	std::cout << "WGS84 geographic coordinate system in WKT format:" << std::endl;
// 	std::cout << pszWGS84WKT << std::endl;

	// 释放内存
	CPLFree(pszUTMWKT);
	CPLFree(pszWGS84WKT);
	return wkts;
}

int MyTiffFile::reProjection(const char* outputFile)
{
	// 获取其中心经度和中心纬度
	double* lon_lat = getImageCenterCoordinates();
	// 获取目标坐标系下的Wkt描述
	std::vector<std::string> wkts = getUTMAndWGS84WKT(lon_lat[0], lon_lat[1]);
	// 图像的重投影, 假设投影到UTM坐标系
	imageReprojection(outputFile, wkts[1].c_str(), InterpolationType::BILINEAR);
	return 0;
}


// void MyTiffFile::splitRegion(std::string dstPath, std::string shpFile) {
// 	// 禁用磁盘空间检查
// 	CPLSetConfigOption("CHECK_DISK_FREE_SPACE", "FALSE");
// 
// 	GDALAllRegister();
// 	OGRRegisterAll();
// 
// 	const char* tiffWKT = this->poDataset->GetProjectionRef();
// 	const char* reprojectedShapefile = "D:\\CodeFile\\C++_FILE\\NetCDF_Code\\TIF2NC\\testFiles\\reprojected.shp";
// 	ReprojectShapefile(shpFile.c_str(), reprojectedShapefile, tiffWKT);
// 
// 	GDALDataset* poShapefileDS = (GDALDataset*)GDALOpenEx(reprojectedShapefile, GDAL_OF_VECTOR, NULL, NULL, NULL);
// 	if (poShapefileDS == NULL) {
// 		std::cerr << "Open failed: " << reprojectedShapefile << std::endl;
// 		return;
// 	}
// 
// 	OGRLayer* poLayer = poShapefileDS->GetLayer(0);
// 	OGRFeature* poFeature = poLayer->GetNextFeature();
// 	OGRGeometry* poGeometry = poFeature->GetGeometryRef();
// 
// 	OGRPolygon* poPolygon = (OGRPolygon*)poGeometry;
// 	OGREnvelope envelope;
// 	poPolygon->getEnvelope(&envelope);
// 
// 	double adfGeoTransform[6];
// 	this->poDataset->GetGeoTransform(adfGeoTransform);
// 
// 	// Debug information
// 	std::cout << "GeoTransform: "
// 		<< adfGeoTransform[0] << ", "
// 		<< adfGeoTransform[1] << ", "
// 		<< adfGeoTransform[2] << ", "
// 		<< adfGeoTransform[3] << ", "
// 		<< adfGeoTransform[4] << ", "
// 		<< adfGeoTransform[5] << std::endl;
// 
// 	int xMin = static_cast<int>((envelope.MinX - adfGeoTransform[0]) / adfGeoTransform[1]);
// 	int xMax = static_cast<int>((envelope.MaxX - adfGeoTransform[0]) / adfGeoTransform[1]);
// 	int yMin = static_cast<int>((envelope.MaxY - adfGeoTransform[3]) / adfGeoTransform[5]);
// 	int yMax = static_cast<int>((envelope.MinY - adfGeoTransform[3]) / adfGeoTransform[5]);
// 
// 	// Ensure correct row and column calculation
// 	if (yMin > yMax) std::swap(yMin, yMax);
// 	if (xMin > xMax) std::swap(xMin, xMax);
// 
// 	// Validate the calculated dimensions
// 	if (xMin < 0) xMin = 0;
// 	if (xMax > this->poDataset->GetRasterXSize()) xMax = this->poDataset->GetRasterXSize();
// 	if (yMin < 0) yMin = 0;
// 	if (yMax > this->poDataset->GetRasterYSize()) yMax = this->poDataset->GetRasterYSize();
// 
// 	int nCols = xMax - xMin;
// 	int nRows = yMax - yMin;
// 
// 	// Correct yMin and yMax if needed
// 	if (nRows < 0) {
// 		std::cerr << "Adjusted yMin and yMax to fix negative row count" << std::endl;
// 		std::swap(yMin, yMax);
// 		nRows = yMax - yMin;
// 	}
// 
// 	// Debug information
// 	std::cout << "xMin: " << xMin << ", xMax: " << xMax << ", yMin: " << yMin << ", yMax: " << yMax << std::endl;
// 	std::cout << "nCols: " << nCols << ", nRows: " << nRows << std::endl;
// 
// 	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
// 	if (poDriver == NULL) {
// 		std::cerr << "GTiff driver not available." << std::endl;
// 		return;
// 	}
// 
// 	GDALDataset* poOutputDS = poDriver->Create(dstPath.c_str(), nCols, nRows, this->poDataset->GetRasterCount(), GDT_Byte, NULL);
// 	if (poOutputDS == NULL) {
// 		std::cerr << "Creation of output TIFF failed: " << dstPath << std::endl;
// 		return;
// 	}
// 
// 	double adfOutputGeoTransform[6] = {
// 		adfGeoTransform[0] + xMin * adfGeoTransform[1],
// 		adfGeoTransform[1], 0,
// 		adfGeoTransform[3] + yMin * adfGeoTransform[5],
// 		0, adfGeoTransform[5]
// 	};
// 	poOutputDS->SetGeoTransform(adfOutputGeoTransform);
// 	poOutputDS->SetProjection(tiffWKT);
// 
// 	for (int i = 1; i <= this->poDataset->GetRasterCount(); i++) {
// 		GDALRasterBand* poBand = this->poDataset->GetRasterBand(i);
// 		GDALRasterBand* poOutputBand = poOutputDS->GetRasterBand(i);
// 		std::vector<float> buffer(nCols * nRows);
// 		CPLErr err = poBand->RasterIO(GF_Read, xMin, yMin, nCols, nRows, buffer.data(), nCols, nRows, GDT_Float32, 0, 0);
// 		if (err != CE_None) {
// 			std::cerr << "RasterIO read failed." << std::endl;
// 			return;
// 		}
// 		err = poOutputBand->RasterIO(GF_Write, 0, 0, nCols, nRows, buffer.data(), nCols, nRows, GDT_Float32, 0, 0);
// 		if (err != CE_None) {
// 			std::cerr << "RasterIO write failed." << std::endl;
// 			return;
// 		}
// 	}
// 
// 	GDALClose(poShapefileDS);
// 	GDALClose(poOutputDS);
// 
// 	std::cout << "Clipping completed successfully." << std::endl;
// }
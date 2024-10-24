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
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");  //����֧������·�����ļ���
	this->poDataset = (GDALDataset*)GDALOpen(this->tiffImgPath.c_str(), GA_ReadOnly);
	if (poDataset == NULL)
	{
		std::cerr << "ָ�����ļ����ܴ�!" << std::endl;
		return;
	}
	this->poDataset->GetGeoTransform(this->tifinfos.adfGeoTransform);   // ��ȡ��������Ϣ
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
	 * �Ը�tif�����ز���
	 * @param dstPath: Ŀ���ļ���·��
	 * @param targetResolution: Ŀ��ֱ���, ������mΪ��λ����
	 * @param interType: �ز����㷨���ƣ�Ĭ����˫���Բ�ֵ
	 */
	 // �����㷨֮ǰ���Ѿ��Զ���ȡ�˸�ͼ�������ϵ��һЩ������Ϣ����˿���ֱ���ã���ͶӰ�Ļ�������ϵ�ǲ��ñ��
	std::cout << "��ʼ�ز�����Ŀ��ֱ���Ϊ" << targetResolution << " m......\n";
	// ����Ŵ���
	double scaleFactor = calculateScaleFactor(poDataset, targetResolution);
	// �����µķֱ���
	double newResX = this->tifinfos.adfGeoTransform[1] / scaleFactor;
	double newResY = this->tifinfos.adfGeoTransform[5] / scaleFactor; 

	// ��Ӱ��Ŀ�Ⱥ͸߶�
	int nWidth = static_cast<int>(this->tifinfos.nSrcXSize * scaleFactor);
	int nHeight = static_cast<int>(this->tifinfos.nSrcYSize * scaleFactor);

	// �����µ�TIFӰ��
	GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* poNewDataset = pDriver->Create(dstPath.c_str(), nWidth, nHeight, 1, this->tifinfos.gdalType, nullptr);

	// �����µĵ���任��Ϣ
	double newGeoTransform[6] = {
		tifinfos.adfGeoTransform[0],    // ���Ͻ�X����
		newResX,                        // �µ����ؿ��
		tifinfos.adfGeoTransform[2],    // ��ת����
		tifinfos.adfGeoTransform[3],    // ���Ͻ�Y����
		tifinfos.adfGeoTransform[4],    // ��ת����
		newResY                         // �µ����ظ߶�
	};
	poNewDataset->SetGeoTransform(newGeoTransform);
	poNewDataset->SetProjection(this->tifinfos.projection_info);

	float* pData = new float[nWidth * nHeight];
	// �������µ�ͶӰ���ز�������
	GDALReprojectImage(static_cast<GDALDatasetH>(poDataset), this->tifinfos.projection_info,
		static_cast<GDALDatasetH>(poNewDataset), this->tifinfos.projection_info,
		InterpType2GDALResampleAlg(interType), 0, 0, GDALTermProgress, nullptr, nullptr);
	
	// �������ͼ��� NoData ֵ
	if (this->tifinfos.bGotNoDataValue) {
		for (int i = 1; i <= this->tifinfos.nBands; i++) {
			GDALSetRasterNoDataValue(GDALGetRasterBand(poNewDataset, i), this->tifinfos.dfNoDataValue);
		}
	}
	GDALClose(poNewDataset);
	std::cout << "�ز�����ɣ��ļ�������" << dstPath << "��!";
}


void MyTiffFile::toMat(cv::Mat& matImg)
{
	/*********************��дGDAL����*****************************
	**  �����ݶ�ȡ��Mat��,����GDAL��������ӳ�䵽OpenCV��������
	**  ��Ҫע�����cv::Mat����ʱ����(height,width)�ĸ�ʽ����GDAL��(width,height)�պ��෴��
	**  ͬʱGADL����ʼͨ������1������0.
	**  �������Mat���͵�ͼ��matImg,����ͼ�������Լ�ͨ��������ȫ�ɴ����е�GDALDataset* poDataset����
	***************************************************************/

	if (this->poDataset == nullptr) {
		std::cerr << "���ݼ�δ��ʼ��" << std::endl;
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
		std::cerr << "��֧�ֵ�GDAL��������" << std::endl;
		return;
	}

	// �����ʵ���cv::Mat
	if (nBands == 1) {
		matImg = cv::Mat(nYSize, nXSize, cvType);
	}
	else {
		matImg = cv::Mat(nYSize, nXSize, CV_MAKETYPE(cvType, nBands));
	}

	// ��ȡÿ�����ε����ݵ�Mat��
	for (int i = 0; i < nBands; ++i)
	{
		GDALRasterBand* poBand = poDataset->GetRasterBand(i + 1);
		void* pData = matImg.data + i * matImg.step[1]; // ÿ�����ε�����ƫ��
		poBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize, pData, nXSize, nYSize, gdalType, 0, 0);
	}
}


GDALResampleAlg MyTiffFile::InterpType2GDALResampleAlg(InterpolationType type)
{
	// ������ö��ֵת��ΪGDALResampleAlg����
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

	// ��ȡӰ���ԭʼ�ֱ��ʣ��ȣ�
	double originalResolutionDegrees = geoTransform[1];

	// ��ȡӰ�����ĵ�γ��
	double centerLatitude = getCenterLatitude(poDataset);

	// �����γ����1�ȵľ��루�ף�
	double metersPerDegree = degreesToMetersAtLatitude(1.0, centerLatitude);

	// ����ԭʼ�ֱ��ʶ�Ӧ����
	double originalResolutionMeters = originalResolutionDegrees * metersPerDegree;

	// ����Ŵ�����ԭʼ�׷ֱ��ʳ���Ŀ���׷ֱ��ʣ�
	double scaleFactor = originalResolutionMeters / targetResolutionMeters;

	return scaleFactor;
}

double MyTiffFile::getCenterLatitude(GDALDataset* poDataset)
{
	double topLatitude = this->tifinfos.adfGeoTransform[3];
	double pixelHeight = this->tifinfos.adfGeoTransform[5]; // ���ظ߶�ͨ��Ϊ����
	int height = this->tifinfos.nSrcYSize;

	// ��������γ��
	double centerLatitude = topLatitude + pixelHeight * height / 2;
	return centerLatitude;
}

void MyTiffFile::imageReprojection(const char* pszDstFile, const char* pszDstWKT, InterpolationType eResampleMethod, 
	double dResX/* = 0.0*/, double dResY/* = 0.0*/, const char* pszFormat/* = "GTiff"*/)
{
	/*
	 --- ���Ŀ��ͶӰpszDstWKT��һ����������ϵ���� WGS84����
	 --- ��ô dResX �� dResY ��ʾÿ�������ھ��Ⱥ�γ�ȷ����ϵķֱ��ʣ���λ�Ƕȡ�
	 --- ���Ŀ��ͶӰpszDstWKT��һ��ͶӰ����ϵ���� UTM����
	 --- ��ô dResX �� dResY ��ʾÿ�����صĿ�Ⱥ͸߶���ͶӰ�����е�ʵ�ʾ��룬��λ���ס�
	*/
	std::cout << "\n��ʼ������ͶӰ...";
	if (pszDstWKT == NULL)
	{
		std::cout << "Ŀ��ͶӰ��ϢΪ��!\n";
		return;
	}

	// ����ͨ��ת��ѡ��, CSLSetNameValue������ṩ�ļ�ֵ�Ը��»򴴽��µļ�ֵ�Բ����ظ��º������
	// void* ��һ��ͨ��ָ�����ͣ�����ָ���κ����͵����ݡ�
	// GDALCreateGenImgProjTransformer2 ���� void* ��Ϊ�����ؾ����ʵ��ϸ�ڣ����ṩ����ԡ�
	// ����ʹ��ʱ������Ҫ��ʽ�ز��� void* ��ֻ��Ҫ�������ݸ� GDAL ����غ�����GDAL �ᴦ���ڲ������ݽṹ��
	char** papszTO = NULL;
	papszTO = CSLSetNameValue(papszTO, "SRC_SRS", this->tifinfos.projection_info);
	papszTO = CSLSetNameValue(papszTO, "DST_SRS", pszDstWKT);
	void* hTransformArg = GDALCreateGenImgProjTransformer2(this->poDataset, NULL, papszTO);
	if (hTransformArg == NULL)
	{
		std::cout << "����ͨ��ת��ѡ��ʧ��!\n";
	}

	/*ʹ�� SuggestedwarpOutput�����������ͼ��������Χ����С������������Ϣ
	* adfGeoTransform[6]������һ�����飬���ڴ洢���ͼ��ķ���任������ͨ����Ϊ "������"��
	* ����任����������ڽ��������꣨�С��У�ת��Ϊ�������꣨���ȡ�γ�ȣ������ʽ���£�
	* adfGeoTransform[0]�����Ͻ����ص�X���꣨ͨ���Ǿ��ȣ���
	* adfGeoTransform[1]��ÿ�����ض�Ӧ��X�����ϵ����ش�С�����ؿ�ȣ���
	* adfGeoTransform[2]����ƫ��������ת����б���������ת��Ϊ0����
	* adfGeoTransform[3]�����Ͻ����ص�Y���꣨ͨ����γ�ȣ���
	* adfGeoTransform[4]����ƫ��������ת����б���������ת��Ϊ0����
	* adfGeoTransform[5]��ÿ�����ض�Ӧ��Y�����ϵ����ش�С�����ظ߶ȣ�ͨ���Ǹ�������ΪY�����������ӣ���
	* adfExtent[4]������һ�����飬�洢���ͼ��ĵ���������Χ���ֱ��ʾͼ������������ϡ����߽磨ͨ���Ǿ��Ⱥ�γ�ȣ���
	* �� adfExtent[0] Ϊ���߽磬adfExtent[1] Ϊ���߽磬adfExtent[2] Ϊ�ϱ߽磬adfExtent[3] Ϊ���߽硣
	* nPixels, nLines���������������ڴ洢���ͼ������ؿ�ȣ��������͸߶ȣ��������������ͼ��ķֱ��ʴ�С��
	*/
	double outGeoTransform[6];
	double adfExtent[4];  // ���ͼ���������Χ
	int nCols, nRows;  // ���ͼ�������������
	if (GDALSuggestedWarpOutput2(this->poDataset, GDALGenImgProjTransform, hTransformArg, outGeoTransform, &nCols, &nRows, adfExtent, 0) != CE_None)
	{
		std::cout << "��ȡͼ���������Χ������������Ϣ����!\n" << std::endl;
		return;
	}
	GDALDestroyGenImgProjTransformer(hTransformArg);

	// ���濪ʼ�����û�ָ���ķֱ������������ͼ��Ĵ�С������������Ϣ
	double dResXSize = dResX;
	double dResYSize = dResY;

	if (dResXSize == 0.0 && dResYSize == 0.0 && this->tifinfos.projection_info[0] == 'P' && pszDstWKT[0] == 'P')
	{
		// ���Դͼ��ͽ��ͼ��ΪͶӰ����ϵͳ������ֱ�����ԭʼӰ��һ��
		dResXSize = ABS(this->tifinfos.adfGeoTransform[1]);
		dResYSize = ABS(this->tifinfos.adfGeoTransform[5]);
	}

	/*===========================�ⲿ���Ǳ���ģ���Ϊ���û�ָ���˷ֱ�����Ҫ���¼���=============================*/
	//����û�ָ�������ͼ��ķֱ���
	if (dResXSize != 0.0 || dResYSize != 0.0)
	{
		//���ָֻ����һ����ʹ���Զ�����Ľ��
		if (dResXSize == 0.0) dResXSize = outGeoTransform[1];
		if (dResYSize == 0.0) dResYSize = outGeoTransform[5];
		//ȷ���ֱ��ʷ�����ȷ
		if (dResXSize < 0.0)
			dResXSize = -dResXSize;
		if (dResYSize > 0.0)
			dResYSize = -dResYSize;
		//�������ͼ��ķ�Χ
		double minX = outGeoTransform[0];
		double maxX = outGeoTransform[0] + outGeoTransform[1] * nCols;
		double maxY = outGeoTransform[3];
		double minY = outGeoTransform[3] + outGeoTransform[5] * nRows;
		//�����û�ָ���ķֱ���������ͼ��������С�Լ���Χ
		nCols = (int)(((maxX - minX) / dResXSize) + 0.5);
		nRows = (int)(((minY - maxY) / dResYSize) + 0.5);
		outGeoTransform[0] = minX;
		outGeoTransform[3] = maxY;
		outGeoTransform[1] = dResXSize;
		outGeoTransform[5] = dResYSize;
	}
	/*===========================�������ͼ�񲢽�����ͶӰ=============================*/
	// �������ͼ������Ҫ����Ϣ��:
	/* ͼ�������nLines, ����nPixels, �����������ͼ����������ͣ� ���ͼ���ͶӰ��Ϣ�����ͼ�����������Ϣ
	 * ������ؼ����ǻ�ȡ���ͼ��� ***���������������ͼ�����������Ϣ***,����Ĵ����Ѿ���ȡ��
	 */
	GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* poNewDataset = pDriver->Create(pszDstFile, nCols, nRows, this->tifinfos.nBands, this->tifinfos.gdalType, nullptr);
	poNewDataset->SetProjection(pszDstWKT);
	poNewDataset->SetGeoTransform(outGeoTransform);

	// �������ͼ��� NoData ֵ
	if (this->tifinfos.bGotNoDataValue) {
		for (int i = 1; i <= this->tifinfos.nBands; i++) {
			GDALSetRasterNoDataValue(GDALGetRasterBand(poNewDataset, i), this->tifinfos.dfNoDataValue);
		}
	}
	/*====================================Part2 : �ú���������ͶӰ===============================================*/
	// ����GDALReprojectImage������ͶӰ�����Ǹú������ᴴ�����ͼ��ֻ��Խ��ͼ�����Ϣ�����޸ģ�������ʹ�õ�ʱ�����Ҫ�Ƚ�
	// ���ͼ�񴴽��ã����ͼ����ԭ����ͶӰ��Ϣ��Ԫ���ݣ���ɫ����Ϣ����ִ����ú���֮�󶼻ᱻ����
	// ������ʹ�����������ʱ����Ҫ�ȴ��������ͼ���������Ĵ��붼�����ڻ�ȡ���ͼ��Ĵ�С���ֱ��ʵ���Ϣ�Ĺ��̣����ڴ���hDstDS
	if (GDALReprojectImage(this->poDataset, this->tifinfos.projection_info, poNewDataset, pszDstWKT,
		InterpType2GDALResampleAlg(eResampleMethod), 128, 0.0, GDALTermProgress, NULL, NULL) != CE_None)
	{
		GDALClose(poNewDataset);
		GDALDeleteDataset(pDriver, pszDstFile);//ɾ���������
	}
	std::cout << "��ͶӰ���!\n" << std::endl;
	GDALClose(poNewDataset);
}


double* MyTiffFile::getImageCenterCoordinates()
{
	// ����ͼ���������ص�����
	int centerX = this->tifinfos.nSrcXSize / 2;
	int centerY = this->tifinfos.nSrcYSize / 2;

	// ͨ������任����������ת��Ϊ��������
	// ��������X = GT(0) + Pixel * GT(1) + Line * GT(2)
	// ��������Y = GT(3) + Pixel * GT(4) + Line * GT(5)
	double geoX = this->tifinfos.adfGeoTransform[0] + centerX * this->tifinfos.adfGeoTransform[1] + 
		centerY * this->tifinfos.adfGeoTransform[2];
	double geoY = this->tifinfos.adfGeoTransform[3] + centerX * this->tifinfos.adfGeoTransform[4] + 
		centerY * this->tifinfos.adfGeoTransform[5];

	double lon_lat[2] = { geoX, geoY };
	// ������ĵĵ�������
	// std::cout << "ͼ�����ĵĵ������� (����, γ��): (" << geoX << ", " << geoY << ")" << std::endl;
	return lon_lat;
}

std::vector<std::string> MyTiffFile::getUTMAndWGS84WKT(double longitude, double latitude)
{
	// ����UTM����
	int zone = static_cast<int>(std::floor((longitude + 180.0) / 6.0)) + 1;

	// �ж��Ƿ��Ǳ�����
	int isNorth = (latitude >= 0) ? TRUE : FALSE;

	// ����WGS84��������ϵ�Ŀռ�ο�
	OGRSpatialReference oWGS84SRS;
	oWGS84SRS.SetWellKnownGeogCS("WGS84");

	// ����WGS84��WKT��ʽ
	char* pszWGS84WKT = nullptr;
	oWGS84SRS.exportToPrettyWkt(&pszWGS84WKT, FALSE);

	// ����UTMͶӰ����ϵ�Ŀռ�ο�
	OGRSpatialReference oUTMSRS;
	oUTMSRS.SetUTM(zone, isNorth); // ����UTM��
	oUTMSRS.SetWellKnownGeogCS("WGS84"); // UTM����WGS84��������ϵ

	// ����UTM��WKT��ʽ
	char* pszUTMWKT = nullptr;
	oUTMSRS.exportToPrettyWkt(&pszUTMWKT, FALSE);

	std::vector<std::string> wkts;
	wkts.push_back(std::string(pszWGS84WKT));
	wkts.push_back(std::string(pszUTMWKT));

	// ������
// 	std::cout << "UTM Zone " << zone << (isNorth ? "N" : "S") << " in WKT format:" << std::endl;
// 	std::cout << pszUTMWKT << std::endl;
// 
// 	std::cout << "WGS84 geographic coordinate system in WKT format:" << std::endl;
// 	std::cout << pszWGS84WKT << std::endl;

	// �ͷ��ڴ�
	CPLFree(pszUTMWKT);
	CPLFree(pszWGS84WKT);
	return wkts;
}

int MyTiffFile::reProjection(const char* outputFile)
{
	// ��ȡ�����ľ��Ⱥ�����γ��
	double* lon_lat = getImageCenterCoordinates();
	// ��ȡĿ������ϵ�µ�Wkt����
	std::vector<std::string> wkts = getUTMAndWGS84WKT(lon_lat[0], lon_lat[1]);
	// ͼ�����ͶӰ, ����ͶӰ��UTM����ϵ
	imageReprojection(outputFile, wkts[1].c_str(), InterpolationType::BILINEAR);
	return 0;
}


// void MyTiffFile::splitRegion(std::string dstPath, std::string shpFile) {
// 	// ���ô��̿ռ���
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
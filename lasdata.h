#ifndef LASDATA_H_INCLUDED
#define LASDATA_H_INCLUDED

#include <string>
#include <vector>

struct float3
{
    float x, y, z;
};

class PointCloud
{
public:
    PointCloud(const std::string &path);

    win_t getVertsCount();
    float3 *getVertsData();

private:

    std::vector<float3> verts;

    #pragma pack(1)
    struct __attribute__ ((packed)) Header
    {
        char magic[4];
        wint16_t fileSourceID;
        wint16_t globalEncoding;
        wint16_t guidData1;
        wint16_t guidData2;
        wint16_t guidData3;
        wint8_t guidData4[8];
        wint8_t versionMaj, versionMin;
        char systemIdentifier[32];
        char genSoftware[32];
        wint16_t creationDay, creationYear;
        wint16_t headerSize;
        wint32_t pointDataOffset;
        wint32_t numVarLenRecords;
        wint8_t pointDataRecordFormat;
        wint16_t pointDataRecordLen;
        wint32_t numberOfPoints;
        wint32_t numPointsByReturn[5];
        double scaleX, scaleY, scaleZ;
        double offX, offY, offZ;
        double minX, minY, minZ;
        double maxX, maxY, maxZ;
    };

    struct __attribute__ ((packed)) PointRecord1
    {
        win32_t x, y, z;
        win16_t intensity;
        win8_t flags;
        win8_t classification;
        win8_t scanAngleRank;
        win8_t userData;
        win16_t pointSourceID;
        double gpsTime;
    };

    void read(const std::string &path);
};

#endif // LASDATA_H_INCLUDED

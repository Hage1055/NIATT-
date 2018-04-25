#include <liblas/liblas.hpp>
#include <fstream>  // std::ifstream
#include <iostream> // std::cout

std::ifstream ifs;
ifs.open("file.las", std::ios::in | std::ios::binary); //Saddam...Change the file pathway to find where the file is for you

liblas::ReaderFactory f;
liblas::Reader reader = f.CreateWithStream(ifs); //Can also use reader, but this accounts for file compression as well

liblas::Header const& header = reader.GetHeader();

std::cout << "Compressed: " << (header.Compressed() == true) ? "true":"false";
std::cout << "Signature: " << header.GetFileSignature() << '\n';
std::cout << "Points count: " << header.GetPointRecordsCount() << '\n'; //The reader needs to access header data

//reader.ReadPointAt(2);
//liblas::Point const& p = reader.GetPoint(); //This reads a random point. Not necessary to read strips of data

//reader.Seek(10); //This will start the reader at a specified point if necessary

while (reader.ReadNextPoint())
{
    liblas::Point const& p = reader.GetPoint();

    std::cout << p.GetX() << ", " << p.GetY() << ", " << p.GetZ() << "\n";

    float xcore[] = p.GetX();
    float ycore[] = p.GetY();
    float zcore[] = p.GetZ(); //Define points as floating point integers

} //Iterate through all points. Use the above code to start at a specified data point

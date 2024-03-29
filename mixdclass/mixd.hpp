#ifndef _MIXD_H_
#define _MIXD_H_

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include <utility>
#include <vector>
#include <algorithm>
#include <limits>

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace mixd
{


class MixdException
{

private:
    std::string _msg;

public:
    MixdException() { }

    MixdException(const std::string & msg) : _msg(msg) { }

    inline const std::string & msg() const { return _msg; }

};



template <typename T>
class MixdFile
{

private:
    std::string _fname;
    size_t _rows, _cols;

    T * _data;

    bool _debug;



    // check if machine is big endian
    inline bool isBigEndian()
    {
        short word = 0x4321;
        if ((* (char*) & word) != 0x21)
            return true;
        else
            return false;
    }

    // swap bytes according to machine's endianness
    inline void swapbytes(char *array, size_t nelem, size_t elsize)
    {
        size_t sizet, sizem, i, j;
        char *bytea, *byteb;
        sizet = elsize;
        sizem = sizet - 1;
        bytea = (char*)malloc(sizet);
        byteb = (char*)malloc(sizet);

        for (i=0; i<nelem; i++)
        {
            memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
            for (j = 0; j < sizet; j++) byteb[j] = bytea[sizem - j];
                memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
        }

        free(bytea);
        free(byteb);
    }


public:

    MixdFile(const std::string & fname, size_t rows, size_t cols=1, bool debug=true)
        :   _fname(fname), _rows(rows), _cols(cols), _debug(debug)
    {
        if(_debug)
        {
            if(_fname == "")
                std::cout << "Allocating memory for scratch file... " << std::flush;
            else
                std::cout << "Allocating memory for file " << _fname << "... " << std::flush;
        }

        if(_rows>0 && _cols>0)
            _data = new T[_rows*_cols];

        else
            throw MixdException("illegal size of MIXD file");

        if(_debug)
            std::cout << "done!" << std::endl;
    }


    ~MixdFile()
    {
        delete[] _data;
    }


    inline T* data() const
    {
        return _data;
    }

    inline void fname(std::string newfname)
    {
        _fname = newfname;
    }

    inline T& at(size_t row, size_t col=0) const
    {
        return _data[ row*_cols + col ];
    }

    inline T& operator() (size_t row, size_t col=0) const
    {
        return at(row,col);
    }

    inline size_t rows() const
    {
        return _rows;
    }

    inline size_t cols() const
    {
        return _cols;
    }


    inline void init(T val=0)
    {
        for(size_t i=0; i<_rows*_cols; i++)
            _data[i] = val;
    }

    inline T max()
    {
        return(*std::max_element(_data, _data + _rows * _cols -1));
    }

    inline T maxcol(size_t iCol)
    {
        T maxValue = std::numeric_limits<T>::min();
        for(size_t i=iCol; i<_rows*_cols; i+=_cols)
        {
            if(_data[i] > maxValue)
            {
                maxValue = _data[i];
            }
        }
        return maxValue;
    }

    inline T mincol(size_t iCol)
    {
        T minValue = std::numeric_limits<T>::max();
        for(size_t i=iCol; i<_rows*_cols; i+=_cols)
        {
            if(_data[i] < minValue)
            {
                minValue = _data[i];
            }
        }
        return minValue;
    }

    void transpose()
    {
        size_t h = _rows;
        _rows = _cols;
        _cols = h;

        if(_rows > 1 && _cols > 1)
        {
            T *temp = new T[_rows * _cols];

            // copy all values
            memcpy(temp, _data, _rows*_cols*sizeof(T));

            for(size_t i=0; i<_rows; i++)
                for(size_t j=0; j<_cols; j++)
                    at(i,j) = temp[j*_rows + i];

            delete[] temp;
        }
    }


    void read(size_t skip=0, bool transpose=false, size_t rowskip=0)
    {
        if(_debug)
        {
            std::cout << "Reading file " << _fname;
            if(transpose) std::cout << " (transposed)... " << std::flush;
            else          std::cout << "... " << std::flush;
        }

        FILE *fid = fopen(_fname.c_str(), "rb");

        if (fid == NULL)
            throw MixdException("could not open file for reading");

        if(fseek(fid, skip*_rows*_cols*sizeof(T)+rowskip*_cols*sizeof(T), SEEK_SET))
            throw MixdException("could not reposition file pointer");

        if(!transpose)
        {
            // read file in one chunk
            size_t itemsread = fread(_data, sizeof(T), _rows*_cols, fid);

            if(itemsread != _rows*_cols)
                throw MixdException("Error reading file!");
        }
        else
        {
            size_t h = _rows;
            _rows = _cols;
            _cols = h;

            // read file one element after the other
            for(size_t i=0; i<_cols; i++)
                for(size_t j=0; j<_rows; j++)
                {
                    size_t itemsread = fread(&(at(j,i)), sizeof(T), 1, fid);

                    if(itemsread != 1)
                        throw MixdException("Error reading file!");
                }
        }

        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        fclose(fid);

        if(_debug)
            std::cout << "done!" << std::endl;
    }


    void write(size_t skip=0, size_t rowskip=0)
    {
        if(_debug)
            std::cout << "Writing file " << _fname << "... " << std::flush;

        FILE *fid = fopen(_fname.c_str(), "wb");

        if (fid == NULL)
            throw MixdException("could not open file for writing");

        if(fseek(fid, skip*_rows*_cols*sizeof(T)+rowskip*_cols*sizeof(T), SEEK_SET))
            throw MixdException("could not reposition file pointer");

        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        size_t itemswritten = fwrite(_data, sizeof(T), _rows*_cols, fid);

        if(itemswritten != _rows*_cols)
            throw MixdException("Error writing file!");

        // swap back for further use!
        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        fclose(fid);

        if(_debug)
            std::cout << "done!" << std::endl;
    }

    void append()
    {
        if(_debug)
            std::cout << "Writing file " << _fname << "... " << std::flush;

        FILE *fid = fopen(_fname.c_str(), "ab");

        if (fid == NULL)
            throw MixdException("could not open file for writing");

        /* if(fseek(fid, 0, SEEK_END)) */
        /*     throw MixdException("could not reposition file pointer"); */

        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        size_t itemswritten = fwrite(_data, sizeof(T), _rows*_cols, fid);

        if(itemswritten != _rows*_cols)
            throw MixdException("Error writing file!");

        // swap back for further use!
        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        fclose(fid);

        if(_debug)
            std::cout << "done!" << std::endl;
    }


#ifdef PARALLEL
    void read(MPI_Comm & comm, size_t skip=0, bool transpose=false, size_t rowskip=0)
    {
        if(_debug)
        {
            std::cout << "Reading file " << _fname;
            if(transpose) std::cout << " (transposed)... " << std::flush;
            else          std::cout << "... " << std::flush;
        }

        MPI_File fid;
        MPI_Status stat;

        int flag;

        flag = MPI_File_open(comm, (char*)_fname.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);

        if (flag)
            throw MixdException("could not open file for reading");

        flag = MPI_File_seek(fid, skip*_rows*_cols*sizeof(T)+rowskip*_cols*sizeof(T), MPI_SEEK_SET);

        if(flag)
            throw MixdException("could not reposition file pointer");

        if(!transpose)
        {
            // read file in one chunk
            flag = MPI_File_read(fid, _data, _rows*_cols*sizeof(T), MPI_CHAR, &stat);

            if(flag)
                throw MixdException("Error reading file!");
        }
        else
        {
            // something is wrong here... don't use!
            throw MixdException("transposed parallel read not implemented yet");

//             size_t h = _rows;
//             _rows = _cols;
//             _cols = h;
//
//             // read file one element after the other
//             for(size_t i=0; i<_cols; i++)
//                 for(size_t j=0; j<_rows; j++)
//                 {
//                     flag = MPI_File_read(fid, &(at(j,i)), sizeof(T), MPI_CHAR, &stat);
//
//                     if(flag)
//                         throw MixdException("Error reading file!");
//                 }
        }

        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        MPI_File_close(&fid);

        if(_debug)
            std::cout << "done!" << std::endl;
    }


    void write(MPI_Comm & comm, size_t skip=0, size_t rowskip=0)
    {
        if(_debug)
            std::cout << "Writing file " << _fname << "... " << std::flush;

        MPI_File fid;
        MPI_Status stat;

        int flag;

        flag = MPI_File_open(comm, (char*)_fname.c_str(), (MPI_MODE_WRONLY|MPI_MODE_CREATE), MPI_INFO_NULL, &fid);

        if (flag)
            throw MixdException("could not open file for writing");

        flag = MPI_File_seek(fid, skip*_rows*_cols*sizeof(T)+rowskip*_cols*sizeof(T), MPI_SEEK_SET);

        if(flag)
            throw MixdException("could not reposition file pointer");


        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        flag = MPI_File_write(fid, _data, _rows*_cols*sizeof(T), MPI_CHAR, &stat);

        if(flag)
            throw MixdException("Error writing file!");

        // swap back for further use!
        if(!isBigEndian())
            swapbytes((char*)_data, _rows*_cols, sizeof(T));

        MPI_File_close(&fid);

        if(_debug)
            std::cout << "done!" << std::endl;
    }
#endif


//     This is for some reason unsafe... not sure why exactly.
//     template< typename CType >
//     void convert(CType *conv) const
//     {
//         for(size_t i=0; i<_rows*_cols; i++)
//             conv[i] = (CType) _data[i];
//     }

    void toFloat(float *conv) const
    {
        for(size_t i=0; i<_rows*_cols; i++)
            conv[i] = (float) _data[i];
    }

};


//====================================================================
// here follow some convenience functions for reading ASCII files...
//====================================================================


// convert numeric value to string
template <typename T>
inline std::string str(const T & val)
{
    std::stringstream stream;
    stream << val;
    return stream.str();
}

// file name suffix
inline std::string suffix(const std::string & s)
{
    return s.substr(s.find_last_of('.') + 1);
}

// left trim:
// delete white spaces and tabs at the beginning of s
// and return a new string object
inline std::string ltrim(const std::string & s)
{
    size_t beg = s.find_first_not_of(" \t");

    if(beg == std::string::npos)
        return "";
    else
        return s.substr(beg);
}

inline bool valid(const std::string & s)
{
    return (s.length()>0 && s.at(0)!='#' && s.at(0)!='%');
}

// retrieve the next non-comment line from file stream
inline std::string nextValidLine(std::ifstream & file)
{
    std::string line("");

    while(std::getline(file,line))
    {
        line = ltrim(line);

        // skip empty and comment lines
        if(valid(line)) return line;
    }

    throw MixdException("ERROR: No more valid lines!");
}


inline std::pair< std::string, std::vector<std::string> > getKeyVal(const std::string & s)
{
    size_t beg_key = 0;
    size_t len_key = s.find_first_of(" \t");

    size_t beg_val = s.find_first_not_of(" \t", beg_key+len_key);
    size_t len_val = s.find_last_not_of(" \t") - beg_val + 1;

    if(beg_key==std::string::npos || len_key==std::string::npos || beg_val==std::string::npos || len_val==std::string::npos)
        throw MixdException("PARSER ERROR: Key-Value pair could not be identified! Bad line:\n" + s);

    std::string key = s.substr(beg_key, len_key);
    std::string val = s.substr(beg_val, len_val);

    std::vector<std::string> vals;

    char *cval = (char*)val.c_str();
    char *pch;

    pch = strtok(cval," \t");
    while (pch != NULL)
    {
        vals.push_back(std::string(pch));
        pch = strtok (NULL, " \t");
    }

    return std::pair< std::string,std::vector<std::string> >( key, vals );
}


inline std::vector<std::string> getValues(const std::string & s)
{
    std::vector<std::string> vals;

    char *cval = (char*)s.c_str();
    char *pch;

    pch = strtok(cval," \t");
    while (pch != NULL)
    {
        vals.push_back(std::string(pch));
        pch = strtok (NULL, " \t");
    }

    return vals;
}


// read ASCII minf file and return number of nodes and elements in nn and ne
inline void readminf(const std::string & fname, long *nn, long *ne)
{
    std::cout << "Reading minf file " << fname << "... " << std::flush;

    if(nn==NULL || ne==NULL)
        throw MixdException("illegal NULL pointers");

    *nn = 0;
    *ne = 0;

    std::ifstream file(fname.c_str(), std::ifstream::in);

    if(!file.is_open())
        throw MixdException("could not open file");

    std::string line;

    while(std::getline(file, line))
    {
        line = ltrim(line);

//        if(nen!=NULL && line.compare(0,3,"nen")==0)
//        {
//            line = ltrim(line.substr(3));
//             std::cout << "nen   | line = " << line << std::endl;
//            *nen = atol(line.c_str());
//        }
//
//        else if(nef!=NULL && line.compare(0,3,"nef")==0)
//        {
//            line = ltrim(line.substr(3));
//             std::cout << "nef   | line = " << line << std::endl;
//            *nef = atol(line.c_str());
//        }

        if(line.compare(0,3,"ne ")==0 || line.compare(0,3,"ne\t")==0)
        {
            line = ltrim(line.substr(2));
//             std::cout << "ne    | line = " << line << std::endl;
            *ne = atol(line.c_str());
        }

        else if(line.compare(0,2,"nn")==0)
        {
            line = ltrim(line.substr(2));
//             std::cout << "nn    | line = " << line << std::endl;
            *nn = atol(line.c_str());
        }

    }

    file.close();

    if(*nn<=0 || *ne<=0)
        throw MixdException("illegal number of nodes or elements: ne = " + str<long>(*ne) + ", nn = " + str<long>(*nn));

    std::cout << "done!" << std::endl;
}

// write ASCII minf file
inline void writeminf(const std::string & fname, long nn, long ne)
{
    std::cout << "Writing minf file " << fname << "... " << std::flush;

    if(nn<=0 || ne<=0)
        throw MixdException("illegal number of nodes or elements");

    std::ofstream file(fname.c_str(), std::ifstream::out);

    if(!file.is_open())
        throw MixdException("could not open file");

    file << "ne" << std::setw(14) << ne << std::endl;
    file << "nn" << std::setw(14) << nn << std::endl;

//    if(nen > 0)
//        file << "nen" << std::setw(13) << nen << std::endl;
//
//    if(nef > 0)
//        file << "nef" << std::setw(13) << nef << std::endl;

    file.close();

    std::cout << "done!" << std::endl;
}


inline void parse_elem_type(const std::string & typestr, int *nsd, int *nen, int *nef=NULL, int *nfn=NULL)
{
    std::string tmp(typestr);

    // remove leading dash (e.g. from command line option)
    if(tmp.at(0) == '-')
        tmp.erase(0,1);

    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

    if(tmp == "tri")
    {
        if(nsd!=NULL) *nsd = 2;
        if(nen!=NULL) *nen = 3;
        if(nef!=NULL) *nef = 3;
        if(nfn!=NULL) *nfn = 2;
    }

    else if(tmp == "tet")
    {
        if(nsd!=NULL) *nsd = 3;
        if(nen!=NULL) *nen = 4;
        if(nef!=NULL) *nef = 4;
        if(nfn!=NULL) *nfn = 3;
    }

    else if(tmp == "qua")
    {
        if(nsd!=NULL) *nsd = 2;
        if(nen!=NULL) *nen = 4;
        if(nef!=NULL) *nef = 4;
        if(nfn!=NULL) *nfn = 2;
    }

    else if(tmp == "hex")
    {
        if(nsd!=NULL) *nsd = 3;
        if(nen!=NULL) *nen = 8;
        if(nef!=NULL) *nef = 6;
        if(nfn!=NULL) *nfn = 4;
    }

    else throw MixdException("Unknown element type: " + typestr);
}


} // namespace mixd

#endif

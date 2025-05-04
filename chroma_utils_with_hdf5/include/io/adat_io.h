// -*- C++ -*-
/*! @file
 * @brief IO support
 */

#ifndef ADAT_IO_H
#define ADAT_IO_H

#include <string>
#include <fstream>
#include <sstream>
#include <xml_array.h>
#include <xml_array2d.h>
#include <xml_array3d.h>
#include "adat/adat_arraynd.h"
#include "adat/adat_arrays.h"
#include "adat_byteorder.h"
#include <vector>

namespace ADATIO 
{
  // namespace composition
//  using XMLArray::Array;


  //--------------------------------------------------------------------------------
  //! Binary reader base class
  class BinaryReader
  {
  public:
    BinaryReader();
    virtual ~BinaryReader();

    //! Return true if some failure occurred in previous IO operation
    virtual bool fail();

    //! Read array of bytes and broadcast to all nodes
    virtual void readArray(char* output, size_t nbytes, size_t nmemb);

    // Overloaded reader functions
    virtual void readDesc(std::string& result);

    //! Read some max number of characters - 1 upto and excluding a newline
    /*! This is the getline function for the underlying stream */
    virtual void read(std::string& result, size_t nbytes);

    virtual void read(char& result);
    virtual void read(int& result);
    virtual void read(unsigned int& result);
    virtual void read(short int& result);
    virtual void read(unsigned short int& result);
    virtual void read(long int& result);
    virtual void read(unsigned long int& result);
    virtual void read(float& result);
    virtual void read(double& result);
    virtual void read(bool& result);

    //! Get the current checksum
    virtual ADATUtil::n_uint32_t getChecksum() const = 0;
    
  protected:
    // The universal data-reader. All the read functions call this
    template< typename T>
    void
    readPrimitive(T& output);

    //! Get the current checksum to modify
    virtual ADATUtil::n_uint32_t& setChecksum() = 0;
  
    // Get the internal ostream
    virtual std::istream& getIstream() = 0;
  };


  // Telephone book of basic primitives
  void readDesc(BinaryReader& bin, std::string& input);
  void read(BinaryReader& bin, std::string& input, size_t maxBytes);
  void read(BinaryReader& bin, char& input);
  void read(BinaryReader& bin, int& input);
  void read(BinaryReader& bin, unsigned int& input);
  void read(BinaryReader& bin, short int& input);
  void read(BinaryReader& bin, unsigned short int& input);
  void read(BinaryReader& bin, long int& input);
  void read(BinaryReader& bin, unsigned long int& input);
  void read(BinaryReader& bin, float& input);
  void read(BinaryReader& bin, double& input);
  void read(BinaryReader& bin, bool& input);

  // Different bindings for same operators
  BinaryReader& operator>>(BinaryReader& bin, char& input);
  BinaryReader& operator>>(BinaryReader& bin, int& input);
  BinaryReader& operator>>(BinaryReader& bin, unsigned int& input);
  BinaryReader& operator>>(BinaryReader& bin, short int& input);
  BinaryReader& operator>>(BinaryReader& bin, unsigned short int& input);
  BinaryReader& operator>>(BinaryReader& bin, long int& input);
  BinaryReader& operator>>(BinaryReader& bin, unsigned long int& input);
  BinaryReader& operator>>(BinaryReader& bin, float& input);
  BinaryReader& operator>>(BinaryReader& bin, double& input);
  BinaryReader& operator>>(BinaryReader& bin, bool& input);

  //! Read a binary XMLArray::Array object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The XMLArray::Array can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, XMLArray::Array<T>& d)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    d.resize(n);

    for(int i=0; i < d.size(); ++i)
      read(bin, d[i]);
  }

  //! Read a binary XMLArray::Array object
  /*!
    This assumes that the number of elements to be read is not written in
    the file, \e i.e. that the data was written with the corresponding write
    code. The number of elements must therefore be supplied by the caller
    \param bin The initialised binary reader
    \param d The data to be filled.
    \param num The number of elements.

    \pre The binary reader must have opened the file.
    \pre The XMLArray::Array must have space for at least \a num elements.  
  */
  template<class T>
  inline
  void read(BinaryReader& bin, XMLArray::Array<T>& d, int num)
  {
    for(int i=0; i < num; ++i)
      read(bin, d[i]);
  }

  //! Read a binary XMLArray::Array2d object
  /*!
    This assumes that the number of elements to be read is not written in
    the file, \e i.e. that the data was written with the corresponding write
    code. The number of elements must therefore be supplied by the caller
    \param bin The initialised binary reader
    \param d The data to be filled.
    \param num1 The first dimension of the array
    \param num2 The second dimension of the array..  

    \pre The binary reader must have opened the file.
    \pre The XMLArray::Array2d must have space for at least \a num elements.  
  */
  template<class T>
  inline
  void read(BinaryReader& bin, XMLArray::Array2d<T>& d, int num1, int num2)
  {
    for(int i=0; i < num2; ++i)
      for(int j=0; j < num1; ++j)
	read(bin, d[j][i]);

  }


  //! Read a binary XMLArray::Array2d element
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The XMLArray::Array2d can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, XMLArray::Array2d<T>& d)
  {
    int n1;
    int n2;
    read(bin, n1);    // the size is always written, even if 0
    read(bin, n2);    // the size is always written, even if 0
    d.resize(n1,n2);
  
    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
      {
	read(bin, d[j][i]);
      }
  }


  //! Read a binary XMLArray::Array3d element
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The XMLArray::Array2d can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, XMLArray::Array3d<T>& d)
  {
    int n1;
    int n2;
    int n3;
    read(bin, n1);    // the size is always written, even if 0
    read(bin, n2);    // the size is always written, even if 0
    read(bin, n3);    // the size is always written, even if 0

    // Destructively resize the array
    d.resize(n3,n2,n1);

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
	for(int k=0; k < d.size3(); ++k)
	{
	  //read(bin, d[k][j][i]);
	  read(bin, d(k, j, i) );
	}
  }


  //! Read a binary ArrayNd element
  template<class T>
  inline
  void read(BinaryReader& bin, ADAT::ArrayNd<T>& d)
  {
    XMLArray::Array<int> siz;
    read(bin, siz); // read the array of the sizes

    d.resize(siz);

    for(int i=0; i < d.numElem(); ++i)
      read(bin, d.getElem(i));
  }


  //! Read a binary std::vector object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The std::vector can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, std::vector<T>& d)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    d.resize(n);

    for(int i=0; i < d.size(); ++i)
      read(bin, d[i]);
  }

  //! Read a binary ADAT::Array1dO object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The ADAT::Array1dO can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, ADAT::Array1dO<T>& d)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    d.resize(n);

    for(int i=1; i <= d.size(); ++i)
      read(bin, d[i]);
  }


  //--------------------------------------------------------------------------------
  //! Binary buffer reader class
  class BinaryBufferReader : public BinaryReader
  {
  public:
    BinaryBufferReader();
    ~BinaryBufferReader();

    //! Construct from a string
    explicit BinaryBufferReader(const std::string& s);

    //! Construct from a string
    void open(const std::string& s);

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Is it empty?
    bool eof() const;
        
    //! Get the current checksum
    ADATUtil::n_uint32_t getChecksum() const {return checksum;}
  
  protected:
    //! Get the current checksum to modify
    ADATUtil::n_uint32_t& setChecksum() {return checksum;}
  
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

  private:
    //! Checksum
    ADATUtil::n_uint32_t checksum;
    std::istringstream f;
  };


  //--------------------------------------------------------------------------------
  //! Binary file reader class
  class BinaryFileReader : public BinaryReader
  {
  public:
    BinaryFileReader();
    ~BinaryFileReader();
    //! Open from a file
    explicit BinaryFileReader(const std::string& p);

    bool is_open();
    void open(const std::string& p);
    void close();

    //! Seek offset from end
    void seekEnd(long int off);

    //! Rewind 
    void rewind();

    //! Get the current checksum
    ADATUtil::n_uint32_t getChecksum() const {return checksum;}
  
  protected:
    //! Get the current checksum to modify
    ADATUtil::n_uint32_t& setChecksum() {return checksum;}
  
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

  private:
    //! Checksum
    ADATUtil::n_uint32_t checksum;

    std::ifstream f;
  };


  //-------------------------------------------------------------------------------------
  //!  Binary output base class
  /*!
    This class is used to write data to a binary file. The data in the file
    is big-endian. If the host nachine is little-endian, the data
    is byte-swapped.
    
    Files need to be opened before any of the write methods are used  
    
    The write methods are also wrapped by externally defined functions
    and << operators,   
  */
  class BinaryWriter
  {
  public:
    BinaryWriter();

    /*!
      Closes the last file opened
    */
    virtual ~BinaryWriter();

    //! Flushes the buffer
    virtual void flush() = 0;

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in previous IO operation
    */
    virtual bool fail();

    //! Write data from the primary node.
    /*!
      \param output The data to write
      \param nbytes The size in bytes of each datum
      \param The number of data.
    */
    virtual void writeArray(const char* output, size_t nbytes, size_t nmemb);

    // Overloaded Writer Functions

    //! Writes a fixed length of characters like an array
    virtual void writeDesc(const std::string& output);

    /*!
      A newline is appended to the written string.
    */
    virtual void write(const std::string& output);
    /*!
      A newline is appended to the written string.
    */
    virtual void write(const char* output);
    virtual void write(const char& output);
    virtual void write(const int& output);
    virtual void write(const unsigned int& output);
    virtual void write(const short int& output);
    virtual void write(const unsigned short int& output);
    virtual void write(const long int& output);
    virtual void write(const unsigned long int& output);
    virtual void write(const float& output);
    virtual void write(const double& output);
    virtual void write(const bool& output);

    //! Get the current checksum
    virtual ADATUtil::n_uint32_t getChecksum() const = 0;
  
  protected:

    //! The universal data-writer.
    /*!
      All the write functions call this.
      \param output The location of the datum to be written.
    */
    template< typename T>
    void
    writePrimitive(const T& output);
 
    //! Get the current checksum to modify
    virtual ADATUtil::n_uint32_t& setChecksum() = 0;
  
    //! Get the internal output stream
    virtual std::ostream& getOstream() = 0;
  };


  // Telephone book of basic primitives
  void writeDesc(BinaryWriter& bin, const std::string& output);
  void write(BinaryWriter& bin, const std::string& output);
  void write(BinaryWriter& bin, const char* output);
  void write(BinaryWriter& bin, char output);
  void write(BinaryWriter& bin, int output);
  void write(BinaryWriter& bin, unsigned int output);
  void write(BinaryWriter& bin, short int output);
  void write(BinaryWriter& bin, unsigned short int output);
  void write(BinaryWriter& bin, long int output);
  void write(BinaryWriter& bin, unsigned long int output);
  void write(BinaryWriter& bin, float output);
  void write(BinaryWriter& bin, double output);
  void write(BinaryWriter& bin, bool output);

  // Different bindings for same operators
  BinaryWriter& operator<<(BinaryWriter& bin, const std::string& output);
  BinaryWriter& operator<<(BinaryWriter& bin, const char* output);
  BinaryWriter& operator<<(BinaryWriter& bin, char output);
  BinaryWriter& operator<<(BinaryWriter& bin, int output);
  BinaryWriter& operator<<(BinaryWriter& bin, unsigned int output);
  BinaryWriter& operator<<(BinaryWriter& bin, short int output);
  BinaryWriter& operator<<(BinaryWriter& bin, unsigned short int output);
  BinaryWriter& operator<<(BinaryWriter& bin, long int output);
  BinaryWriter& operator<<(BinaryWriter& bin, unsigned long int output);
  BinaryWriter& operator<<(BinaryWriter& bin, float output);
  BinaryWriter& operator<<(BinaryWriter& bin, double output);
  BinaryWriter& operator<<(BinaryWriter& bin, bool output);

  //! Write all of a binary XMLArray::Array object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const XMLArray::Array<T>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size
    for(int i=0; i < d.size(); ++i)
      write(bin, d[i]);
  }

  //! Write some or all of a binary XMLArray::Array object
  /*!
    This does not write the number of elements to the file.
    \param bin The initialised binary writer
    \param d The data to be filled.
    \param num The number of elements to write.

    \pre The binary writer must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const XMLArray::Array<T>& d, int num)
  {
    for(int i=0; i < num; ++i)
      write(bin, d[i]);
  }


  //! Write a binary XMLArray::Array2d element
  template<class T>
  inline
  void write(BinaryWriter& bin, const XMLArray::Array2d<T>& d)
  {
    int n2 = d.size2();
    int n1 = d.size1();
    write(bin, n2);    // always write the size
    write(bin, n1);    // always write the size

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
      {
	write(bin, d[j][i]);
      }

  }

  //! Write a fixed number of binary XMLArray::Array2d element - no element count written
  template<class T>
  inline
  void write(BinaryWriter& bin, const XMLArray::Array2d<T>& d, int num1, int num2)
  {
    for(int i=0; i < num2; ++i)
      for(int j=0; j < num1; ++j)
	write(bin, d[j][i]);
  }



  //! Write a binary XMLArray::Array3d element
  template<class T>
  inline
  void write(BinaryWriter& bin, const XMLArray::Array3d<T>& d)
  {
    int n3 = d.size3();
    int n2 = d.size2();
    int n1 = d.size1();
    write(bin, n3);    // always write the size
    write(bin, n2);    // always write the size
    write(bin, n1);    // always write the size

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
	for(int k=0; k < d.size3(); ++k)
	  write(bin, d(k,j,i) );
		//write(bin, d[k][j][i]);

  }

  //! Write a fixed number of binary XMLArray::Array3d element - no element count written
  template<class T>
  inline
  void write(BinaryWriter& bin, const XMLArray::Array3d<T>& d, 
	     int num1, int num2, int num3)
  {
    for(int k=0; k < num3 ; ++k)
      for(int j=0; j < num2; ++j)
	for(int i=0; i < num1; ++i)
	  write(bin, d[i][j][k]);

  }


  //! Write a binary multiNd element
  template<class T>
  inline
  void write(BinaryWriter& bin, const ADAT::ArrayNd<T>& d)
  {
    write(bin, d.size()); // write the array of the sizes

    for(int i=0; i < d.numElem(); ++i)
      write(bin, d.getElem(i));
  }


  //! Write all of a binary std::vector object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const std::vector<T>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size
    for(int i=0; i < d.size(); ++i)
      write(bin, d[i]);
  }

  //! Write all of a binary ADAT::Array object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const ADAT::Array1dO<T>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size
    for(int i=1; i <= d.size(); ++i)
      write(bin, d[i]);
  }



  //-------------------------------------------------------------------------------------
  //!  Binary buffer writer output class
  /*!
    This class is used to write data to a binary file. The data in the file
    is big-endian. If the host nachine is little-endian, the data
    is byte-swapped.
    
    Files need to be opened before any of the write methods are used  
    
    The write methods are also wrapped by externally defined functions
    and << operators,   
  */
  class BinaryBufferWriter : public BinaryWriter
  {
  public:
    BinaryBufferWriter();

    //! Construct from a string
    explicit BinaryBufferWriter(const std::string& s);

    //! Closes the buffer
    ~BinaryBufferWriter();

    //! Construct from a string
    void open(const std::string& s);

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Flushes the buffer
    void flush() {}

    //! Get the current checksum
    ADATUtil::n_uint32_t getChecksum() const {return checksum;}
  
  protected:
    //! Get the current checksum
    ADATUtil::n_uint32_t& setChecksum() {return checksum;}
  
    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

  private:
    //! Checksum
    ADATUtil::n_uint32_t checksum;
    std::ostringstream f;
  };


  //-------------------------------------------------------------------------------------
  //!  Binary file output class
  /*!
    This class is used to write data to a binary file. The data in the file
    is big-endian. If the host nachine is little-endian, the data
    is byte-swapped.
    
    Files need to be opened before any of the write methods are used  
    
    The write methods are also wrapped by externally defined functions
    and << operators,   
  */
  class BinaryFileWriter : public BinaryWriter
  {
  public:
    BinaryFileWriter();

    /*!
      Closes the last file opened
    */
    ~BinaryFileWriter();

    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    explicit BinaryFileWriter(const std::string& p);

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();
    
    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    void open(const std::string& p);

    //! Closes the last file opened   
    void close();

    //! Flushes the buffer
    void flush();

    //! Get the current checksum
    ADATUtil::n_uint32_t getChecksum() const {return checksum;}
  
  protected:
    //! Get the current checksum to modify
    ADATUtil::n_uint32_t& setChecksum() {return checksum;}
  
    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

  private:
    //! Checksum
    ADATUtil::n_uint32_t checksum;

    std::ofstream f;
  };


  /*! @} */   // end of group io

}  // namespace ADATIO

#endif

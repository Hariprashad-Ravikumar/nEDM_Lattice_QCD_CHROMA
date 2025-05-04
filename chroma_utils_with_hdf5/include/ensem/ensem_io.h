// -*- C++ -*-
// $Id: ensem_io.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! @file
 * @brief IO support
 */

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace ENSEM {


/*! @defgroup io IO
 *
 * File input and output operations on ENSEM types
 *
 * @{
 */

//--------------------------------------------------------------------------------
//! Text input class
/*!
  This class is used to read data from a text file. Input is done on the
  primary node and all nodes end up with the same data.

  The read methods are also wrapped by externally defined >> operators,
*/

class TextReader
{
public:
    TextReader();
    /*!
      Closes the last file opened
    */
    ~TextReader(); 

    /*!
      Opens a file for reading.
      \param p The name of the file
    */
    explicit TextReader(const std::string& p);

    //! Opens a file for reading.
    /*!
      \param p The name of the file
    */
    void open(const std::string& p);

    //! Closes the last file opened
    void close();

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in the previous IO operation
    */
    bool fail();

    // Readers for builtin types
    void read(std::string& result);
    void read(char& result);
    void read(int& result);
    void read(unsigned int& result);
    void read(short int& result);
    void read(unsigned short int& result);
    void read(long int& result);
    void read(unsigned long int& result);
    void read(float& result);
    void read(double& result);
    void read(bool& result);

    //! Get the internal input stream
    std::istream& getIstream() {return f;}

protected:

    //! The universal data-reader.
    /*!
      All the read functions call this.
      \param output The location to which the datum is read.
    */
    template< typename T>
    void
    readPrimitive(T& output);

private:
  std::ifstream f;
};


// Different bindings for same operators
TextReader& operator>>(TextReader& txt, std::string& input);
TextReader& operator>>(TextReader& txt, char& input);
TextReader& operator>>(TextReader& txt, int& input);
TextReader& operator>>(TextReader& txt, unsigned int& input);
TextReader& operator>>(TextReader& txt, short int& input);
TextReader& operator>>(TextReader& txt, unsigned short int& input);
TextReader& operator>>(TextReader& txt, long int& input);
TextReader& operator>>(TextReader& txt, unsigned long int& input);
TextReader& operator>>(TextReader& txt, float& input);
TextReader& operator>>(TextReader& txt, double& input);
TextReader& operator>>(TextReader& txt, bool& input);


//-----------------------------------------
//! Text output class
/*!
  This class is used to write data to a text file.
  Output is done from the primary node only..
  
  The write methods are also wrapped by externally defined >> operators,
*/

class TextWriter
{
public:
    TextWriter();

    /*!
      Closes the last file opened
    */
    ~TextWriter();

    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    explicit TextWriter(const std::string& p);

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

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in previous IO operation
    */
    bool fail();

    // Overloaded Writer Functions
    void write(const std::string& output);
    void write(const char* output);
    void write(const char& output);
    void write(const int& output);
    void write(const unsigned int& output);
    void write(const short int& output);
    void write(const unsigned short int& output);
    void write(const long int& output);
    void write(const unsigned long int& output);
    void write(const float& output);
    void write(const double& output);
    void write(const bool& output);

    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

protected:
    //! The universal data-writer.
    /*!
      All the write functions call this.
      \param output The location of the datum to be written.
    */
    template< typename T>
    void
    writePrimitive(const T& output);

private:
  std::ofstream f;
};


// Different bindings for same operators
TextWriter& operator<<(TextWriter& txt, const std::string& output);
TextWriter& operator<<(TextWriter& txt, const char* output);
TextWriter& operator<<(TextWriter& txt, char output);
TextWriter& operator<<(TextWriter& txt, int output);
TextWriter& operator<<(TextWriter& txt, unsigned int output);
TextWriter& operator<<(TextWriter& txt, short int output);
TextWriter& operator<<(TextWriter& txt, unsigned short int output);
TextWriter& operator<<(TextWriter& txt, long int output);
TextWriter& operator<<(TextWriter& txt, unsigned long int output);
TextWriter& operator<<(TextWriter& txt, float output);
TextWriter& operator<<(TextWriter& txt, double output);
TextWriter& operator<<(TextWriter& txt, bool output);


/*! @} */   // end of group io

} // namespace ENSEM

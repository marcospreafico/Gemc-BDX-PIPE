// evioFileChannel.hxx

//  Author:  Elliott Wolin, JLab, 18-feb-2010


#ifndef _evioFileChannel_hxx
#define _evioFileChannel_hxx


#include <iostream>
#include <stdint.h>
#include "evioChannel.hxx"
#include "evioUtil.hxx"
#include "evio.h"


using namespace std;


namespace evio {


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 * Implements evioChannel functionality for I/O to and from files.
 * Basically a wrapper around the original evio C library.
 */
class evioFileChannel : public evioChannel {

public:
  evioFileChannel(const string &fileName, const string &mode = "r", int size = 1000000) noexcept(false);
  evioFileChannel(const string &fileName, evioDictionary *dict,
                  const string &mode = "r", int size = 1000000) noexcept(false);
  evioFileChannel(const string &fileName, evioDictionary *dict, const uint32_t *firstEvent,
                  const string &mode = "w", int size = 1000000) noexcept(false);
  virtual ~evioFileChannel(void);


  void open(void) noexcept(false);

  bool read(void) noexcept(false);
  bool read(uint32_t *myEventBuf, int length) noexcept(false);
  bool readAlloc(uint32_t **buffer, uint32_t *bufLen) noexcept(false);
  bool readNoCopy(void) noexcept(false);
  bool readRandom(uint32_t bufferNumber) noexcept(false);

  void write(void) noexcept(false);
  void write(const uint32_t *myEventBuf) noexcept(false);
  void write(const evioChannel &channel) noexcept(false);
  void write(const evioChannel *channel) noexcept(false);
  void write(const evioChannelBufferizable &o) noexcept(false);
  void write(const evioChannelBufferizable *o) noexcept(false);

  void close(void) noexcept(false);

  int ioctl(const string &request, void *argp) noexcept(false);

  const uint32_t *getBuffer(void) const noexcept(false);
  int getBufSize(void) const;
  const uint32_t *getNoCopyBuffer(void) const noexcept(false);
  const uint32_t *getRandomBuffer(void) const noexcept(false);
  void getRandomAccessTable(uint32_t *** const table, uint32_t *len) const noexcept(false);

  string getFileName(void) const;
  string getMode(void) const;
  string getFileXMLDictionary(void) const;


private:
  string filename;            /**<Name of evio file.*/
  string mode;                /**<Open mode, "r" for read, "ra" for random access read,
                                * "w" for write, "a" for append, "s" for splitting while writing.*/
  int handle;                 /**<Internal evio handle.*/
  uint32_t *buf;              /**<Pointer to internal event buffer.*/
  int bufSize;                /**<Size of internal event buffer.*/
  const uint32_t *firstEvent; /**<Pointer first event buffer.*/
  const uint32_t *noCopyBuf;  /**<Pointer to no copy event buffer.*/
  const uint32_t *randomBuf;  /**<Pointer to random read buffer.*/
  string fileXMLDictionary;   /**<XML dictionary in file.*/
  bool createdFileDictionary; /**<true if internally created new dictionary from file.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace evio


#endif

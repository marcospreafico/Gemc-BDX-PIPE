// evioSocketChannel.hxx

//  Author:  Elliott Wolin, JLab, 12-Apr-2012


#ifndef _evioSocketChannel_hxx
#define _evioSocketChannel_hxx


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
 * Implements evioChannel functionality for I/O to and from socket.
 */
class evioSocketChannel : public evioChannel {

public:
  evioSocketChannel(int socFd, const string &mode = "r", int size=100000) noexcept(false);
  evioSocketChannel(int socFd, evioDictionary *dict, const string &mode = "r", int size=100000) noexcept(false);
  virtual ~evioSocketChannel(void);

  void open(void) noexcept(false);

  bool read(void) noexcept(false);
  bool read(uint32_t *myEventBuf, int length) noexcept(false);
  bool readAlloc(uint32_t **buffer, uint32_t *bufLen) noexcept(false);
  bool readNoCopy(void) noexcept(false);

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

  string getMode(void) const;
  int getSocketFD(void) const noexcept(false);
  string getSocketXMLDictionary(void) const;


private:
  int sockFD;                     /**<Socket file descriptor.*/
  string mode;                    /**<Open mode, "r" or "w".*/
  int handle;                     /**<Internal evio handle.*/
  uint32_t *buf;                  /**<Pointer to internal event socket.*/
  int bufSize;                    /**<Size of internal socket.*/
  const uint32_t *noCopyBuf;      /**<Pointer to no copy buffer.*/
  string socketXMLDictionary;     /**<XML dictionary in socket.*/
  bool createdSocketDictionary;   /**<true if created dictionary from socket.*/
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace evio


#endif

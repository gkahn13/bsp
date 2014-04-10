#ifndef SERIAL_HPP
#define SERIAL_HPP

#include <termios.h>
#include <stdio.h>
#include <sys/time.h>
#include <string>

class SerialDevice {
public:
  SerialDevice(const std::string &port, speed_t baud);
  ~SerialDevice();
  // Connect to the device.
  void Connect();
  // Close the connection to the device and restore original settings.
  void Disconnect();
  // Store nbytes in result, assumes memory has been allocated.
  // Return number of bytes read, return -1 if error
  int Read(int nbytes, unsigned char *data);
  // Read from the device until '\r\n' is encountered; return -1 on error
  int Readline(std::string *result);
  // Send a command to the serial port, line should end with '\r\n';
  int Writeline(unsigned char *data, int nbytes);
private:
  void DieIfTrue(bool error, const std::string &msg);

  std::string port_;
  int fd_;
  bool connected_;
  speed_t baud_;
  
  fd_set read_fds_, write_fds_;
  struct termios oldtio_, newtio_;
};

#endif

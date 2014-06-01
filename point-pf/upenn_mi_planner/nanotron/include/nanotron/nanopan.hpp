#ifndef NANO_PAN_HPP
#define NANO_PAN_HPP

#include <ros/ros.h>

#include <nanotron/serial.hpp>

namespace nanotron {
  // See nanoPAN 5375 Development Kit User Guide for description of error
  // status messages.
  enum status_t {
    RG_STAT_DEFAULT = 0x00,
    RG_STAT_NO_ERROR = 0x0F,
    RG_STAT_T1 = 0x01,
    RG_STAT_T2 = 0x02,
    RG_STAT_T3 = 0x04,
    RG_STAT_T4 = 0x08
  };
  
  class NanoPAN {
  public:
    NanoPAN(SerialDevice *s);

    // Get a message from the device and call the appropriate method.  If no
    // message is received within timeout.  Returns -1 on error.
    int ProcessMessage();

    /* Start / stop making range measurements*/
    int Start();
    int Stop();
    /* Set the MAC of the current connected device*/
    int SetMac(int mac);
    /* Set a list of destination MAC addresses */
    int SetDests(const std::vector<int> &macs);
    /* Set the period between distance measurements */
    int SetPeriod(int per);
    /* Print status of radio */
    int Print();
    
    /* Callbacks */
    virtual int Distance(float distance, int src_id, int dest_id, status_t status) {
      return 0;
    }
  private:
    SerialDevice *ser_;
  };

  class Echo : public NanoPAN {
  public:
    Echo(SerialDevice *s) : NanoPAN(s) {}

    int Distance(float distance, int src_id, int dest_id, status_t status) {
      ROS_INFO("distance %f src_id %i dest_id %i status: %i",
               distance, src_id, dest_id, status);
      return 0;
    }
  };
};
#endif

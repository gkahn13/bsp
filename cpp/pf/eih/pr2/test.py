"""
- add base_laser3d
- get depth data from it

- dynamics function uses openrave model

- observation function uses depth data from base_laser3d
"""

import IPython

import openravepy as rave
import time

import rospy

import pr2

def openrave_test():
    env = rave.Environment()
    env.Load('robots/pr2-beta-sim.robot.xml')
    r = env.GetRobots()[0]
    
    #sensors = env.GetAttachedSensors()
    
    env.SetViewer('qtcoin')
    IPython.embed()
    
    ienablesensor = 0
    while True:
        sensors = env.GetSensors()
        for i,sensor in enumerate(sensors):
            if i==ienablesensor:
                sensor.Configure(rave.Sensor.ConfigureCommand.PowerOn)
                sensor.Configure(rave.Sensor.ConfigureCommand.RenderDataOn)
            else:
                sensor.Configure(rave.Sensor.ConfigureCommand.PowerOff)
                sensor.Configure(rave.Sensor.ConfigureCommand.RenderDataOff)
        print 'showing sensor %s, try moving obstacles'%(sensors[ienablesensor].GetName())
        if sensors[ienablesensor].Supports(rave.Sensor.Type.Laser):
            # if laser, wait for the sensor data to be updated and then print it
            olddata = sensors[ienablesensor].GetSensorData(rave.Sensor.Type.Laser)
            while True:
                data = sensors[ienablesensor].GetSensorData(rave.Sensor.Type.Laser)
                if data.stamp != olddata.stamp:
                    break
                time.sleep(0.1)
            print 'sensor data: ',data.ranges                        
        time.sleep(5)
        ienablesensor = (ienablesensor+1)%len(sensors)

def gazebo_test():
    rospy.init_node('talker',anonymous=True)
    r = pr2.PR2()
    
    IPython.embed()

if __name__ == '__main__':
    gazebo_test()
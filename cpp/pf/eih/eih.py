import IPython

import openravepy as rave
from openravepy import Sensor
import time

def eih():
    env = rave.Environment()
    env.Load('data/testwamcamera.env.xml')
    
    env.SetViewer('qtcoin')
    
    ienablesensor = 0
    while True:
        sensors = env.GetSensors()
        for i,sensor in enumerate(sensors):
            if i==ienablesensor:
                sensor.Configure(Sensor.ConfigureCommand.PowerOn)
                sensor.Configure(Sensor.ConfigureCommand.RenderDataOn)
            else:
                sensor.Configure(Sensor.ConfigureCommand.PowerOff)
                sensor.Configure(Sensor.ConfigureCommand.RenderDataOff)
        print 'showing sensor %s, try moving obstacles'%(sensors[ienablesensor].GetName())
        if sensors[ienablesensor].Supports(Sensor.Type.Laser):
            # if laser, wait for the sensor data to be updated and then print it
            olddata = sensors[ienablesensor].GetSensorData(Sensor.Type.Laser)
            while True:
                data = sensors[ienablesensor].GetSensorData(Sensor.Type.Laser)
                if data.stamp != olddata.stamp:
                    break
                time.sleep(0.1)
            print 'sensor data: ',data.ranges                        
        time.sleep(5)
        ienablesensor = (ienablesensor+1)%len(sensors)

    
    IPython.embed()

def pr2_sensors():
    env = rave.Environment()
    #env.Load('robots/pr2-beta-sim.robot.xml')
    env.Load('envs/testpr2sensors.env.xml')
    r = env.GetRobots()[0]
    
    env.SetViewer('qtcoin') 
    
    ienablesensor = 0
    while True:
        start_time = time.time()
        sensors = env.GetSensors()
        for i,sensor in enumerate(sensors):
            if i==ienablesensor:
                sensor.Configure(Sensor.ConfigureCommand.PowerOn)
                sensor.Configure(Sensor.ConfigureCommand.RenderDataOn)
            else:
                sensor.Configure(Sensor.ConfigureCommand.PowerOff)
                sensor.Configure(Sensor.ConfigureCommand.RenderDataOff)
        print 'showing sensor %s, try moving obstacles'%(sensors[ienablesensor].GetName())
        if sensors[ienablesensor].Supports(Sensor.Type.Laser):
            # if laser, wait for the sensor data to be updated and then print it
            olddata = sensors[ienablesensor].GetSensorData(Sensor.Type.Laser)
            while True:
                data = sensors[ienablesensor].GetSensorData(Sensor.Type.Laser)
                if data.stamp != olddata.stamp:
                    break
                time.sleep(0.1)
            print 'sensor data: ',data.ranges                        
        #time.sleep(5)
        ienablesensor = (ienablesensor+1)%len(sensors)
        print('Elapsed: {0}'.format(time.time() - start_time))

    
    IPython.embed()

def pr2_flashlidar():
    env = rave.Environment()
    env.Load('envs/testpr2sensors.env.xml')
    
    env.SetViewer('qtcoin') 
    time.sleep(1)
    
    start_time = time.time()
    sensors = [s for s in env.GetSensors() if s.GetName().find("flashlidar") != -1]
    lidar = sensors[0]
    
    lidar.Configure(Sensor.ConfigureCommand.PowerOn)
    #lidar.Configure(Sensor.ConfigureCommand.RenderDataOn)
            
    while True:
        start_time = time.time()
        olddata = lidar.GetSensorData(Sensor.Type.Laser)
        while True:
            data = lidar.GetSensorData(Sensor.Type.Laser)
            if data.stamp != olddata.stamp:
                break
            time.sleep(0.1)
        print('Elapsed: {0}'.format(time.time() - start_time))
        break

    
    IPython.embed()


def pr2_links():
    env = rave.Environment()
    #env.Load('robots/pr2-beta-sim.robot.xml')
    env.Load('envs/testpr2sensors.env.xml')
    r = env.GetRobots()[0]
    
    # put sensors at wide_stereo_link
    
    IPython.embed()

if __name__ == '__main__':
    #eih()
    #pr2_sensors()
    pr2_flashlidar()
    #pr2_links()
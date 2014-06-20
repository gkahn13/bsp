import time
import numpy as np

import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx

import utils

import IPython

def env_mesh():
    env = rave.Environment()
    env.Load('envs/pr2-test.env.xml')
    env.SetViewer('qtcoin')
    time.sleep(1)
    robot = env.GetRobots()[0]
    
    handles = list()
    
    cam_pose = tfx.pose(robot.GetSensor('head_cam').GetTransform())
    handles += utils.plot_transform(env, cam_pose.matrix)
    
    mug = env.GetKinBody('mug')
    mug_pose = tfx.pose(mug.GetTransform())
    
    #start, end = tfx.pose(robot.GetLink('wide_stereo_link').GetTransform()).position.array, mug_pose.position.array
    start, end = cam_pose.position.array, mug_pose.position.array
    origin = start
    direction = (end - start) / np.linalg.norm(end-start)
    handles += utils.plot_segment(env, origin, origin+5*direction)
    
    #rays = np.hstack((origin+5*direction, origin))
    rays = np.hstack((origin, origin+5*direction))
    with env:
        t, hits = env.CheckCollisionRays(np.array([rays]), None)
    
    print((t, hits))
    
    #handles.append(utils.plot_point(env, hits1[0,0:3], size=.04, color=(1,0,0)))
    
    
    IPython.embed()

def rosen_test():
    env = rave.Environment()
    
    #env.Load('envs/pr2-test.env.xml')
    #env.Load('data/box3.kinbody.xml')
    body1 = env.ReadKinBodyXMLFile('data/mug1.kinbody.xml')
    env.Add(body1)
    
    env.SetViewer('qtcoin')
    time.sleep(1)
    with env:
        #start = np.array([0,0,0])
        #end = np.array([0,0,1])
        start = np.array([1,0,0])
        end = np.array([0,0,.05])
        
        ray = np.hstack((start, end-start))
        t, hits = env.CheckCollisionRays(np.array([ray]),None)
        print((t, hits))
        
    handles = []
    handles.append(utils.plot_segment(env, start, end, color=(0,1,0)))
    IPython.embed()
    
if __name__ == '__main__':
    #env_mesh()
    rosen_test()
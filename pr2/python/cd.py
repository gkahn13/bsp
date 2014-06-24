import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx

import utils

import IPython

handles = list()

def test_cd():
    global handles
    env = rave.Environment()
    env.Load('../envs/pr2-test.env.xml')
    env.SetViewer('qtcoin')
    
    mug = env.GetKinBody('mug')
    mug_pos = tfx.pose(mug.GetTransform()).position.array 
    mugcd = rave.databases.convexdecomposition.ConvexDecompositionModel(mug)
    if not mugcd.load():
        mugcd.autogenerate()
    
    mugcd_trimesh = mugcd.GenerateTrimeshFromHulls(mugcd.linkgeometry[0][0][1])
    
    new_mug = rave.RaveCreateKinBody(env,'')
    new_mug.SetName('new_mug')
    new_mug.InitFromTrimesh(mugcd_trimesh)
    new_mug.SetTransform(mug.GetTransform())
    #env.Add(new_mug, True)
    env.Remove(mug)
    
    I, V = mugcd_trimesh.indices, mugcd_trimesh.vertices
    for indices in I:
        v0 = mug_pos + V[indices[0],:]
        v1 = mug_pos + V[indices[1],:]
        v2 = mug_pos + V[indices[2],:]
        
        handles += utils.plot_segment(env, v0, v1)
        handles += utils.plot_segment(env, v1, v2)
        handles += utils.plot_segment(env, v2, v0)
    
    IPython.embed()

if __name__ == '__main__':
    test_cd()
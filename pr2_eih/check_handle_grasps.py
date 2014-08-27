"""
subscribe to
-- /kinfu/graspable_points (base_link)
-- /kinfu/table_pose (base_link)
-- /handle_detector/avg_handle_poses

convex decomposition of graspable points
add table
for each handle pose
  call trajopt
  if no collision and reaches handle
    publish pose trajectory to /check_handle_grasps/trajectory
"""

import rospy, roslib
roslib.load_manifest('pr2_utils')
from pr2 import convexify_pointcloud, planner
from pr2_sim import simulator, arm, camera
roslib.load_manifest('tfx')
import tfx

import geometry_msgs.msg as gm
import sensor_msgs.msg as sm

import numpy as np

import IPython

class CheckHandleGrasps:
    def __init__(self):
        self.graspable_points = None
        self.table_pose = None
        self.avg_handle_poses = None
        
        self.graspable_points_sub = rospy.Subscriber('/handle_detector/point_cloud', sm.PointCloud2, self._graspable_points_callback)
        self.table_pose_sub = rospy.Subscriber('/kinfu/table_pose', gm.PoseStamped, self._table_pose_callback) # TODO: change to /kinfu/plane_bounding_box
        self.avg_handle_poses_sub = rospy.Subscriber('/handle_detector/avg_handle_poses',
                                                gm.PoseArray, self._avg_handle_poses_callback)
        
        self.trajectory_pub = rospy.Publisher('/check_handle_grasps/trajectory', gm.PoseArray)
        
        self.sim = simulator.Simulator(view=True)
        self.arm = arm.Arm('right', sim=self.sim)
        self.planner = planner.Planner('right', sim=self.sim)
        
        rospy.sleep(0.5)
        self.sim.update()
        
    ############################
    # callbacks and retreivers #
    ############################
        
    def _graspable_points_callback(self, msg):
        self.graspable_points = msg
        
    def _table_pose_callback(self, msg):
        self.table_pose = tfx.pose(msg)
        
    def _avg_handle_poses_callback(self, msg):
        self.avg_handle_poses = [tfx.pose(p, header=msg.header) for p in msg.poses]
        
    def get_most_recent_callbacks(self):
        """
        Clears all callback variables and blocks
        until they are all updated
        """
        self.graspable_points = None
        self.table_pose = tfx.pose([0,0,0]) # TODO: temp
        self.avg_handle_poses = None
        
        while not rospy.is_shutdown() and \
            (self.graspable_points is None or
             self.table_pose is None or
             self.avg_handle_poses is None):
            self.avg_handle_poses = None # in case publishing handles from old cloud
            if self.graspable_points is None:
                print('Waiting for graspable points')
            if self.table_pose is None:
                print('Waiting for table pose')
            if self.avg_handle_poses is None:
                print('Waiting for avg handle poses')
            print('')
            rospy.sleep(1)
    
    #####################
    # interface methods #
    #####################
    
    def run(self):
        while not rospy.is_shutdown():
            print('Getting most recent callbacks')
            self.get_most_recent_callbacks()
            
            graspable_points = self.graspable_points
            table_pose = self.table_pose
            avg_handle_poses = self.avg_handle_poses
            
            #assert graspable_points.header.frame_id.replace('/','').count(table_pose.frame.replace('/','')) > 0
            
            print('Convexifying point cloud')
            self.sim.update()
            self.sim.clear_kinbodies()
            convexify_pointcloud.add_convexified_pointcloud_to_env(self.sim, graspable_points)
            
            print('Adding table mesh')
            ########
            # TODO #
            ########
            
            self.sim.update()
            arm_pose = self.arm.get_pose()
            avg_handle_poses = sorted(avg_handle_poses, key=lambda p: np.linalg.norm(p.position.array - arm_pose.position.array))
            
            for i, handle_pose in enumerate(avg_handle_poses):
                self.sim.clear_plots()
                self.sim.update()
                self.sim.plot_transform(self.sim.transform_from_to(handle_pose, handle_pose.frame, 'world'))
                print('Checking if handle pose {0} is collision free and reachable'.format(i))
                print('Press enter to continue')
                raw_input()
                
                joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), handle_pose, ignore_orientation=False)
                
                is_collision_free = joint_traj is not None
                if is_collision_free:
                    print('Trajectory is collision free')
                    
                    poses = [tfx.pose(self.arm.fk(joints)) for joints in joint_traj]
                    dist = np.linalg.norm(poses[-1].position.array - handle_pose.position.array)
                    print('Distance of final pose to handle: {0}'.format(dist))
                    is_within_reach = dist < .001
                
                    if is_within_reach:
                        print('Handle pose is within reach at final timestep')
                        trajectory = gm.PoseArray()
                        trajectory.header.stamp = rospy.Time.now()
                        trajectory.header.frame_id = 'base_link'
                        
                        for pose in poses:
                            trajectory.poses.append(pose.msg.Pose())
                            self.sim.plot_transform(self.sim.transform_from_to(pose, pose.frame, 'world'))
                        self.arm.set_pose(poses[-1])
                        
                        self.trajectory_pub.publish(trajectory)
                        print('Press enter to continue')
                        raw_input()
                    else:
                        print('Handle pose is not within reach at final timestep')
                else:
                    print('Trajectory is not collision free')
                print('')
                    
            rospy.sleep(0.1)
    
if __name__ == '__main__':
    rospy.init_node('check_handle_grasps', anonymous=True)
    
    chg = CheckHandleGrasps()
    chg.run()
    
    
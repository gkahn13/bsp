#!/usr/bin/env python

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

class CheckHandleGrasps:
    def __init(self):
        self.graspable_points = None
        self.table_pose = None
        self.avg_handle_poses = None
        
        self.graspable_points_sub = rospy.Subscriber('/kinfu/graspable_points', sm.PointCloud2, self._graspable_points_callback)
        self.table_pose_sub = rospy.Subscriber('/kinfu/table_pose', gm.PoseStamped, self._table_pose_callback)
        self.avg_handle_poses_sub = rospy.Subscriber('/handle_detector/avg_handle_poses',
                                                gm.PoseArray, self._avg_handle_poses_callback)
        
        self.trajectory_pub = rospy.Publisher('/check_handle_grasps/trajectory', gm.PoseArray)
        
        self.sim = simulator.Simulator(view=True)
        self.arm = arm.Arm('right', sim=self.sim)
        self.planner = planner.Planner('right', sim=self.sim)
        
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
        self.table_pose = None
        self.avg_handle_poses = None
        
        while not rospy.is_shutdown() and \
            (self.graspable_points is None or
             self.table_pose is None or
             self.avg_handle_poses is None):
            self.avg_handle_poses = None # in case publishing handles from old cloud
            rospy.sleep(0.1)
    
    #####################
    # interface methods #
    #####################
    
    def run(self):
        while not rospy.is_shutdown():
            self.get_most_recent_callbacks()
            
            graspable_points = self.graspable_points
            table_pose = self.table_pose
            avg_handle_poses = self.avg_handle_poses
            
            self.sim.clear_kinbodies()
            convexify_pointcloud.add_convexified_pointcloud_to_env(self.sim, graspable_points)
            
            for i, handle_pose in enumerate(avg_handle_poses):
                print('Checking if handle pose {0} is collision free and reachable'.format(i))
                self.sim.clear_plots()
                self.sim.plot_transform(self.sim.transform_from_to(handle_pose, handle_pose.frame, 'world'))
                
                self.sim.update()
                joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), handle_pose, ignore_orientation=True)
                poses = [tfx.pose(self.arm.fk(joints)) for joints in joint_traj]
                
                is_within_reach = (poses[-1].position - handle_pose).norm < .001
                is_collision_free = True # TODO
                
                if is_within_reach and is_collision_free:
                    trajectory = gm.PoseArray()
                    trajectory.header.stamp = rospy.Time.now()
                    trajectory.header.frame_id = 'base_link'
                    
                    for pose in poses:
                        trajectory.poses.append(pose.msg.Pose())
                        
                    self.trajectory_pub.publish(trajectory)
                    
            rospy.sleep(0.1)
    
if __name__ == '__main__':
    rospy.init_node('check_handle_grasps', anonymous=True)
    
    chg = CheckHandleGrasps()
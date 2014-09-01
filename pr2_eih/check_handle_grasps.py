"""
subscribe to
-- /kinfu/graspable_points (base_link)
-- /kinfu/plane_bounding_box (base_link)
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
roslib.load_manifest('pcl_utils')
import pcl_utils.msg
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
        self.table_extents = None
        self.avg_handle_poses = None
        
        self.graspable_points_sub = rospy.Subscriber('/handle_detector/point_cloud', sm.PointCloud2, self._graspable_points_callback)
        self.table_sub = rospy.Subscriber('/kinfu/plane_bounding_box', pcl_utils.msg.BoundingBox, self._table_callback)
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
        
    def _table_callback(self, msg):
        assert(msg.header.frame_id.count('base_link') > 0)
        assert(len(msg.vectors) == 3)
        
        # calculate rotation matrix and extent lengths
        extents = [tfx.point(e) for e in msg.vectors]
        R = np.zeros((3,3))
        for i, e in enumerate(extents):
            R[:,i] = e.array/e.norm
            
        extent_lengths = np.array([e.norm for e in extents])
            
        self.table_pose = tfx.pose(msg.center, R).matrix
        self.table_extents = extent_lengths
        
    def _avg_handle_poses_callback(self, msg):
        self.avg_handle_poses = [tfx.pose(p, header=msg.header) for p in msg.poses]
        
    def get_most_recent_callbacks(self):
        """
        Clears all callback variables and blocks
        until they are all updated
        """
#         self.table_pose = None 
#         self.table_extents = None 
#         print('Waiting for table pose and extents')
#         while not rospy.is_shutdown() and (self.table_pose is None or self.table_extents is None):
#             rospy.sleep(.01)
#         
#         self.graspable_points = None            
#         print('Waiting for graspable points')
#         while not rospy.is_shutdown() and self.graspable_points is None:
#             rospy.sleep(.01)
#         
#         self.avg_handle_poses = None
#         print('Waiting for avg handle poses')
#         while not rospy.is_shutdown() and self.avg_handle_poses is None:
#             rospy.sleep(.01)

        self.table_pose = None
        self.table_extents = None
        printed_table = False
        self.graspable_points = None
        printed_graspable_points = False
        self.avg_handle_poses = None
        printed_avg_handle_poses = False
        while not rospy.is_shutdown() and \
            (self.table_pose is None or self.graspable_points is None or self.avg_handle_poses is None):
            if self.table_pose is not None and not printed_table:
                print('Found table pose')
                printed_table = True
            if self.graspable_points is not None and not printed_graspable_points:
                print('Found graspable points')
                printed_graspable_points = True
            if self.avg_handle_poses is not None and not printed_avg_handle_poses:
                print('Found avg handle poses')
                printed_avg_handle_poses = True
            rospy.sleep(0.5)
    
    #####################
    # interface methods #
    #####################
    
    def run(self):
        grasp_frame = 'r_gripper_center_frame'
        
        while not rospy.is_shutdown():
            print('Getting most recent callbacks')
            rospy.sleep(0.05)
            self.sim.update()
            self.get_most_recent_callbacks()
            
            graspable_points = self.graspable_points
            table_pose = self.table_pose
            table_extents = self.table_extents
            avg_handle_poses = self.avg_handle_poses
            
             # fails if too few points
            if graspable_points.width < 50: 
                print('Not enough points: {0}'.format(graspable_points.width))
                continue
            if len(avg_handle_poses) == 0:
                print('No handles found')
                continue
            
            #assert graspable_points.header.frame_id.replace('/','').count(table_pose.frame.replace('/','')) > 0
            
            print('Convexifying point cloud')
            self.sim.update()
            self.sim.clear_kinbodies()
            try:
                convexify_pointcloud.add_convexified_pointcloud_to_env(self.sim, graspable_points)
            except Exception as e:
                print('Error convexifying point cloud: {0}'.format(e))
                continue
            
            
            print('Adding table mesh')
            self.sim.add_box(table_pose, table_extents, check_collision=False)
            
            arm_pose = self.arm.get_pose()
            #avg_handle_poses = sorted(avg_handle_poses, key=lambda p: np.linalg.norm(p.position.array - arm_pose.position.array))
            avg_handle_poses = sorted(avg_handle_poses, key=lambda p: -p.position.z)
            
            for i, handle_pose in enumerate(avg_handle_poses):
                self.sim.clear_plots()
                self.sim.update()
                self.sim.plot_transform(self.sim.transform_from_to(handle_pose, handle_pose.frame, 'world'))
                print('Checking if handle pose {0} is collision free and reachable'.format(i))
                
                joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), handle_pose,
                                                               ignore_orientation=True, link_name=grasp_frame)
                
                is_collision_free = joint_traj is not None
                if is_collision_free:
                    print('Trajectory is collision free')
                    
                    poses = [tfx.pose(self.arm.fk(joints)) for joints in joint_traj]
                    
                    self.arm.set_joints(joint_traj[-1])
                    T = self.sim.transform_from_to(np.eye(4), grasp_frame, 'base_link')
                    final_grasp_pose = tfx.pose(T, frame='base_link')
                    
                    dist = np.linalg.norm(final_grasp_pose.position.array - handle_pose.position.array)
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
#                         print('Press enter to continue')
#                         raw_input()
                        break # only publish closest successful grasp
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
    
    
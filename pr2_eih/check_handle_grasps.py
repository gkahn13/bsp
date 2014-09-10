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
import std_msgs.msg
import trajectory_msgs.msg as tm

import numpy as np

# import warnings
# warnings.filterwarnings('error')

import IPython

class CheckHandleGrasps:
    def __init__(self):
        self.graspable_points = None
        self.table_pose = None
        self.table_extents = None
        self.avg_handle_poses = None
        self.is_reset_kinfu = False
        self.home_pose = None
        
        min_positions = np.array([0.25,-1,0.3])
        max_positions = np.array([1,1,1.5])
        self.point_cloud_filter = lambda point: min(min_positions < point) and min(point < max_positions)
#         self.point_cloud_filter = lambda point: True
        
        self.graspable_points_sub = rospy.Subscriber('/kinfu/graspable_points', sm.PointCloud2, self._graspable_points_callback)
        self.table_sub = rospy.Subscriber('/kinfu/plane_bounding_box', pcl_utils.msg.BoundingBox, self._table_callback)
        self.avg_handle_poses_sub = rospy.Subscriber('/handle_detector/avg_handle_poses',
                                                gm.PoseArray, self._avg_handle_poses_callback)
        self.reset_kinfu_sub = rospy.Subscriber('/reset_kinfu', std_msgs.msg.Empty, self._reset_kinfu_callback)
        self.home_pose_sub = rospy.Subscriber('/bsp/home_pose', gm.PoseStamped, self._home_pose_callback)
        
        self.grasp_joint_traj_pub = rospy.Publisher('/check_handle_grasp/grasp_joint_trajectory', tm.JointTrajectory)
        self.return_grasp_joint_traj_pub = rospy.Publisher('/check_handle_grasp/return_grasp_joint_trajectory', tm.JointTrajectory)
        
        self.sim = simulator.Simulator(view=True)
        self.arm = arm.Arm('right', sim=self.sim)
        self.planner = planner.Planner('right', sim=self.sim)
        
        self.arm.set_posture('mantis')
        self.home_pose = self.arm.get_pose()
        
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
        extent_lengths[-1] *= 1
            
        self.table_pose = tfx.pose(msg.center, R).matrix
        self.table_extents = extent_lengths
        
    def _avg_handle_poses_callback(self, msg):
        self.avg_handle_poses = [tfx.pose(p, header=msg.header) for p in msg.poses]
    
    def _reset_kinfu_callback(self, msg):
        self.is_reset_kinfu = True
        
    def _home_pose_callback(self, msg):
        print('Home pose received')
        self.home_pose = tfx.pose(msg)
    
    def get_most_recent_callbacks(self):
        """
        Clears all callback variables and blocks
        until they are all updated
        """
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
        while not rospy.is_shutdown():
            print('Getting most recent callbacks')
            rospy.sleep(0.05)
            self.sim.update()
            self.get_most_recent_callbacks()
            self.is_reset_kinfu = False
            
            graspable_points = self.graspable_points
            table_pose = self.table_pose
            table_extents = self.table_extents
            avg_handle_poses = self.avg_handle_poses
            
             # fails if too few points
            print('Point cloud size: {0}'.format(graspable_points.width))
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
                convexify_pointcloud.add_convexified_pointcloud_to_env(self.sim, graspable_points,
                                                                       point_cloud_filter=self.point_cloud_filter)
            except Warning as e:
                print('Error convexifying point cloud: {0}'.format(e))
                continue
            except Exception as e:
                print('Error convexifying point cloud: {0}'.format(e))
                continue
            
            
            print('Adding table mesh')
            self.add_table(table_pose, table_extents)
            
            arm_pose = self.arm.get_pose()
            #avg_handle_poses = sorted(avg_handle_poses, key=lambda p: np.linalg.norm(p.position.array - arm_pose.position.array))
            avg_handle_poses = sorted(avg_handle_poses, key=lambda p: -p.position.z)
            
            for i, handle_pose in enumerate(avg_handle_poses):
                if self.is_reset_kinfu:
                    break
                
                print('Checking if handle pose {0} is collision free and reachable'.format(i))
                grasp_joint_traj = self.get_grasp_trajectory(handle_pose)
                
                if grasp_joint_traj is not None:
                    grasp_joints = grasp_joint_traj[-1]
                    return_grasp_joint_traj = self.get_return_from_grasp_trajectory(grasp_joints, table_pose)
                    
                    if return_grasp_joint_traj is None:
                        return_grasp_joint_traj = list()
                    
                    self.plot_joint_traj(grasp_joint_traj)
                    self.plot_joint_traj(return_grasp_joint_traj)
                    
                    self.publish_joint_traj(self.grasp_joint_traj_pub, grasp_joint_traj)
                    self.publish_joint_traj(self.return_grasp_joint_traj_pub, return_grasp_joint_traj)
                    break
                                    
            rospy.sleep(0.1)
            
    ##################
    # Helper methods #
    ##################
    
    def add_table(self, table_pose, table_extents):
        """
        Adds table to environment. If colliding with robot, scoots in +x direction until not colliding
        :param table_pose np.ndarray 4x4 in base_link
        :param table_extents np.ndarray 3d in base_link
        """
        while not self.sim.add_box(table_pose, table_extents, check_collision=True):
            table_pose[0,3] += 0.05 # move forward until not in collision
            
    def get_grasp_trajectory(self, handle_pose, grasp_frame='r_gripper_center_frame'):
        """
        :param handle_pose tfx.pose desired grasp pose
        :param grasp_frame what trajopt plans for
        :return None if failure, else list of joints
        """
        self.sim.clear_plots()
        self.sim.update()
        self.sim.plot_transform(self.sim.transform_from_to(handle_pose, handle_pose.frame, 'world'))
        
        joint_traj = self.planner.get_grasp_joint_trajectory(self.arm.get_joints(), handle_pose,
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
                return joint_traj
            else:
                print('Handle pose is not within reach at final timestep')
        else:
            print('Trajectory is not collision free')
        print('')
        
        return None
    
    def get_return_from_grasp_trajectory(self, grasp_joints, table_pose):
        """
        Get return grasp trajectory by going up, then planning a high trajectory
        :param grasp_joints joints of the grasp pose
        :param table_pose
        :return None if failure, else list of joints
        """
        print('Getting trajectory from grasp back to home')
        table_pose = tfx.pose(table_pose)
        return_grasp_traj = list()
        
        print('Getting above grasp pose')
        self.arm.set_joints(grasp_joints)
        grasp_pose = self.arm.get_pose()
        
        above_grasp_pose = grasp_pose
        above_grasp_pose.position.z += grasp_pose.position.z - table_pose.position.z
        above_grasp_joints = self.arm.ik(above_grasp_pose)
        
        if above_grasp_joints is None:
            above_grasp_joints = grasp_joints
            
        self.arm.set_joints(above_grasp_joints)
        self.sim.remove_colliding_kinbodies() # so we start out feasible
        
        if self.home_pose is None:
            return None
        
        joint_traj = self.planner.get_return_from_grasp_joint_trajectory(above_grasp_joints, self.home_pose)
        if joint_traj is None or len(joint_traj) == 0:
            return None
        return joint_traj
            
    def publish_joint_traj(self, pub, joint_traj):
        joint_traj_msg = tm.JointTrajectory()
        joint_traj_msg.header.stamp = rospy.Time.now()
        joint_traj_msg.points = [tm.JointTrajectoryPoint(positions=joints) for joints in joint_traj]
        
        pub.publish(joint_traj_msg)
        
    def plot_joint_traj(self, joint_traj):
        poses = [tfx.pose(self.arm.fk(joints)) for joints in joint_traj]
        for pose in poses:
            self.sim.plot_transform(self.sim.transform_from_to(pose.matrix, 'base_link', 'world'))
        
if __name__ == '__main__':
    rospy.init_node('check_handle_grasps', anonymous=True)
    
    chg = CheckHandleGrasps()
    chg.run()
    
    
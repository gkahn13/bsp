import rospy, roslib, rosbag
import tf
roslib.load_manifest('sensor_msgs')
roslib.load_manifest('handle_detector')
roslib.load_manifest('tfx')
import sensor_msgs.msg as sm
import geometry_msgs.msg as gm
import handle_detector.msg as hd_msg
import tfx

import sys
import argparse
import threading
import numpy as np

import utils

import IPython

# don't forget
# roslaunch localization_sensor.launch

class Playback:
    def __init__(self, bag_name, pc_topic, pc_publish_topic='/camera/depth_registered/points'):
        self.pc_topic = pc_topic
        
        self.num_pcs = 0
        self.bag = rosbag.Bag(bag_name + '.bag')
        self.bag_msgs = list()
        for topic, msg, t in self.bag.read_messages():
            self.bag_msgs.append((topic, msg, t))
            if topic.count(self.pc_topic) == 1:
                self.num_pcs += 1
                
        print('Number of messages in {0}.bag: {1}'.format(bag_name, len(self.bag_msgs)))
        print('Number of point-clouds read: {0}'.format(self.num_pcs))
        print('Publishing to point-cloud topic: {0}'.format(pc_publish_topic))

        self.handles_sub = rospy.Subscriber('/handle_detector/handle_list', hd_msg.HandleListMsg, self._handles_callback)
        self.handles_pose_pub = rospy.Publisher('/handle_detector/handle_poses', gm.PoseArray)
        self.avg_handles_pose_pub = rospy.Publisher('/handle_detector/avg_handle_poses', gm.PoseArray)
        self.handle_pose_0_pub = rospy.Publisher('/handle_detector/handle_pose_0', gm.PoseStamped)
        
        self.pc_pub = rospy.Publisher(pc_publish_topic, sm.PointCloud2)
        self.br = tf.TransformBroadcaster()
        
        self.exited = False
        
        self.tfs = dict()
        self.tfs_thread = threading.Thread(target=self._tfs_loop)
        self.tfs_thread.start()
        
        self.handle_list_msg = hd_msg.HandleListMsg()
        self.handles_thread = threading.Thread(target=self._handles_loop)
        self.handles_thread.start()
        
    def _tfs_loop(self):
        """
        For each transform in self.tfs,
        continuously publishes the last one
        """
        try:
            while not rospy.is_shutdown() and not self.exited:
                tfs = dict(self.tfs)
                for t_list in tfs.values():
                    t = t_list[-1]
                    trans = [t.transform.translation.x, t.transform.translation.y, t.transform.translation.z]
                    rot = [t.transform.rotation.x, t.transform.rotation.y, t.transform.rotation.z, t.transform.rotation.w]
                    self.br.sendTransform(trans,
                                          rot,
                                          rospy.Time.now(),
                                          t.child_frame_id,
                                          t.header.frame_id)
                    
                rospy.sleep(.001)
        except (KeyboardInterrupt, SystemExit):
            sys.exit()
            
    def _handles_loop(self):
        """
        For each handle in HandleListMsg,
        calculate average pose
        """
        rospy.sleep(5)
        
        while not rospy.is_shutdown() and not self.exited:
            rospy.sleep(.01)
            
            handle_list_msg = self.handle_list_msg
            
            pose_array = gm.PoseArray()
            pose_array.header.frame_id = handle_list_msg.header.frame_id
            pose_array.header.stamp = rospy.Time.now()
            
            avg_pose_array = gm.PoseArray()
            avg_pose_array.header.frame_id = handle_list_msg.header.frame_id
            avg_pose_array.header.stamp = rospy.Time.now()
    
            cam_to_base = tfx.lookupTransform('base_link', handle_list_msg.header.frame_id).matrix[:3,:3]
            switch = np.matrix([[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]])        
            for handle in handle_list_msg.handles:
                all_poses = [tfx.pose(cylinder.pose, stamp=rospy.Time.now(), frame=handle_list_msg.header.frame_id) for cylinder in handle.cylinders]
                
                rotated_poses = [tfx.pose(p.position, tfx.tb_angles(p.orientation.matrix*switch)) for p in all_poses]
                filtered_poses = list()
                for rot_pose in rotated_poses:
                    r_base = cam_to_base*rot_pose.orientation.matrix
                    if r_base[0,0] > 0:
                        if r_base[2,2] > 0:
                            rot_pose.orientation = tfx.tb_angles(rot_pose.orientation.matrix*tfx.tb_angles(0,0,180).matrix) 
                        filtered_poses.append(rot_pose)
                
                pose_array.poses += [pose.msg.Pose() for pose in filtered_poses]
                
                if len(filtered_poses) > 0:
                    avg_position = sum([p.position.array for p in filtered_poses])/float(len(filtered_poses))
                    avg_quat = sum([p.orientation.quaternion for p in filtered_poses])/float(len(filtered_poses))
                    avg_pose_array.poses.append(tfx.pose(avg_position, avg_quat).msg.Pose())
                
                
            if len(pose_array.poses) > 0:
                self.handles_pose_pub.publish(pose_array)
                self.avg_handles_pose_pub.publish(avg_pose_array)
                
            
    def _handles_callback(self, msg):
        self.handle_list_msg = msg
                        
    def play(self):
        """
        Switches between point-clouds in the bag file
        while continuously publishing the transforms
        """
        curr_index, topic = 0, ''
        while topic.count(self.pc_topic) == 0:
            curr_index += 1
            topic, msg, t = self.bag_msgs[curr_index]
            
            if topic.count('tf') == 1:
                for t in msg.transforms:
                    key = (t.child_frame_id, t.header.frame_id)
                    if key in self.tfs.keys():
                        self.tfs[key].append(t)
                    else:
                        self.tfs[key] = [t]
            elif topic.count(self.pc_topic) == 1:
                msg.header.stamp = rospy.Time.now()
                self.pc_pub.publish(msg)
        
        curr_pc = 0
        input = ''
        while not rospy.is_shutdown():
            print('At point-cloud {0} (of range [{1},{2}])'.format(curr_pc, 0, self.num_pcs-1))
            print('Enter -/+ to proceed')
            input = utils.Getch.getch()
            print('')
            
            self.tfs = dict()
            
            if input == '-' or input == '+':
                incr = -1 if input == '-' else 1
                
                if (input == '-' and curr_pc <= 0) or (input == '+' and curr_pc >= self.num_pcs -1):
                    print('Point-cloud {0} out of range'.format(curr_pc + incr))
                    continue
                else:
                    curr_pc += incr
                    
                topic = ''
                while topic.count(self.pc_topic) == 0:
                    curr_index += incr
                    topic, msg, t = self.bag_msgs[curr_index]
                    
                    if topic.count('tf') == 1:
                        for t in msg.transforms:
                            key = (t.child_frame_id, t.header.frame_id)
                            if key in self.tfs.keys():
                                self.tfs[key].append(t)
                            else:
                                self.tfs[key] = [t]
                    elif topic.count(self.pc_topic) == 1:
                        msg.header.stamp = rospy.Time.now()
                        self.pc_pub.publish(msg)
                        
                    rospy.sleep(.001)
            elif input == 'q':
                print('Exiting')
                self.exited = True
                break
            
        self.tfs_thread.join()
        self.handles_thread.join()
  
if __name__ == '__main__':
    rospy.init_node('playback', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--bag',type=str,default='bags/spray_can')
    parser.add_argument('--topic',type=str,default='save_pc')
    args = parser.parse_args()
    
    playback = Playback(args.bag, args.topic)
    playback.play()
    sys.exit()

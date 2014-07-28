import rospy, roslib, rosbag
import tf
roslib.load_manifest('sensor_msgs')
roslib.load_manifest('handle_detector')
import sensor_msgs.msg as sm
import geometry_msgs.msg as gm
import handle_detector.msg as hd_msg

import sys
import argparse
import threading

import utils

import IPython

# don't forget
# roslaunch localization_sensor.launch

class Playback:
    def __init__(self, bag_name, pc_publish_topic='/camera/depth_registered/points'):
        self.num_pcs = 0
        self.bag = rosbag.Bag(bag_name + '.bag')
        self.bag_msgs = list()
        for topic, msg, t in self.bag.read_messages():
            self.bag_msgs.append((topic, msg, t))
            if topic == 'save_pc':
                self.num_pcs += 1
                
        print('Number of point-clouds read: {0}'.format(self.num_pcs))

        self.handles_sub = rospy.Subscriber('/localization/handle_list', hd_msg.HandleListMsg, self._handles_callback)
        self.handles_pose_pub = rospy.Publisher('/localization/handle_poses', gm.PoseArray)
        
        self.pc_pub = rospy.Publisher(pc_publish_topic, sm.PointCloud2)
        self.br = tf.TransformBroadcaster()
        
        self.tfs = dict()
        self.exited = False
        self.tfs_thread = threading.Thread(target=self._tfs_loop)
        self.tfs_thread.start()
        
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
            
    def _handles_callback(self, msg):
        pose_array = gm.PoseArray()
        
        pose_array.header = msg.header
        for handle in msg.handles:
            for cylinder in handle.cylinders:
                pose_array.poses.append(cylinder.pose)
                
        if len(pose_array.poses) > 0:
            self.handles_pose_pub.publish(pose_array)
            
    def play(self):
        """
        Switches between point-clouds in the bag file
        while continuously publishing the transforms
        """
        curr_index, topic = 0, ''
        while topic != 'save_pc':
            curr_index += 1
            topic, msg, t = self.bag_msgs[curr_index]
            
            if topic == 'tf':
                for t in msg.transforms:
                    key = (t.child_frame_id, t.header.frame_id)
                    if key in self.tfs.keys():
                        self.tfs[key].append(t)
                    else:
                        self.tfs[key] = [t]
            elif topic == 'save_pc':
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
                while topic != 'save_pc':
                    curr_index += incr
                    topic, msg, t = self.bag_msgs[curr_index]
                    
                    if topic == 'tf':
                        for t in msg.transforms:
                            key = (t.child_frame_id, t.header.frame_id)
                            if key in self.tfs.keys():
                                self.tfs[key].append(t)
                            else:
                                self.tfs[key] = [t]
                    elif topic == 'save_pc':
                        msg.header.stamp = rospy.Time.now()
                        self.pc_pub.publish(msg)
                        
                    rospy.sleep(.001)
            elif input == 'q':
                print('Exiting')
                self.exited = True
                break
  
if __name__ == '__main__':
    rospy.init_node('playback', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--bag',type=str,default='bags/spray_can')
    args = parser.parse_args()
    
    playback = Playback(args.bag)
    playback.play()
    sys.exit()

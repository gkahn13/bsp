import rospy, roslib, rosbag
import tf
roslib.load_manifest('sensor_msgs')
import sensor_msgs.msg as sm

import argparse

import IPython

# don't forget
# roslaunch localization_sensor.launch

class Playback:
    def __init__(self, bag_name, pc_publish_topic='/camera/depth_registered/points'):
        self.pcs = list()
        self.tfs = list()
        
        self.bag = rosbag.Bag(bag_name)
        for topic, msg, t in self.bag.read_messages():
            if topic == 'tf':
                self.tfs.append(msg)
            if topic == 'save_pc':
                self.pcs.append(msg)
        
        rospy.loginfo('Number of point-clouds read: {0}'.format(len(self.pcs)))
        rospy.loginfo('Number of tfs saved: {0}'.format(len(self.tfs)))
        
        self.pub = rospy.Publisher(pc_publish_topic, sm.PointCloud2)
        self.br = tf.TransformBroadcaster()
    
    def run_all(self):
        num_pcs_published = 0
        tfs = dict()
        
        for topic, msg, t in self.bag.read_messages():
            
            if topic == 'tf':
                for t in msg.transforms:
                    key = (t.child_frame_id, t.header.frame_id)
                    if key in tfs.keys():
                        tfs[key].append(t)
                    else:
                        tfs[key] = [t]
                    # translation, rotation, time, child, parent
                    trans = [t.transform.translation.x, t.transform.translation.y, t.transform.translation.z]
                    rot = [t.transform.rotation.w, t.transform.rotation.x, t.transform.rotation.y, t.transform.rotation.z]
                    self.br.sendTransform(trans,
                                          rot,
                                          t.header.stamp,
                                          t.child_frame_id,
                                          t.header.frame_id)
            if topic == 'save_pc':
                rospy.loginfo('publishing loop')
                while not rospy.is_shutdown(): # put this loop on a separate thread
                    msg.header.stamp = rospy.Time.now()
                    self.pub.publish(msg)
                    
                    #self.br.sendTransform([-0.106925711716, -0.00652973239027, -0.0566985277547],
                    #              [0.501389436235, 0.466364575859, 0.509987956093, 0.520600614921],
                    #              rospy.Time.now(),
                    #              'camera_rgb_optical_frame',
                    #              'r_gripper_tool_frame')
                    
                    for t_list in tfs.values():
                        t = t_list[-1]
                        trans = [t.transform.translation.x, t.transform.translation.y, t.transform.translation.z]
                        rot = [t.transform.rotation.w, t.transform.rotation.x, t.transform.rotation.y, t.transform.rotation.z]
                        self.br.sendTransform(trans,
                                              rot,
                                              rospy.Time.now(),
                                              t.child_frame_id,
                                              t.header.frame_id)
                        
                    rospy.sleep(.01)
                    
                
                if num_pcs_published > 0:
                    print('Press enter')
                    raw_input()
                
                rospy.loginfo('Published pc {0}'.format(num_pcs_published))
                num_pcs_published += 1
            
            if rospy.is_shutdown():
                break
            
            rospy.sleep(.001)
        #for i, pc in enumerate(self.pcs):
        #    rospy.loginfo('Press enter to publish pc {0}'.format(i))
        #    raw_input()
        #    
        #    self.pub.publish(pc)
    
if __name__ == '__main__':
    rospy.init_node('playback', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--bag',type=str,default='bags/spray_can.bag')
    args = parser.parse_args()
    
    playback = Playback(args.bag)
    #playback.run_all()
    
    IPython.embed()
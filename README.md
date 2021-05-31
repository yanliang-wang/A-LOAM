## 0. INTRODUCTION

**This repository is forked from A-LOAM, and I added some Chinese comments and the demo of RS-Ruby Lite(80-beam) in the code.** 

LOAM is a very simple but wonderful SLAM algorithm for line scan LiDAR. It only requires the x,y,z channel of the point cloud and it can support any pose of LiDAR theoretically without any modification. I think the only modification to apply LOAM to your LiDAR is to add the info(vertical resolution, number of scans, etc.) of your LiDAR to the code, which guarantees LOAM can know the scan id of every point. **The demo of RS-Ruby Lite** introduces how to add the info of Lidar to the code.

##  How to add the info of LiDAR

First, print the angle of every point.(comment 214~229 and  uncomment line 231~240 in `scanRegistration.cpp` ).We can get the angle of every point.

```
current angle:  -24.5978        -27.969
current angle:  -23.6955        0.902229
current angle:  -22.7933        0.902225
current angle:  -21.8911        0.902229
current angle:  -20.9889        0.902227
current angle:  -20.0866        0.902225
current angle:  -19.1844        0.902227
current angle:  -18.2822        0.902227
current angle:  -17.38  		0.902229
current angle:  -16.4777        0.902227
...
```

Next, design the computation of `scanId` according to the angles. (line 231~240)

That's all.

The result of RS-Ruby Lite is as follows.

![](./picture/rs-ruby_demo.png)



## The original README.md of A-LOAM is as follows.

# A-LOAM

## Advanced implementation of LOAM

A-LOAM is an Advanced implementation of LOAM (J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time), which uses Eigen and Ceres Solver to simplify code structure. This code is modified from LOAM and [LOAM_NOTED](https://github.com/cuitaixiang/LOAM_NOTED). This code is clean and simple without complicated mathematical derivation and redundant operations. It is a good learning material for SLAM beginners.

<img src="./picture/kitti.png" width = 55% height = 55%/>

**Modifier:** [Tong Qin](http://www.qintonguav.com), [Shaozu Cao](https://github.com/shaozu)


## 1. Prerequisites
### 1.1 **Ubuntu** and **ROS**
Ubuntu 64-bit 16.04 or 18.04.
ROS Kinetic or Melodic. [ROS Installation](http://wiki.ros.org/ROS/Installation)


### 1.2. **Ceres Solver**
Follow [Ceres Installation](http://ceres-solver.org/installation.html).

### 1.3. **PCL**
Follow [PCL Installation](http://www.pointclouds.org/downloads/linux.html).


## 2. Build A-LOAM
Clone the repository and catkin_make:

```
    cd ~/catkin_ws/src
    git clone https://github.com/HKUST-Aerial-Robotics/A-LOAM.git
    cd ../
    catkin_make
    source ~/catkin_ws/devel/setup.bash
```

## 3. Velodyne VLP-16 Example
Download [NSH indoor outdoor](https://drive.google.com/file/d/1s05tBQOLNEDDurlg48KiUWxCp-YqYyGH/view) to YOUR_DATASET_FOLDER. 

```
    roslaunch aloam_velodyne aloam_velodyne_VLP_16.launch
    rosbag play YOUR_DATASET_FOLDER/nsh_indoor_outdoor.bag
```


## 4. KITTI Example (Velodyne HDL-64)
Download [KITTI Odometry dataset](http://www.cvlibs.net/datasets/kitti/eval_odometry.php) to YOUR_DATASET_FOLDER and set the `dataset_folder` and `sequence_number` parameters in `kitti_helper.launch` file. Note you also convert KITTI dataset to bag file for easy use by setting proper parameters in `kitti_helper.launch`. 

```
    roslaunch aloam_velodyne aloam_velodyne_HDL_64.launch
    roslaunch aloam_velodyne kitti_helper.launch
```
<img src="./picture/kitti_gif.gif" width = 720 height = 351 />

## 5. Docker Support
To further facilitate the building process, we add docker in our code. Docker environment is like a sandbox, thus makes our code environment-independent. To run with docker, first make sure [ros](http://wiki.ros.org/ROS/Installation) and [docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/) are installed on your machine. Then add your account to `docker` group by `sudo usermod -aG docker $YOUR_USER_NAME`. **Relaunch the terminal or logout and re-login if you get `Permission denied` error**, type:
```
cd ~/catkin_ws/src/A-LOAM/docker
make build
```
The build process may take a while depends on your machine. After that, run `./run.sh 16` or `./run.sh 64` to launch A-LOAM, then you should be able to see the result.


## 6.Acknowledgements
Thanks for LOAM(J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time) and [LOAM_NOTED](https://github.com/cuitaixiang/LOAM_NOTED).


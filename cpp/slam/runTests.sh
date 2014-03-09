landmarkDir=landmarkTextFiles

for num_landmarks in {3,4,5,6};
do
    cp ${landmarkDir}/landmarks-${num_landmarks}.txt landmarks.txt
    cd ..
    for slam_type in {"traj-plan","belief","state","control"};
    do
	./bin/release/slam-${slam_type}-${num_landmarks}
    done
    cd slam
done
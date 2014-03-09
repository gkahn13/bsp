landmarkDir=landmarkTextFiles

# {10,15,20,25,30,35,40,45,50}
for num_landmarks in {30,35,40,45,50};
do
    cp ${landmarkDir}/landmarks-${num_landmarks}.txt landmarks.txt
    cd ..
    for slam_type in {"state",};
    do
	./bin/release/slam-${slam_type}-${num_landmarks}
    done
    cd slam
done
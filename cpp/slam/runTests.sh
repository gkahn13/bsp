landmarkDir=landmarkTextFiles

ulimit -s 32768

for num_landmarks in {3,4,5};
do
    cp ${landmarkDir}/landmarks-${num_landmarks}.txt landmarks.txt
    cd ..
    for slam_type in {"belief",};
    do
	./bin/release/slam-${slam_type}-${num_landmarks}
    done
    cd slam
done
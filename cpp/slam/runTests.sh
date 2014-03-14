landmarkDir=landmarkTextFiles

ulimit -s 32768

for num_landmarks in {3,4,5,6,10,15,20,25,30,35,40,45,50};
do
    cp ${landmarkDir}/landmarks-${num_landmarks}.txt landmarks.txt
    cd ..
    for slam_type in {"ilqg",};
    do
	./bin/release/slam-${slam_type}-${num_landmarks}
    done
    cd slam
done
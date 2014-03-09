num_runs_per_test=10
landmarkDir=landmarkTextFiles
mkdir -p $landmarkDir

for num_landmarks in {10,15,20,25,30,35,40,45,50};
do
    python generate_landmarks.py $num_landmarks $num_runs_per_test
    cp landmarks.txt ${landmarkDir}/landmarks-${num_landmarks}.txt
    cd ..
    make clean
    for slam_type in {"traj-plan","state","control"};
    do
	make slam-${slam_type} BUILD=release NUM_LANDMARKS=${num_landmarks}
    done
    cd slam
done
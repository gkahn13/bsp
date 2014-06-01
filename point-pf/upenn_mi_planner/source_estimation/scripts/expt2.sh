
rosrun rvo_move move_client __name:=bar  _x:=1.0 _y:=4.6   move_server:=/scarab26/move_server &
rosrun rvo_move move_client __name:=quux _x:=-0.75 _y:=16.5   move_server:=/scarab32/move_server &
rosrun rvo_move move_client __name:=foo  _x:=-6.0 _y:=28.0   move_server:=/scarab29/move_server &


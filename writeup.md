# Model predictive control

## State
The internal state of the MPC is composed of an arbritrary set of state variables that model the system to be controlled. Chosing the right state variable for the control task is crucial to obtain the desired behaviour of the control loop. For the task of controlling the speed and steering angle of a vehicle driving around a track the following state vector is choosen:
**[x,y,ψ,v]**
- x and y are the position of the vehicle
- ψ is the vehicle orientation with respect to the center lane
- v is the velocity of the vehicle

## Actuators
In order to control the system the MPC needs some output on the system. The available actuators are the steering angle and acceleration. These are summarized in the vector:
**[δ,a]**
- δ is the steering angle, which influences the orientation as well as the position of the vehicle
- a is the acceleration, which influences the vehicle speed as well as the position of the vehicle

## Update equations
The MPC predicts a future state in time (dt) by updating the current state with the planned actuation. Therefore we need the following update equations for the future state [x,y,ψ,v] based on the current state and the actuators [δ,a]:
- x(t+dt) = x(t) + v(t)\*cos(ψ(t))\*dt
- y(t+dt) = y(t) + v(t)\*sin(ψ(t))\*dt
- ψ(t+dt) = ψ(t) + (v(t)/L_f)\*δ(t)\*dt
- v(t+dt) = v(t) + a(t)\*dt

## Cost function
The MPC controls the vehicle by predicting the impact of actuation for N timesteps with a delta-T dt. With a cost function implemented as operator() of the class FG_eval in MPC.cpp the optimal actuations for minimizing the cost are solved. The cost functions and weighting are crucial tuning factors for the MPC algorithm. As contributors to the cost the following was chosen:
- The difference from the reference state where we would like the vehicle to be (center of lane, oriented along lane, driving with set speed)
  - ```double cost_cte;```  Cross-Track error
  - ```double cost_eps;```  Orientation (angle difference to lane)
  - ```double cost_v;```    Velocity difference to set speed
- Usage of actuators and difference between actuations for concurrent timesteps is also associated with a cost to counteract erratic actuation
  - ```double cost_delta_cur;```   Steering actuation (δ)
  - ```double cost_delta_diff;```  Difference of steering actuation (δ(t+dt) - δ(t))
  - ```double cost_a_cur;```       Acceleration actuation (a)
  - ```double cost_a_diff;```      Difference of acceleration actuation (a(t+dt) - a(t))
To further penalize large difference the cost values are squared. With the Ipopt library the equation system for minimizing the cost is implemented in Solve() of the class MPC in mpc.cpp. 

## MPC loop
Everytime the main.cpp receives a new state from the simulator the MPC solver is called. The cost minimal actuations for N timesteps are determined and the actuation values for the first timestep are passed back to the simulator. With the next received state the cycle starts from the beginning.

## Paramter tuning
In order to successfully drive around the track the right values for the horizon (N, dt) and the cost functions (cost_...) have to be chosen. In order for efficient tuning the parameters are stored as config file in /data/parameters.json so that no recompiling is required. The parameters were determined by starting with a relatively small horizon (N=9, dt=0.1). For this horizon an initial set of cost paramters was chosen (cost_cte=100, cost_eps=10, cost_v=10), cost_delta_cur=0, cost_delta_diff=0, cost_a_cur=0, cost_a_diff=0). Starting from there the final set was empirically determined:
```
{
	"reference" : {
		"cte" : 0,		  --> center of lane
		"eps" : 0,			--> orianted along lane
		"v" : 40				--> 40km/h as target speed
		},
	"horizon" : {
		"N" : 15,				--> larger values don't improve driving behaviour, but become more computing intesive
		"dt" : 0.1			--> the starting value of 0.1 was a good fit
		},
	"cost" : {
		"cte" : 3000,   --> Minimizing the cte to keep the vehicle is most important, therefore a very high cost
		"eps" : 400,    --> Minimizing \eps to keep the vehicle in the center for future states is quite important, cost therefore somewhere in the middle
		"v" : 1,   			--> Low cost to keep vehicle speed, for some curves it becomes crucial to deccelerate heavily, cost therefore very low
		"delta_cur" : 1,	--> Current δ needs to have a low cost to keep vehicle on track
		"delta_diff" : 300,  --> A higher cost for the δ difference reduces the erratic driving behaviour of the vehicle
		"a_cur": 1,     --> Low cost for acceleration to deccelerate in curves
		"a_diff" : 10   --> A slightly higher cost for acceleration difference
		}
}
```

With this final set of parameters the MPC can successfully steer the vehicle around the track while mostly keeping the target velocity of 40km/h.

## Delay
In a real vehicle the actuators will have a delay before affecting the vehicle behaviour. To incorporate the delay a vehicle model is used to predict the state where the vehicle will be in 100ms (line 124-140 in main.cpp). When enabling the delay the MPC was still able to drive around the track by using this method, but performance was quite reduced and the deviation from the central lane was higher.
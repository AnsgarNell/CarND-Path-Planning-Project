[//]: # (Image References)
[image1]: ../img/FSM.png

# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

---

## Reflection

### Waypoints
To feed the `next_x_vals` and `next_y_vals` values the Q & A video has been used, and the same code was copied. After doing some tests, and using this approximation as a beginning, further options were implemented.

###Finite States Machine
A FSM has been defined with the following states:



1. `keep_lane`
2. `prepare_change_left`
3. `prepare_change_right`
4. `change_left`
5. `change_right`

This FSM corresponds with the one used in Lesson 4 lecture 8, as it can be seen int he next figure:

![alt text][image1]

### Transition function

In Lesson 4 lecture 10 we can find a pseudocode to implement a transition function between different FSM states. This has been converted to C++ code and used to calculate the `next_state` variable, used for transition.

### Cost function

The transition function needs a cost function in order to know which state would be the best to change to. In few words, it can be described as follows: If there is no car in front of us, keeping the lane must always be the best option. If we find someone who goes slower than us, we check other lanes status. If there are cars that go faster, or no cars at all, we prepare to change to the best of them, also taken into account lanes that are two changes away.

On the other side, if there is a risk involved in changing a lane, we will keep trying until we find a gap, or our lane again is the best option. To avoid this dangerous lane changing, a very high cost is associated in order to never take this action.

Also, we wait until the lane change has been completed to avoid high jerk or the car doing continuous lane changes without finishing any of them.

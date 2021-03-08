# Abstract
In this project, we present a traditional graph based SLAM algorithm on a well-recognized online SLAM dataset called Victoria Park.

A bunch of innovative ideas were raised during the implementation. Some naturally appears when we tried to solve problems and some was because of our curiosity towards the functioning of
the alogirthm. These ideas can be summarized as three main point:

- To see how data association affects the algorithm
In general, graph based SLAM under a simulation environment always assumes known and correct data association, which simplifies the whole process. However,
in this project, a maximum likelihood algorithm based on the feature of landmarks is used to acquire the data association. 

- To increase the speed of our algorithm
The amount of odometry data in the dataset is multiple times greater than the radar data. It means that the information on movement and observation is uneven. If optimization is done
through all nodes, the amount of odometry data will dominate the complexity. However, it is possible to reduce the time using our algorithm which decreases the nodes to be optimized.
This is achieved by merging a series of nodes without observation data provided.

- To add signature for landmark
Signature is a unique symbol that helps tell us the right landmark. If all the landmarks looks exactly the same, it is impossible to find such kind of feature. However, Victoria park
is recorded in a forest with woods and used them as landmarks. It is possible to rougly estimate the size of the wood based on radar data. Hence we add one more dimension, which is the
radius of the wood, as the signture. This shall increase the precision of data association.



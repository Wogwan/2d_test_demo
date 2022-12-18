# svm_irl
Apprenticeship based learning (Inverse Reinforcement Learning)

Implementation based on http://ai.stanford.edu/~ang/papers/icml04-apprentice.pdf

Requirements:
-------------

Convex Optimization Toolbox for Matlab by Stephen P. Boyd


Run Command:
------------
    matlab svm_irl.m
    
Output:
-------
After Doing value iteration, we get a policy (left side) for the original reward (right side), from the below figure

**Note:** While Running you won't see this because ploton is 0 in mdp class because there is another plot of learned reward which is shown below, it changes as I am taking random demonstration from the policy.

#### Policy with Value Iteration, Original Reward: ####
![alt text](https://github.com/srikanthmalla/svm_irl/blob/master/value_iteration%20and%20original%20reward.png )

#### Learnt Reward: #### 
![alt text](https://github.com/srikanthmalla/svm_irl/blob/master/learned_reward.png)


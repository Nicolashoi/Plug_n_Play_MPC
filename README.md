# SP_PnP control for interconnected systems

Plug-and-play (PnP) control schemes are developed where the interconnected
system is decomposed into smaller subsystems and a local controller is developed for each.
When the network topology changes, only the local controllers of a few subsystems are modified depending
on how close the subsystems are to these changes. Particularly, PnP model predictive control
(MPC) has gained attention due to its ability to take state and input constraints into account. PnP MPC
schemes are usually divided into two phases; an offline redesign phase which ensures stability (using passivity theory) and an online phase which ensures feasibility.

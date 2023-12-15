# Gemc-BDX-PIPE
GEMC repository for BDX PIPE

To compile evio2root you need to set the environment BANKS to the evio2root folder: 
e.g. in the .cshrc you need to add 
seteng BANKS /opt/sim/evio2root/ 

To compile everything you need a starting installation of gemc/5.2 to have all the dependencies required. 
In order to compile with root >6.30 you need to change the standard evio dependency by modifying the libsrc++ folder. This folder is usially in the /evio/1.9/src/ folder in the gemc 5.2 installation


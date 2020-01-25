# contact-lines-on-an-elastic-membrane
This  project contains key codes for our paper **"Static interface profiles for contact lines on an elastic membrane"**.

 Folder **key_codes** contains four matlab files combined to be the whole code for our project: **main.m reshape_droplet.m update_membrane.m move_contactline.m**.

**main.m** sets initial parameters of the numerical model and gives the whole structure of the numerical method.

**reshape_droplet.m** is the function to compute the location of the current droplet using the area constraint and contact angles.

**update_membrane.m** is the function to form and solve the main linear system in our numerical method to update the configuration of the current membrane.

**move_contactline.m** is the function to update the contact line position along the membrane.

Folder **data** contains all data used to produce **Figure 3-5** in our paper, please refer to the folder for details.

Folder **tools** contains programs using data to produce figures.

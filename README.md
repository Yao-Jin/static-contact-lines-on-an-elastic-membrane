# contact-lines-on-an-elastic-membrane
This website provides the codes and data for the paper  **"Static interface profiles for contact lines on an elastic membrane"**.

The folder **codes** contains the matlab functions/scripts for the numerical method proposed in the paper:

  **main.m** sets initial parameters of the numerical model and gives the whole structure of the numerical method.

**reshape_droplet.m** is the function to compute the location of the current droplet using the area constraint and contact angles.

**update_membrane.m** is the function to form and solve the main linear system in our numerical method to update the configuration of the current membrane.

**move_contactline.m** is the function to update the contact line position along the membrane.

Folder **data** contains all data used to produce **Figure 3-5** in our paper, please refer to the folder for details.

Folder **tools** contains programs using data to produce figures.

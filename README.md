# venus
Venus' Circumsolar Ring
The code presented here was created and used by me during a summer project at the institute of astronomy, cambridge. Data obtained from STEREO satelliites was reduced to detect the ring, following a recipe given by Jones et al (2013). The following gives a brief description of what each code file does:

1) Jones Recipe_A and Jones Recipe_B : The code that reduces the input data from STEREO images and yields the circumsolar ring. The A and B correspond to STEREO A and STEREO B because there are some small differences in the code.

2) functions : just a collection of some simplified functions that can be imported into a given module and used

3) positions : useful code that indicates whether Venus, Earth, or the Milky Way might be in the field of view of either satellite.

4) orb-elements : a plot of the fractional variation in the orbital elements of the two STEREO satellites.

5) ring-image : a piece of idl code that yields line-of-sight integrated images of Andrew's models of circumsolar disks from the perspective of the two STEREO satellites

6) idl_image : converts the output of the idl ring-image code to graphs and useful images that can be used to compare the images obtained from the model with the data

7) Disk_model_view : An initial naive approach to viewing Andrew's code's output from an edge-on perspective. Does not integrate along the line of sight.

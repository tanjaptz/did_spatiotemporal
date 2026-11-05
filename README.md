# did_spatiotemporal

R script for conducting a spatial and temporal difference-in-differences (DiD) analysis, in this case on dam construction and protest events in Canada.
The script takes a bounding box (e.g. country) divides it in cells for control groups and makes equally sized treatment zones. 
Then it uses this data to call `att_gt()` from the `did` package. 

Example of cells created by the script.
<img width="1200" height="1100" alt="figure_3" src="https://github.com/user-attachments/assets/c173d3f9-09fe-43b2-affe-15a7e8057550" />

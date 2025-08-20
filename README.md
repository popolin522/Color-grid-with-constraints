# Research question
Color a square grid (7 by 7) with three colors while satisfying some rules and ensuring some level of uniformity / symmetry. This stems from one of my side projects about assigning hydrophilic (green) and hydrophobic (blue) chains on a nanocube surface. 
This problem is challenging because it's freaking hard to sample.
# Rules in the script
See the pic attached below for some general ideal.
1. Curvature selectivity. Blue color (or hydrophobic species) should be on the vertex first, edge, then face.
2. Region-specific density. Some regions can only have a certain number of colors (green + blue).
3. Blue ratio. total number of blue on the grid cannot exceed some number.
4. Be symmetric-ish. For each coloring move, aim for symmetry—color the next point as opposite as possible to the last (like opposite edge).
<img width="468" height="267" alt="image" src="https://github.com/user-attachments/assets/0705a32d-91d0-41be-9b08-068b2835f114" />

# Improper (or engineer) coloring solution
I am smart enough to realize this is basically a coloring problem. Apparently not smart enough to solve it in a mathematically glorious way. So here's the coding solution most likely to offend mathematicians. 
For my math friends please check out this note to learn more about proper q-coloring https://www.math.tau.ac.il/~peledron/homepage_files/Proper%20colorings%20IMU%20meeting.pdf
You will not learn it here.

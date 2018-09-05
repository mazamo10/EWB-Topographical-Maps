# EWB-Topographical-Maps

Hello, and welcome to the irrigation team. The purpose of this code is to convert data gathered in the field into Cartesian coordinates, and then to use those coordinates to draw contour lines in AutoCAD so that the resulting maps can be edited.

In order to locate the instrument in space, two reference points located at the same height are needed. The distance between the instrument and these reference points is measured along with the angle which the instrument sweeps through to each point. This data is then used to triangulate the instrument’s location.

Once the instruments location has been found with respect to the reference points, the remaining data (distance, angle, height) is imported from excel and converted into (x,y,z) tuples. Finally, these points are triangulated through the use of a Delaunay triangulation algorithm and contour lines are created iteratively with the eventual goal of these lines being drawn in AutoCAD (note that in the initial version of the code, the commands which move the mouse will be incorrect, as every display has different coordinate points).

Currently this code is written in Python version 3.6 and includes the following libraries
•	math
•	time
•	pandas
•	pyautogui
•	numpy
•	scipy.optimize
•	matplotlib.pyplot
•	scipy.spatial import Delaunay

All of these libraries, except pyautogui, can be accessed by downloading Anaconda https://www.anaconda.com/download/, otherwise use the pip command https://docs.python.org/3/installing/index.html. Additionally, please install AutoCAD since the final drawing will be made there.

Thank you,
Matthew

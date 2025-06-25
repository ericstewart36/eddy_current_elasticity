// Gmsh project created on Tue Apr  8 10:18:40 2025
SetFactory("OpenCASCADE");
//
Mesh.RecombineAll = 0;     // Prevent quads
Mesh.Recombine3DAll = 0;   // Prevent hexes
Mesh.Algorithm = 6;        // 2D triangle meshing (e.g., frontal-Delaunay)
Mesh.Algorithm3D = 1;      // 3D tetrahedral meshing (Delaunay)
Mesh.ElementOrder = 1;     // Linear elements
//
//+
lf = 15;
lc = 15;
//+
Point(1) = {0, 0, 0, lf};
//+
Point(2) = {76, 0, 0, lf};
//+
Point(3) = {487, 0, 0, lf};
//+
Point(4) = {487, 0, 100, lf};
//+
Point(5) = {76, 0, 100, lf};
//+
Point(6) = {0, 0, 100, lf};
//+
Point(7) = {243.5, 0, 0, lf};
//+
Point(8) = {243.5, 0, 50, lf};
//+
Point(9) = {243.5, 0, 100, lf};
//+
// // //+
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 1};
//+
Line(5) = {2, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 5};
//+
Line(9) = {7, 3};
//+
Line(10) = {3, 4};
//+
Line(11) = {4, 9};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {9, 10, 11, -7, -6};
//+
Plane Surface(3) = {3};
//+
Extrude {0, 3.175, 0} {
  Surface{3}; Surface{2}; Surface{1};  Layers {3}; 
}
// Define a box for the air domain
Box(9)={-50,-50,-50, (637-50),100,200};
//+
BooleanDifference(10) = { Volume{9};}{ Volume{1}; Volume{2}; Volume{3}; };
//+
MeshSize {19,20,21,22,23,24,25,26} = lc;
//+
Physical Volume("Free_vol", 44) = {1, 2};
//+
Physical Volume("Clamp_vol", 45) = {3};
//+
Physical Volume("Air_vol", 46) = {10};
//+
Physical Surface("Clamp_surfs", 47) = {17, 15, 14, 16, 1};
//+
Physical Surface("Air_surfs", 48) = {21, 18, 22, 23, 19, 20};
//+
Physical Surface("current_surf", 49) = {7};
//+
Physical Surface("maxwell_surf", 50) = {13, 2, 11, 10, 9, 3, 6, 4, 5};
// //+
// Gmsh project created on Tue Mar 11 15:43:45 2025
SetFactory("OpenCASCADE");
lf = 3;
lc = 50;
//+
Point(2) = {27, 0, 0, lf};
//+
Point(3) = {55, 0, 0, lf};
//+
Point(4) = {80, 0, 0, lf};
//+
Point(5) = {95, 0, 0, lf};
//+
Point(6) = {95, -52, 0, lf};
//+
Point(7) = {80, -52, 0, lf};
//+
Point(8) = {55, -52, 0, lf};
//+
Point(9) = {27, -52, 0, lf};
//+
Point(11) = {0, 3.8, 0, lf};
//+
Point(12) = {0, 6.8, 0, lf};
//+
Point(13) = {65, 6.8, 0, lf};
//+
Point(14) = {65, 3.8, 0, lf};
//+
Point(15) = {0, 200, 0, lc};
//+
Point(16) = {300, 200, 0, lc};
//+
Point(17) = {300, -250, 0, lc};
//+
Point(18) = {0, -250, 0, lc};
//+
Line(1) = {11, 12};
//+
Line(2) = {12, 13};
//+
Line(3) = {14, 13};
//+
Line(4) = {14, 11};
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 2};
//+
Line(9) = {4, 5};
//+
Line(10) = {5, 6};
//+
Line(11) = {6, 7};
//+
Line(12) = {7, 4};
//+
Line(13) = {12, 15};
//+
Line(14) = {15, 16};
//+
Line(15) = {16, 17};
//+
Line(16) = {17, 18};
//+
Line(17) = {18, 11};
//+
Curve Loop(1) = {14, 15, 16, 17, -4, 3, -2, 13};
//+
Curve Loop(2) = {6, 7, 8, 5};
//+
Curve Loop(3) = {10, 11, 12, 9};
//+
Plane Surface(1) = {1, 2, 3};
//+
Curve Loop(4) = {6, 7, 8, 5};
//+
Plane Surface(2) = {4};
//+
Curve Loop(5) = {9, 10, 11, 12};
//+
Plane Surface(3) = {5};
//+
Curve Loop(6) = {4, 1, 2, -3};
//+
Plane Surface(4) = {6};
//+
Physical Surface(18) = {1};
//+
Physical Surface(19) = {4};
//+
Physical Surface(20) = {3};
//+
Physical Surface(21) = {2};
//+
//Physical Curve(22) = {14};
//+
//Physical Curve(23) = {15};
//+
//Physical Curve(24) = {16};
//+
//Physical Curve(25) = {5, 8, 6, 7};
//+
//Physical Curve(26) = {9, 12, 10, 11};
//+
Physical Curve(27) = {2, 3, 4};
//+
Physical Curve(28) = {13, 17, 1};
//+
Physical Curve(29) = {14, 15, 16};
//+
Physical Curve(30) = {5, 8, 7, 6, 9, 12, 11, 10};
//+
Field[1] = Box;
//
Field[1].VIn = lf/3;
//+
Field[1].VOut = lc;
//+
Field[1].XMax = 70;
//+
Field[1].XMin = 0;
//+
Field[1].YMax = 15;
//+
Field[1].YMin =0;
//+
Field[1].ZMax = 0;
//+
Field[1].ZMin = 0;
//+
Background Field = 1;
//+
// MeshSize {12,11,13,14} = lf/2;

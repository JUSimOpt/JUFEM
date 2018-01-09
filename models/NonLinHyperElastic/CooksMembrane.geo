/* Script to generate Cooks Membrane */


Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {1, 0.5, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};

Point(5) = {0, 0, 1, 1.0};
Point(6) = {0, 1, 1, 1.0};
Point(7) = {1, 0.5, 1, 1.0};
Point(8) = {1, 1, 1, 1.0};
//+
Line(1) = {5, 7};
//+
Line(2) = {7, 8};
//+
Line(3) = {8, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {1, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 2};
Line(8) = {2, 1};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {4, 8};
Line(12) = {7, 3};


Transfinite Line{14,16,18,20,8,10} = 3 Using Progression 1; // Create 3 points on these lines
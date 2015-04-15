
include <a_piece.scad>

//////////////////////////////////////////////
// Piece B
//////////////////////////////////////////////
Platform         = 10;
MountX           = RailX+Platform;
ScrewGuard       = 2;
ScrewHeadDiameter= 4;

module B_Profile(tolerance=0)
{
	rotate([0,-90,0])
	linear_extrude(height=HolderSpace+tolerance,center=true)
	polygon(
		[
			[0,0],
			[MountX,0],
			[MountX,HolderBulk+HolderY+RailY],
			[RailX,HolderBulk+HolderY+RailY],
			[RailX,HolderBulk+HolderY],
			[0,HolderBulk+HolderY]
		]	
	);
}

// To stick with A_Piece()
module B_Screw1()
{
	translate([0,HolderBulk/2,-1])
	cylinder(h=ScrewGuard+2,d=ScrewHoleDiameter+0.5,$fn=Resolution);
	translate([0,HolderBulk/2,ScrewGuard])
	cylinder(h=MountX-ScrewGuard+1,d=ScrewHeadDiameter,$fn=Resolution);
}	

// To Attach servo
ServoScrewDiameter=1;

module B_Screw2()
{
	translate([0,HolderBulk+HolderY+1,2])
	rotate([90,0,0]) cylinder(h=5+1,d=ServoScrewDiameter,$fn=Resolution);
}

// To Attach BackBone
AxisDiameter=1;

module B_Screw3()
{
	translate([0,HolderBulk+HolderY+RailY+1,RailX+Platform/2])
	rotate([90,0,0]) 
	cylinder(h=5+1,d=AxisDiameter,$fn=Resolution);
}


// To be printed with tolerance=-1
module B_Piece(tolerance=0)
{
	difference()
	{
		B_Profile(tolerance);
		B_Screw1();
		B_Screw2();
		B_Screw3();
	}
}




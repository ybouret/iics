Resolution=60;

module FishHead()
{
	union()
	{
		import("fish_head.stl");
		import("fish_junc.stl");
	}	
}

module FishTail()
{
	difference()
	{
		render() import("fish_tail.stl");
		import("fish_junc.stl");
	}	
}

include <a_piece.scad>
include <b_piece.scad>
include <servo.scad>


module AB_Hull_Old()
{
	rotate([0,-90,0])
	linear_extrude(height=HolderSpace,center=true)
	polygon(
		[
			[MountX,0],
			[MountX,HolderBulk+ServoY],
			[-(HolderX+HolderBulk),HolderBulk+ServoY],
			[-(HolderX+HolderBulk),0]
		]	
	);
}

module AB_Hull()
{
	hull()
	{
		A_Piece();
		rotate([0,-90,0])
		linear_extrude(height=HolderSpace,center=true)
		polygon(
			[
				[MountX,0],
				[MountX,HolderBulk+ServoY],
				[0,HolderBulk+ServoY],
				[0,0]
			]	
		);
	}
}


//////////////////////////////////////////////
// BackBone FingerPrint
//////////////////////////////////////////////
BackBoneDiameter=18;
BackBoneAttachDiameter=BackBoneDiameter+1;
module BackBoneAttach()
{
	color("yellow")
	translate([0,HolderBulk+HolderY+RailY+BackBoneAttachDiameter/2,RailX])
	cylinder(h=Platform,d=BackBoneAttachDiameter,$fn=Resolution);
}


//////////////////////////////////////////////
// Move In Face of the backbone
//////////////////////////////////////////////
ToOrigin=[0,-(HolderBulk+HolderY+RailY+BackBoneDiameter/2),-MountX];

//////////////////////////////////////////////
// Negative Space
//////////////////////////////////////////////
module BigHole()
{
	translate(ToOrigin)
	union()
	{
		AB_Hull();

		hull()
		{
			ServoSpace();
			BackBoneAttach();
		}
	}

}




TubaDiameter   = 5;
TubaTorus      = HolderBulk;
TubaHeight     = 100;

module TubaCut()
{
	TubaX = 2*(TubaTorus+TubaDiameter/2)+2;
	TubaY = TubaX/2;
	TubaZ = TubaDiameter+2;
	translate([0,TubaY/2,0]) cube([TubaX,TubaY,TubaZ],center=true);
}

// Start Tuba@Origin
module TubaLink()
{
	translate([0,TubaTorus,-TubaTorus])
	rotate([-90,180,0])
	union()
	{
		cylinder(h=TubaHeight,d=TubaDiameter,$fn=Resolution);
		translate([0,TubaTorus,0])
		rotate([0,-90,0])
		difference()
		{
			rotate_extrude(convexity=10,$fn=Resolution)
			translate([TubaTorus,0,0]) circle(d=TubaDiameter,$fn=Resolution);
			TubaCut();
			rotate([0,0,-90]) TubaCut();
		}
	}
}

module Tuba()
{
	translate([0,HolderBulk+TubaDiameter/2+0.5,-HolderX]) TubaLink();
}


//////////////////////////////////////////////
// Carved Head
//////////////////////////////////////////////
module CarvedHead(zoom=1)
{
	difference()
	{
		union()
		{
			difference()
			{
				scale(zoom) FishHead();
				BigHole();
			}
			translate(ToOrigin) A_Piece();
		}
		translate(ToOrigin) Tuba();
	}
}

//CarvedHead(2.4);
/*
%scale(1.6)
{ 
	//FishHead();
	FishTail();
}
translate([0,0,MountX])
{
translate(ToOrigin)
{
	A_Piece();
	color("red") B_Piece();
	BackBoneAttach();
	Tuba();
	Servo();
	ServoSpace();
}
	cylinder(h=50,r1=5,r2=2,$fn=60);
}
*/
difference()
{
	scale(1.6) FishTail();
	translate([0,0,MountX]) BigHole();
	cylinder(h=50,r1=5,r2=2,$fn=60);
}


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


module AB_Hull()
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




//////////////////////////////////////////////
// BackBone FingerPrint
//////////////////////////////////////////////
BackBoneDiameter=18;
BackBoneAttachDiameter=BackBoneDiameter+1;
module BackBoneAttach()
{
	color("yellow")
	translate([0,HolderBulk+HolderY+RailY+BackBoneAttachDiameter/2,RailX])
	cylinder(h=Platform,d=BackBoneAttachDiameter,$fn=30);
}


//////////////////////////////////////////////
// Move In Face of the backbone
//////////////////////////////////////////////
MoveToOrigin=[0,-(HolderBulk+HolderY+RailY+BackBoneDiameter/2),-MountX];

module BigHole()
{
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

translate(MoveToOrigin)
{
//BigHole();
A_Piece();
B_Piece(tolerance=-1);
BackBoneAttach();
Servo();

}
scale(2.4) 
{
	FishTail();
   %FishHead();
}

v=[0,-40,-77];

module Toto()
{

union()
{
	difference()
	{
		scale(2.5) FishHead();
		translate(v)
		BigHole();
	}

	translate(v)
	{
		A_Piece();
	}
}
}

//%Toto(); 

translate(v) 
{ 
	//Servo(); 
	//translate([0,0,50]) B_Piece();
}

//translate([0,0,20]) scale(2.5) FishTail();

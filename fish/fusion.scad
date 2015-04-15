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


//////////////////////////////////////////////
// Carved Head
//////////////////////////////////////////////
module CarvedHead(zoom=1)
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
}

A_Piece();
Tuba=5;


//CarvedHead(zoom=2.4);
//translate(ToOrigin) Servo();
//translate(ToOrigin) B_Piece(tolerance=-1);


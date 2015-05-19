Resolution = 30;

CraddleBodyDepth  = 24;
CraddleBodyHeight = 16;
CraddleBodyWidth  = 13;

CraddleSlabHeight = 2;
CraddleSlabEdge   = 5;
CraddleSlabDepth  = CraddleBodyDepth+2*CraddleSlabEdge;
CraddleSlabWidth  = CraddleBodyWidth+1;

CraddleBaseHeight = CraddleBodyHeight+CraddleSlabHeight;

CraddleHoodWidth  = CraddleBodyWidth;
CraddleHoodHeight = 12;
CraddleHoodDepth  = CraddleBodyDepth;

CraddleTopHeight  = CraddleBaseHeight + CraddleHoodHeight;

ArmBase   = CraddleTopHeight;
ArmRadius = 13;
ArmHeight = 8.5;
ArmOffset = -17.5+CraddleBodyDepth/2;

//////////////////////////////////////////////////////////////////////
//
// Fish Head/Tail and Junction
//
//////////////////////////////////////////////////////////////////////

//____________________________________________________________________
//
// unscaled Fish Junction
//____________________________________________________________________
module FishJunction()
{
	render()
	translate([0,0,-1])
	scale([10,25,1])
	linear_extrude(height=10,scale=0.95) 
	{
		difference()
		{
			circle(r=1,$fn=Resolution);
			circle(r=0.95,$fn=Resolution);
		}
	}
}

//____________________________________________________________________
//
// Fish Head: union of STL and Junction (then scaled)
//____________________________________________________________________
module FishHead(zoom=1)
{
	render()
	scale(zoom)
	{
		union()
		{
			import("fish_head.stl");
			FishJunction();
		}
	}
}

//____________________________________________________________________
//
// Fish Tail
//____________________________________________________________________
module FishTail(zoom=1)
{
	scale(zoom)
	import("fish_tail.stl");
}

//FishHead(1.42);
//%FishTail(1.42);

//////////////////////////////////////////////////////////////////////
//
// Servo
//
//////////////////////////////////////////////////////////////////////
module Servo()
{
	color("orchid") rotate([0,90,0]) scale(10) import("data/microservo.stl");
}

//////////////////////////////////////////////////////////////////////
//
// Arm
//
//////////////////////////////////////////////////////////////////////
module Arm()
{
	translate([0,ArmBase,ArmOffset])
	rotate([-90,0,0])
	import("data/custom_arm.stl");
}

module ArmSpace()
{
	render()
	translate([0,ArmBase-1,ArmOffset])
	hull()
	{
		rotate([-90,0,0]) cylinder(h=ArmHeight+1,r=ArmRadius,$fn=Resolution);
		translate([-ArmRadius,0,-ArmOffset+CraddleBodyDepth/2-2])
		cube([2*ArmRadius,ArmHeight+1,2]);
	}
	//render()
//	translate([0,ArmBase,0])
//	hull()
//	{
//		translate([0,0,ArmOffset])
//		rotate([-90,0,0])
//		cylinder(h=ArmHeight,r=ArmRadius,$fn=Resolution);
//		translate([0,ArmHeight/2,0]) cube([2*ArmRadius,ArmHeight,2],center=true);
//	}
}

//////////////////////////////////////////////////////////////////////
//
// Craddle
//
//////////////////////////////////////////////////////////////////////
module CraddleBody()
{
	render()
	translate([-CraddleBodyWidth/2,0,-CraddleBodyDepth/2])
	cube([CraddleBodyWidth,CraddleBodyHeight+CraddleSlabHeight/4,CraddleBodyDepth]);
}

module CraddleSlab()
{
	render()
	translate([-CraddleSlabWidth/2,CraddleBodyHeight,-CraddleSlabDepth/2])
	cube([CraddleSlabWidth,CraddleSlabHeight,CraddleSlabDepth]);
}

module CraddleHood()
{
	translate([-CraddleHoodWidth/2,CraddleBaseHeight-CraddleSlabHeight/4,-CraddleHoodDepth/2])
	cube([CraddleHoodWidth,CraddleHoodHeight+CraddleSlabHeight/4,CraddleHoodDepth]);
}

module Craddle()
{
	render()
	union()
	{
		CraddleBody();
		CraddleSlab();
		CraddleHood();
		ArmSpace();
	}
}



%FishHead(1.42);
Servo();
%Craddle();
//Arm();
ArmSpace();


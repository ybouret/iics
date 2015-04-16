Resolution=60;

////////////////////////////////////////////////////////////////////////////////
//
// Computed Data: default=100mm, 
// Z-Center at the Head/Tail Ratio
//
////////////////////////////////////////////////////////////////////////////////

module FishHead(zoom=1)
{
	scale(zoom)
	{
		union()
		{
			import("fish_head.stl");
			import("fish_junc.stl");
		}	
	}
}

module FishTail(zoom=1)
{
	scale(zoom)
	{
		difference()
		{
			import("fish_tail.stl");
			import("fish_junc.stl");
		}
	}
}



////////////////////////////////////////////////////////////////////////////////
//
// Negative Space for the Servo: Craddle Body
//
////////////////////////////////////////////////////////////////////////////////

CraddleBodyDepth  = 24;
CraddleBodyHeight = 16;
CraddleBodyWidth  = 13;

module CraddleBody()
{
	craddle_body_extra=1;
	
	translate(
		[
		0,
		CraddleBodyHeight/2,
		-CraddleBodyDepth/2+craddle_body_extra/2
		])
    cube([
		CraddleBodyWidth,
		CraddleBodyHeight+1,
		CraddleBodyDepth+craddle_body_extra],center=true);		
}

////////////////////////////////////////////////////////////////////////////////
//
// A Servo to play with
//
////////////////////////////////////////////////////////////////////////////////
module Servo()
{
	translate([0,0,-CraddleBodyDepth/2])
	rotate([0,90,0]) scale(10) import("data/microservo.stl");
}

////////////////////////////////////////////////////////////////////////////////
//
// Negative Space for the Servo: Craddle Rail
//
////////////////////////////////////////////////////////////////////////////////
CraddleRailHeight = 2;
CraddleRailExpand = 5;
CraddleRailDepth  = CraddleBodyDepth+2*CraddleRailExpand;

module CraddleRail()
{
	translate([0,CraddleBodyHeight,-CraddleBodyDepth/2])
	translate([0,CraddleRailHeight/2,0])
	cube([CraddleBodyWidth+0.5,CraddleRailHeight,CraddleRailDepth],center=true);
}

////////////////////////////////////////////////////////////////////////////////
//
// Negative Space for the Servo: Craddle Top, to be able to insert servo
//
////////////////////////////////////////////////////////////////////////////////
CraddleTopBase    = CraddleBodyHeight+CraddleRailHeight;
CraddleTopHeight = 12;
CraddleTopDepth  = CraddleBodyDepth+1;

module CraddleTop()
{
	craddle_top_extra=1;
	translate(
		[
			0,
			CraddleTopHeight/2+CraddleTopBase-0.5,
			-CraddleBodyDepth/2+craddle_top_extra/2
		]
		)
	cube([CraddleBodyWidth,CraddleTopHeight+1,CraddleTopDepth+craddle_top_extra],center=true);
}


////////////////////////////////////////////////////////////////////////////////
//
// Negative Space for the Arm
//
////////////////////////////////////////////////////////////////////////////////
ArmBase   = CraddleBodyHeight+CraddleRailHeight+CraddleTopHeight-4;
ArmRadius = 10;
ArmHeight = 8;
ArmOffset = -16;

module ArmSpace()
{
	translate([0,ArmBase,0])
	hull()
	{
	translate([0,0,ArmOffset])
	rotate([-90,0,0])
	cylinder(h=ArmHeight,r=ArmRadius,$fn=Resolution);
	translate([0,ArmHeight/2,0]) cube([2*ArmRadius,ArmHeight,2],center=true);
	}
}

////////////////////////////////////////////////////////////////////////////////
//
// Screw Hole for Platform
//
////////////////////////////////////////////////////////////////////////////////
PlatformScrewDiameter=2;
PlatformScrewDepth   =5;
PlatformScrewOffset  =PlatformScrewDiameter+2;

module PlatformScrew()
{
	translate([0,-PlatformScrewOffset,-PlatformScrewDepth])
	cylinder(h=PlatformScrewDepth+1,d=PlatformScrewDiameter,$fn=Resolution);
}

////////////////////////////////////////////////////////////////////////////////
//
// Tuba For Wiring
//
////////////////////////////////////////////////////////////////////////////////

TubaDiameter   = 5;
TubaTorus      = 20;
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
	translate([0,TubaDiameter/2+1,-CraddleBodyDepth+TubaTorus/3])
	color("green") TubaLink();
}

module Craddle()
{
	union()
	{
		CraddleBody();
		CraddleRail();
		CraddleTop();
		ArmSpace();
		PlatformScrew();
	}
}

Craddle();


////////////////////////////////////////////////////////////////////////////////
//
// Platform To Fix BackBone And Attach Servo
//
////////////////////////////////////////////////////////////////////////////////



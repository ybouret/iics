Resolution=60;

////////////////////////////////////////////////////////////////////////////////
//
// Computed Data: default=100mm, 
// Z-Center at the Head/Tail Ratio
//
////////////////////////////////////////////////////////////////////////////////

module FishHead(zoom=1)
{
	render()
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
	render()
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
	render()
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
	render()
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

module CraddleRail(expand=0)
{
	render()
	translate([0,CraddleBodyHeight,-CraddleBodyDepth/2])
	translate([0,CraddleRailHeight/2+expand/2,0])
	cube([CraddleBodyWidth+0.5,CraddleRailHeight+expand,CraddleRailDepth],center=true);
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
	render()
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
	render()
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
	render()
	translate([0,-PlatformScrewOffset,-PlatformScrewDepth])
	cylinder(h=PlatformScrewDepth+1,d=PlatformScrewDiameter,$fn=Resolution);
}

////////////////////////////////////////////////////////////////////////////////
//
// Tuba For Wiring
//
////////////////////////////////////////////////////////////////////////////////

TubaDiameter   = 5;
TubaTorus      = 18;
TubaHeight     = 60;

module TubaCut()
{

	TubaX = 2*(TubaTorus+TubaDiameter/2)+2;
	TubaY = TubaX/2;
	TubaZ = TubaDiameter+2;
	render()
	translate([0,TubaY/2,0]) cube([TubaX,TubaY,TubaZ],center=true);
}

// Start Tuba@Origin
module TubaLink()
{
	render()
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


// The Tuba, to be subtracted from FishHead-Craddle

module Tuba()
{
	render()
	translate([0,TubaDiameter/2+1,-CraddleBodyDepth+TubaTorus/3])
	color("green") TubaLink();
}


////////////////////////////////////////////////////////////////////////////////
//
// The Craddle, to be subtracted from FishHead
//
////////////////////////////////////////////////////////////////////////////////
module Craddle()
{
	render()
	union()
	{
		CraddleBody();
		CraddleRail();
		CraddleTop();
		ArmSpace();
		PlatformScrew();
	}
}

//%Craddle();


////////////////////////////////////////////////////////////////////////////////
//
// Platform To Fix BackBone And Attach Servo => EdgeBody
//
////////////////////////////////////////////////////////////////////////////////

PlatformDepth = 10;
EdgeDepth     = PlatformDepth+CraddleRailExpand;
EdgeHeight    = 8;


module EdgeContour(extra=0)
{
		render()
		hull()
		{
			union()
			{
				translate([0,-PlatformScrewOffset,-extra])
				cylinder(h=EdgeDepth+extra,d=CraddleBodyWidth,$fn=Resolution);
				
				translate([0,CraddleBodyHeight+CraddleRailHeight-EdgeHeight,-extra])
				translate([0,EdgeHeight/2,(EdgeDepth+extra)/2])
				cube([CraddleBodyWidth,EdgeHeight,EdgeDepth+extra],center=true);
			}
		}
}

module EdgeBody()
{
	render()
	difference()
	{
		EdgeContour(extra=0);
		CraddleRail(expand=1);
	}
}

////////////////////////////////////////////////////////////////////////////////
//
// Hole to attach EdgeBody to Craddle => EdgeAttach
//
////////////////////////////////////////////////////////////////////////////////

EdgeGuardDepth=2;
EdgeGuardDiameter=2;
ScrewHeadDiameter=5;


module EdgeAttach()
{
	render()
	union()
	{
		translate([0,-PlatformScrewOffset,-1])
		cylinder(h=EdgeGuardDepth+2,d=EdgeGuardDiameter,$fn=Resolution);
		translate([0,-PlatformScrewOffset,EdgeGuardDepth])
		cylinder(h=EdgeDepth+1-EdgeGuardDepth,d=ScrewHeadDiameter,$fn=Resolution);
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// Hole to attach EdgeBody to Servo
//
////////////////////////////////////////////////////////////////////////////////

ServoAttachDiameter=1;
ServoAttachDepth   =8;

module ServoAttach()
{
	render()
	translate([0,CraddleBodyHeight-ServoAttachDepth,2])
	rotate([-90,0,0])
	cylinder(h=ServoAttachDepth+1,d=ServoAttachDiameter,$fn=Resolution);
}


////////////////////////////////////////////////////////////////////////////////
//
// Hole to attach EdgeBody to BackBone
//
////////////////////////////////////////////////////////////////////////////////
BBAttachDiameter=2;
BBAttachDepth   =10;

module BBAttach()
{
	render()
	translate([0,CraddleTopBase-BBAttachDepth,CraddleRailExpand+PlatformDepth/2])
	rotate([-90,0,0])
	cylinder(h=BBAttachDepth+1,d=BBAttachDiameter,$fn=Resolution);
}


////////////////////////////////////////////////////////////////////////////////
//
// The Edge, to be printed
//
////////////////////////////////////////////////////////////////////////////////
module Edge()
{
	render()
	difference()
	{
		EdgeBody();
		EdgeAttach();
		ServoAttach();
		BBAttach();
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// BackBone parameters
//
////////////////////////////////////////////////////////////////////////////////
BBDiameter=12;
BBTip     =4;

BBSpace   = 1;
BBDepth   = 45;

BBHeight  = CraddleTopBase+BBDiameter/2;

ToOrigin  = [0,-BBHeight,0];

module BBSocket()
{
	d_ini = BBDiameter+BBSpace;
	d_end = BBTip     +BBSpace;
	render()
	hull()
	{
		
		translate([0,BBHeight,-1]) cylinder(h=EdgeDepth+1,d=d_ini);
		
		translate([0,BBHeight,EdgeDepth])
		 linear_extrude(height=BBDepth+BBSpace,scale=d_end/d_ini,convexity=10)
		 circle(d=d_ini,$fn=Resolution);
	}
}

module BBCraddle()
{
	render()
	union()
	{
		EdgeContour(extra=1);
		BBSocket();
	}
}

////////////////////////////////////////////////////////////////////////////////
//
// Carved Tail
//
////////////////////////////////////////////////////////////////////////////////

module CarvedTail(zoom=1)
{
	render()
	difference()
	{
		FishTail(zoom);
		translate(ToOrigin) BBCraddle();
	}
}

module TryTail(zoom=1)
{
	%FishTail(zoom);
	translate(ToOrigin) BBCraddle();
}


////////////////////////////////////////////////////////////////////////////////
//
// Carved Head
//
////////////////////////////////////////////////////////////////////////////////

module CarvedHead(zoom=1)
{
	difference()
	{
		render()
		difference()
		{
			FishHead(zoom);
			translate(ToOrigin) Craddle();
		}
		translate(ToOrigin) Tuba();
	}
}


module TryHead(zoom=1)
{
	{
		%FishHead(1.5);
		translate(ToOrigin) color("orange")     Craddle();
		translate(ToOrigin) color("lightgreen") Tuba();
		translate(ToOrigin) color("red")      Edge();
	}
}


TryHead(1.5);



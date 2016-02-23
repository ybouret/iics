Resolution=60;

Height = 15;

SPACE = 0.5;

PlugDiameter  = 10;
PlugRadius    = PlugDiameter/2;

SocketDiameter= 20;
SocketRadius  = SocketDiameter/2;
SocketLid     = 3;

Alpha=45;

HoleOffset=SocketLid;
HoleHeight=Height-(2*SocketLid);

BeamWidth = 0.8 * PlugRadius;
BeamDepth = 30;

GuideDiameter=2;
AttachDiameter=3;


CraddleBodyDepth  = 24;
CraddleBodyHeight = 16;
CraddleBodyWidth  = 13;

CraddleRailHeight = 2;
CraddleRailExpand = 5;
CraddleRailDepth  = CraddleBodyDepth+2*CraddleRailExpand;

CraddleTopBase    = CraddleBodyHeight+CraddleRailHeight;

ArmBase   = CraddleBodyHeight+CraddleRailHeight+CraddleTopHeight-4;
ArmRadius = 13;
ArmHeight = 8.5;
ArmOffset = -17.5;

////////////////////////////////////////////////////////////////////////////////
//
// A Servo to play with
//
////////////////////////////////////////////////////////////////////////////////
module Servo()
{
	{
		translate([0,0,-CraddleBodyDepth/2])
		rotate([0,90,0]) scale(10) import("data/microservo.stl");
		
		translate([0,CraddleTopBase+8,ArmOffset])
		translate([0,4,0])
		rotate([-90,0,0])
		import("data/custom_arm.stl");
	}
}


//Make a hole in a Head
module SocketCut()
{
	Forward=SocketDiameter/cos(Alpha);
	ExtendB=(BeamWidth+SPACE)/2;

	translate([0,HoleOffset,0])
	rotate([-90,-90,0])
	linear_extrude(height=HoleHeight)
	{
		
		union()
		{
			polygon(
				[
					[0,-ExtendB],
					[0,ExtendB],
					[Forward*cos(Alpha),ExtendB+Forward*sin(Alpha)],
					[Forward*cos(Alpha),-(ExtendB+Forward*sin(Alpha))],
				]);
			circle(d=PlugDiameter+SPACE,$fn=Resolution);
		}
	}

}


module OneGuide()
{
	render() translate([0,Height/2,-SocketDiameter/2-1])
	cylinder(h=SocketDiameter+2,d=GuideDiameter,$fn=Resolution);
}

module Guides()
{
	GuidePos = SocketDiameter/2-GuideDiameter/2-1;
	render()
	union()
	{
		render() translate([ GuidePos,0,0]) OneGuide();
		render() translate([-GuidePos,0,0]) OneGuide();
	}
}

//The Pocket
module SocketPocket()
{
	render()
	difference()
	{
		render() rotate([-90,0,0])
		cylinder(h=Height,d=SocketDiameter,$fn=Resolution);
	
		render() SocketCut();
		//render() Guides();
	}
}

module DefaultSocket()
{
	render()
	difference()
	{
		render() SocketPocket();
		render() Guides();
	}
}


// Link The pocket with a beam
module Bone()
{
	render()
	difference()
	{
		render()
		union()
		{
			render()
			translate([0,0,BeamDepth])
			DefaultSocket();
	
			render()
			translate([-BeamWidth/2,0,0])
			cube([BeamWidth,Height,BeamDepth-(PlugRadius+SocketRadius)/2]);
	
			render()
			rotate([-90,0,0])
			cylinder(h=Height,d=PlugDiameter,$fn=Resolution);
		}
		
				
		render()
		translate([0,-1,0])
		rotate([-90,0,0])
		cylinder(h=SocketLid+1+SPACE,d=SocketDiameter+2*SPACE,$fn=Resolution);
		
		render()
		translate([0,Height-(SocketLid+SPACE),0])
		rotate([-90,0,0])
		cylinder(h=SocketLid+1+SPACE,d=SocketDiameter+2*SPACE,$fn=Resolution);
	
		render() Guides();
	}
}

module EndBone()
{
	render()
	union()
	{
		render() Bone();
		render() translate([0,0,BeamDepth+SocketRadius-AttachDiameter/2-2])
		translate([0,SocketLid/2,0])
		rotate([-90,0,0])
		cylinder(h=Height-SocketLid,d=AttachDiameter,$fn=Resolution);
	}
}

module Collar()
{
	DefaultSocket();
}

//Collar();
//
//Bone();
//
//EndBone();


//translate([0,-CraddleBodyHeight,-SocketRadius]) Servo();

translate(ToOrigin)
{
SocketPocket();
rotate([0,30,0])
{
	color("orchid") Bone();
	translate([0,0,BeamDepth]) rotate([0,-30,0])
	{
		EndBone();
	}
}
}

//render() difference()
//{
//union()
//{
//	SocketPocket();
//	rotate([0,20,0])
//	{
//		Bone();
//		//translate([0,0,BeamDepth])
//		//rotate([0,40,0]) color("orchid") EndBone();
//	}
//}
//
//translate([-50,Height/2,-20])
//cube([100,Height,100]);
//}


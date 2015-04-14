module Holder()
{
	translate([7,-50-20,-118])
	rotate([0,-90,0])
		import("ServoHolderRound.stl");
}

/*
union()
{
	difference()
	{
		//scale([2,2,2])  import("fish_head.stl");
		//Holder();
	}
	//Holder();
}
*/

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

//////////////////////////////////////////////
// Piece A
//////////////////////////////////////////////
HolderBulk=10;
HolderX   =24;
HolderY   =16;
RailX     =5;
RailY     =2;
ServoY    =30;
HolderSpace=13;

ScrewHoleDiameter=1;
ScrewHoleDepth   =10;

module A_Profile()
{
	rotate([0,-90,0])
	linear_extrude(height=HolderSpace,center=true)
	polygon([ 
		[0,0], 
		[0,HolderBulk], 
		[-HolderX,HolderBulk], 
		[-HolderX,HolderBulk+HolderY], 
		[-(HolderX+RailX),HolderBulk+HolderY],
		[-(HolderX+RailX),HolderBulk+HolderY+RailY],
		[-(HolderX), HolderBulk+HolderY+RailY],
		[-(HolderX), HolderBulk+ServoY],
		[-(HolderX+HolderBulk), HolderBulk+ServoY],
		[-(HolderX+HolderBulk), HolderBulk],
		[-HolderX,0]
	
	]);
}

module A_Screw()
{
	translate([0,HolderBulk/2,-ScrewHoleDepth]) cylinder(h=ScrewHoleDepth,d=ScrewHoleDiameter,$fn=30);
}

module A_Piece()
{
	difference()
	{
		A_Profile();
		A_Screw();
	}
}

//A_Piece();

//////////////////////////////////////////////
// Piece B
//////////////////////////////////////////////
Platform         = 10;
MountX           = RailX+Platform;
ScrewGuard       = 3;
ScrewHeadDiameter=3;

module B_Profile()
{
	rotate([0,-90,0])
	linear_extrude(height=HolderSpace,center=true)
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
	translate([0,HolderBulk/2,0])
	cylinder(h=ScrewGuard,d=ScrewHoleDiameter+0.5,$fn=30);
	translate([0,HolderBulk/2,ScrewGuard])
	cylinder(h=MountX-ScrewGuard,d=ScrewHeadDiameter,$fn=30);
}	

// To Attach servo
ServoScrewDiameter=1;

module B_Screw2()
{
	translate([0,HolderBulk+HolderY,2])
	rotate([90,0,0]) cylinder(h=5,d=ServoScrewDiameter,$fn=30);
}

// To Attach BackBone
AxisDiameter=1;

module B_Screw3()
{
	translate([0,HolderBulk+HolderY+RailY,RailX+Platform/2])
	rotate([90,0,0]) 
	cylinder(h=5,d=AxisDiameter,$fn=30);
}

module B_Piece()
{
	difference()
	{
		B_Profile();
		B_Screw1();
		B_Screw2();
		B_Screw3();
	}
}

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
// trial Servo
//////////////////////////////////////////////
module Servo()
{
	translate([0,HolderBulk+1,-12])
	 rotate([0,90,0]) 
	scale(10) color("green") import("data/microservo.stl");
}


//////////////////////////////////////////////
// BackBone FingerPrint
//////////////////////////////////////////////

BackBoneDiameter=18+0.5;
module BackBoneAttach()
{
	color("yellow")
	translate([0,HolderBulk+HolderY+RailY+BackBoneDiameter/2,RailX])
	cylinder(h=Platform,d=BackBoneDiameter,$fn=30);
}

//////////////////////////////////////////////
// Servo Arm
//////////////////////////////////////////////

ServoArmDiameter =25;
ServoArmThickness=5;
ServoArmCenter   =16;

module ServoArm()
{
	translate([0,ServoY+HolderBulk,-ServoArmCenter])
	rotate([-90,0,0])
	cylinder(h=ServoArmThickness,d=ServoArmDiameter,$fn=30);
}

module ServoSpace()
{
	hull()
	{
		ServoArm();
		translate([0,0,ServoArmCenter+RailX+Platform/2])
		ServoArm();
	}
}


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

//BigHole();
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

Toto(); 
/*
translate(v) 
{ 
	Servo(); 
	translate([0,0,50]) B_Piece();
}
*/


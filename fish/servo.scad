
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


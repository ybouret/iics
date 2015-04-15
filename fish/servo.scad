
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
ServoArmThickness=8;
ServoArmCenter   =16;

module ServoArm()
{
	translate([0,ServoY+HolderBulk-ServoArmThickness/2,-ServoArmCenter])
	rotate([-90,0,0])
	cylinder(h=ServoArmThickness,d=ServoArmDiameter,$fn=Resolution);
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


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
	
//FishTail();
//FishHead();
//translate([10,0,0]) rotate([0,90,0]) translate([0,0,20]) FishTail();

// real length Holder
module Holder()
{

	DX=10;
	DY=20;
	DZ=50;
	//scale([DX,DY,DZ])	cube(height=1,center=true);
	rotate([0,90,0]) linear_extrude(height=DX,center=true) 
	{
		difference()	
		{
			circle(r=DZ/2,center=true,$fn=30);
			translate([0,DZ/2]) square(DZ,center=true);
		}
	}
}

//%scale([3,3,3]) FishHead();
//Holder();
scale(10) import("data/microservo.stl");
//translate([0,8,0]) cube([24,16,14],center=true);




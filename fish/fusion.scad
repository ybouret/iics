module fish()
{
%import("fish_head.stl");
//%import("fish_tail.stl");
}

//scale([2,2,2]) fish();

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
FishHead();
translate([10,0,0]) rotate([0,90,0]) translate([0,0,20]) FishTail();

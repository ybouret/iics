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
	
//scale([2,2,2])  import("fish_head.stl");
scale([2,2,2])  import("fish_tail.stl");

//Holder();

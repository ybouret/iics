Resolution = 30;

//____________________________________________________________________
//
// unscaled Fish Junction
//____________________________________________________________________
module FishJunction()
{
	render()
	color("red")
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

module FishTail(zoom=1)
{
	scale(zoom)
	import("fish_tail.stl");
}

FishHead(1.42);
%FishTail(1.42);
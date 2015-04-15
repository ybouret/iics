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
	translate([0,HolderBulk/2,-ScrewHoleDepth]) cylinder(h=ScrewHoleDepth,d=ScrewHoleDiameter,$fn=Resolution);
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

// num links
n=20;

// side
s=10;

R=s/(2*sin(180/n));


module one()
{
    translate([s/2,0,0])
    import("maillon.stl");
}


//%one();

module maillon(j)
{
    angle = (j*360)/n;
    xpos  = R*cos(angle);
    ypos  = R*sin(angle);
    //color("red") translate([xpos,ypos,0])  cube([1,1,30],center=true);
    translate([xpos,ypos,0])  rotate([0,0,90+angle+180/n]) one();
}

for(i = [0:n-1])
{
    maillon(i);
}
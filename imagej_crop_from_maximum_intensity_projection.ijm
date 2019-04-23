//Sets input and output
input = "/.../";
output = "/.../"

//Loop
list = getFileList(input);
for (i=0; i < list.length; i++) {
	if (startsWith(list[i], "hyb4"))
	if (endsWith(list[i], "39_MIP.tif"))
	action(input, output, list[i]);

function action(input, output, filename) {
	
open(input + filename); {

t=File.nameWithoutExtension;	


// Crop and save. For this, need to use 'Record' function built in on imagej, crop cell you want manually and obtain coordinates from the 'record' window.

setTool("polygon");
makePolygon(//coordinates//);
run("Duplicate...", "duplicate");
saveAs("Tiff", output+t+"//name//.tif");
close();

close();
}
}
}
